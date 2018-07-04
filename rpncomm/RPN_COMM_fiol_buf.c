/* RPN_COMM - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 2018  Division de Recherche en Prevision Numerique
 *                          Environnement Canada
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include <stdint.h>
#include <unistd.h>
#include <stdlib.h>

// circular buffer management functions

/*
 *     circular data buffer layout
 * 
 *     first : points to first element of buffer            (never updated)
 *     limit : points one past the last element in buffer   (never updated)
 *     in    : points to insertion position                 (only updated by thread inserting data)
 *     out   : points to extraction position                (only updated by thread extracting data)
 * 
 *     this will then work without having to use locks when updating in and out
 *     (provided only ONE player inserts and ONE player extracts)
 *     buffer may contain between 0 and limit-first-1 items
 * 
 *     empty buffer, in = out
 * 
 *     +-----------------------------------------------------------------------------------------+
 *     | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | +
 *     +-----------------------------------------------------------------------------------------+
 *      ^                   ^                                                                     ^
 *      | first          in | out                                                                 | limit
 *     ( 0 data elements in buffer, room for (limit - first - 1)
 * 
 *     buffer neither empty nor full
 * 
 *     +-----------------------------------------------------------------------------------------+
 *     | | | | | | | | | | |d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d| | | | | | | | | | | | | | | | +
 *     +-----------------------------------------------------------------------------------------+
 *      ^                   ^                                     ^                               ^
 *      | first             | out                                 | in                            | limit
 *     (in - out) data elements in buffer, room for (limit - in) + (out - in -1) more
 * 
 * 
 *     +-----------------------------------------------------------------------------------------+
 *     |d|d|d|d|d|d|d|d|d|d| | | | | | | | | | | | | | | | | | | |d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d+
 *     +-----------------------------------------------------------------------------------------+
 *      ^                   ^                                     ^                               ^
 *      | first             | in                                  | out                           | limit
 *     (in - first) + limit - out) data elements in buffer, room for (out - in + 1)
 * 
 * 
 *     full buffer, in = (out -1) (modulo limit)
 *     +-----------------------------------------------------------------------------------------+
 *     |d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d| |d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d+
 *     +-----------------------------------------------------------------------------------------+
 *      ^                                                       ^ ^                               ^
 *      | first                                              in | | out                           | limit
 *     +-----------------------------------------------------------------------------------------+
 *     |d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d| |d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d+
 *     +-----------------------------------------------------------------------------------------+
 *      ^                                                                                       ^ ^
 *      | first,out                                                                          in | | limit
 */


// this struct contains the circular buffer indexes
typedef struct{   // first, in, out, limit are in "origin 0"
  uint32_t first;  // index in data of first point of circular buffer
  uint32_t in;     // index in data of insertion point into circular buffer
  uint32_t out;    // index in data of extraction point into circular buffer
  uint32_t limit;  // index in data of point after last point of circular buffer
  uint32_t data[1];
} fiol ;

typedef struct{
  fiol *inbound;
  fiol *outbound;
}fioltabentry;

#if defined(SELF_TEST)
#include <mpi.h>
#endif

#if defined(SELF_TEST) || defined(DIAG)

#include <stdio.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <sys/stat.h>

static int64_t stalled_insert  = 0;   // nb of times waiting on buffer full
static int64_t block_insert    = 0;   // number of block insertions
static int64_t stalled_extract = 0;   // nb of times waiting on buffer empty
static int64_t block_extract   = 0;   // number of block extractions

#endif

static fioltabentry *fioltable = NULL;
static fioltabentry my_fiol = {NULL, NULL};
static uint32_t my_buffer_id = -1;
static int32_t max_id = -1;

// initialize fiol struct
// b        : pointer to control structure for this circular buffer
// datasize : number of 4 byte items in this circular buffer
void FiolBufInit(fiol *b, int32_t datasize){
  b->first = 0;
  b->in    = 0;
  b->out   = 0;
  b->limit = datasize;
}

// create and initialize a table for npairs channel pairs
// using size Bytes from base address mem (addresses mem -> mem+size-1)
// i am thread/process my_id ( 0 <= my_id < npairs)
// return 0 if unsuccessful, npairs if successful
// my_id == 0 identifies the "master" thread/process
int32_t FiolTableInit(unsigned char *mem, uint32_t size, uint32_t npairs, uint32_t my_id){
  int32_t perchannel = size / (2 * npairs);   // number of bytes available per channel (2 channels per "player")
  int32_t ndata;
  int i;
  unsigned char *base = mem;
  fiol *fi, *fo;

  if(perchannel < sizeof(fiol) + 15*sizeof(uint32_t) ) return 0;       // less than 16 data items per channel
  if(my_id < 0 || my_id >= npairs) return 0;                           // invalid id 
  if(fioltable != NULL) return 0;                                      // table already initialized

  fioltable = (fioltabentry *) malloc (2 * npairs * sizeof(fiol *));     // allocate table with 2 * npairs entries
  if(fioltable == NULL) return 0;                               // table allocation failed

  ndata = 1 + (perchannel - sizeof(fiol)) / sizeof(uint32_t);   // size of data array (per channel)
  for(i=0 ; i<npairs ; i++) {
    fioltable[i].inbound  = (fiol *) base;                      // base address for inbound channel
    if(my_id == 0) FiolBufInit(fioltable[i].inbound,ndata);     // initialize inbound channel only on "master"
    base = base + perchannel;
    fioltable[i].outbound = (fiol *) base;                      // base address for outbound channel
    if(my_id == 0) FiolBufInit(fioltable[i].outbound,ndata);    // initialize outbound channel  only on "master"
    base = base + perchannel;
  }
  my_fiol.inbound  = fioltable[my_id].inbound;                  // my inbound circular buffer
  my_fiol.outbound = fioltable[my_id].outbound;                 // my outbound circular buffer
  my_buffer_id = my_id;
  max_id = npairs - 1;

#if defined(DIAG)
  printf("     PE  frst    in   out limit      address      frst    in   out limit      address  (%d)\n",max_id);
  for(i=0 ; i<npairs; i++){
    fi = fioltable[i].inbound ; fo = fioltable[i].outbound;
    printf("%c %5d %5d %5d %5d %5d %16p %5d %5d %5d %5d %16p \n",(my_id == i) ? '*' : ' ', i,
           fi->first, fi->in, fi->out, fi->limit, &(fi->data[0]),
           fo->first, fo->in, fo->out, fo->limit, &(fo->data[0]));
  }
#endif

  return npairs;    // return number of channel pairs allocated
}

// newer version, insertion by blocks as large as possible
// target   : -1 insert into MY outbound channel
//            tid insert into inbound channel of thread/process with id == tid
//            (-1 <= target < max_id)
// b        : pointer to control structure for targeted circular buffer
// userdata : data from user
// nw       : number of 4 byte items from user
// return   : number of items succesfully inserted into buffer
// Blocking version
int32_t FiolBufInsert(int32_t target, void *userdata, int32_t nw){
  fiol volatile *b;
  uint32_t status = 0;        // number of items processed
  uint32_t in ;               // in is not volatile for this routine as it is the only one updating it
  uint32_t volatile *out;
  uint32_t limit ;            // never changes, no point in keeping it volatile via *b
  uint32_t first ;            // never changes, no point in keeping it volatile via *b
  uint32_t temp, temp2;
  int32_t caninsert;          // how many values we can insert in one pass 

#if defined(SELF_TEST) || defined(DIAG)
  uint32_t timeout = 1000000;
  printf("mass insert of %d items\n",nw);
#endif
  if(target > max_id) return 0;
  b = (target < 0) ?  fioltable[my_buffer_id].outbound : fioltable[target].inbound  ;
  in    = b->in;
  out = &(b->out);
  limit = b->limit;
  first = b->first;
  while(nw > 0){
    temp = *out;                        // capture out
    temp2 = (temp == first) ? 1 : 0;    // if out = first, one may not insert beyond b->data[limit-2]
    caninsert = (temp > in) ? (temp - in - 1) : (limit - in - temp2) ;
    caninsert = (caninsert > nw) ? nw : caninsert ;   // no more than nw anyway
    if(caninsert > 0){
      nw = nw - caninsert ;
      while(caninsert--){               // insert caninsert items
	b->data[in++] = ((uint32_t *)userdata)[status++];
      }
      in = (in < limit) ? in : first ;  // handle possible wraparound at limit
      b->in = in;                       // end of block insert, update in pointer
#if defined(DIAG)
      block_insert++;
#endif
    }else{
#if defined(DIAG)
      stalled_insert++;
      timeout--;
      if(timeout <= 0){ 
	printf("I1 first, in, out, limit = %d %d %d %d\n\n",first, b->in, b->out, limit);
	return status; 
      }
#endif
    }
  }
  return status;     // return number of items inserted (should be equal to nw)
}

// older version, one by one insertion into circular buffer
int32_t FiolBufInsertBy1(int32_t target, void *userdata, int32_t nw){
  fiol volatile *b;
  int status = 0;
  int nextin;
  int stalled = 0;
  uint32_t limit;
  uint32_t first;

  if(target > max_id) return 0;
  b = (target < 0) ?  fioltable[my_buffer_id].outbound : fioltable[target].inbound  ;
  nextin = b->in;
  limit = b->limit;
  first = b->first;
  while(nw--){                                                    // not all data inserted
    nextin = ((nextin + 1) < limit) ? nextin + 1 : first;         // next insertion point
    while(nextin == b->out){                                      // circular buffer is full
#if defined(DIAG)
      stalled_insert++;
#endif
    }
    b->data[b->in] = ((int32_t *)userdata)[status]; // insert data
    status = status + 1;                                          // one more token done
    b->in = nextin;                                               // update in after successful insertion
  }
  return status;     // return number of items inserted
}

// Blocking insert of a single item
void FiolBufInsert1(int32_t target, int32_t userdata){
  fiol volatile *b;
  int nextin;

printf("FiolBufInsert1 target = %d\n",target);
  if(target > max_id) return;
  b = (target < 0) ?  fioltable[my_buffer_id].outbound : fioltable[target].inbound  ;
  nextin = b->in ;
  nextin = ((nextin + 1) < b->limit) ? nextin + 1 : b->first;   // next insertion point, possible wraparound
  while(nextin == b->out){                                      // circular buffer is full
#if defined(DIAG)
    stalled_insert++;
#endif
  }
  b->data[b->in] = userdata ;                                   // insert data
  b->in = nextin;                                               // update in after successful insertion
}

// Non Blocking version of FiolBufInsertBy1
// 
int32_t FiolBufInsertNB(int32_t target, void *userdata, int32_t nw){
  fiol volatile *b;
  int status = 0;
  int nextin ;

  if(target > max_id) return 0;
  b = (target < 0) ?  fioltable[my_buffer_id].outbound : fioltable[target].inbound  ;
  nextin = b->in ;
  while(nw--){                                                    // not all data inserted
    nextin = ((nextin + 1) < b->limit) ? nextin + 1 : b->first;   // next insertion point
    if(nextin == b->out)return status;                            // circular buffer is full
    b->data[b->in] = ((int32_t *)userdata)[status];               // insert data
    status = status + 1;                                          // one more token done
    b->in = nextin;                                               // update in after successful insertion
  }
  return status;                                                  // return number of items inserted
}

// extract items from circular buffer by blocks as large as possible
// target   : -1 extract from MY inbound channel
//            tid extract from outbound channel of thread/process with id == tid
//            ( -1 <= target < max_id )
// b        : pointer to control structure for targeted circular buffer
// userdata : user buffer that will receive data
// nw       : number of 4 byte items to user
// return   : number of items succesfully extracted from buffer
// Blocking version
uint32_t FiolBufExtract(int32_t target, void *userdata, int32_t nw){
  volatile fiol *b;
  uint32_t status = 0;
  uint32_t out ;              // not volatile, as this function is the only one updating it
  uint32_t first ;            // never changes, no point in keeping it volatile via *b
  uint32_t limit ;            // never changes, no point in keeping it volatile via *b
  uint32_t canextract;        // number of items that can be extracted in one pass
  uint32_t temp;

#if defined(SELF_TEST) || defined(DIAG)
  uint32_t timeout = 10000000;
  printf("mass extract of %d items\n",nw);
#endif
  if(target > max_id) return 0;
  b = (target < 0) ?  fioltable[my_buffer_id].inbound : fioltable[target].outbound  ;
  out   = b->out; 
  first = b->first;
  limit = b->limit;
  while(nw > 0){                              // not all data inserted
    while(b->in == out){                      // circular buffer is empty
#if defined(DIAG)
      stalled_extract++;
      timeout--;
      if(timeout <= 0){ 
	printf("EB1 first, in, out, limit = %d %d %d %d\n\n",first, b->in, b->out, limit);
	return status; 
      }
#endif
    }
    temp = b->in;
    canextract = (out > temp) ? (limit - out) : (temp - out) ;
    canextract = (canextract > nw) ? nw : canextract;          // max nw items
    if(canextract > 0){
      nw = nw - canextract;                                    // update item count
      while(canextract--){
	((int32_t *)userdata)[status++] = b->data[out++];      // extract canextract items
      }
      out = (out < limit) ? out : first;
      b->out = out;
#if defined(DIAG)
      block_extract++;
    }else{
      stalled_extract++;
      timeout--;
      if(timeout < 0){ 
	printf("EB2 first, in, out, limit = %d %d %d %d\n\n",first, b->in, b->out, limit);
	return status; 
      }
#endif
    }
  }
  return status;   // number of items extracted
}

// older version, extract items one by one
int32_t FiolBufExtractBy1(fiol volatile *b, void *userdata, int32_t nw){
  int32_t target;
  int status = 0;
  uint32_t first;  // never changes, no point in keeping it volatile via *b
  uint32_t limit;  // never changes, no point in keeping it volatile via *b

  if(target > max_id) return 0;
  b = (target < 0) ?  fioltable[my_buffer_id].inbound : fioltable[target].outbound  ;
  first = b->first;
  limit = b->limit;
  while(nw--){                                                    // not all data inserted
    while(b->in == b->out){                                       // circular buffer is empty
#if defined(DIAG)
      stalled_extract++;
#endif
    }
    ((int32_t *)userdata)[status] = b->data[b->out]; // extract data
    status = status + 1;                                          // one more token done
    b->out = (b->out+1 < limit) ? b->out+1 : first;         // next extraction point
  }
  return status;
}

// Blocking extract of single item (item = 4 Bytes)
// returns extracted item
int32_t FiolBufExtract1(int32_t target){
  fiol volatile *b;
  int temp = 0;

  if(target > max_id) return 0;
  b = (target < 0) ?  fioltable[my_buffer_id].inbound : fioltable[target].outbound  ;
  while(b->in == b->out){                                       // circular buffer is empty
#if defined(DIAG)
    stalled_extract++;
#endif
  }
  temp = b->data[b->out];                                     // extract datum
  b->out = (b->out+1 < b->limit) ? b->out+1 : b->first;         // next extraction point
  return temp;  // return extracted item
}

// Non Blocking version of FiolBufExtractBy1
int32_t FiolBufExtractNB(int32_t target, void *userdata, int32_t nw){
  fiol volatile *b;
  int status = 0;

  if(target > max_id) return 0;
  b = (target < 0) ?  fioltable[my_buffer_id].inbound : fioltable[target].outbound  ;
  while(nw--){                                                    // not all data inserted
    if(b->in == b->out) return status;                            // circular buffer is empty
    ((int32_t *)userdata)[status] = b->data[b->out];              // extract data
    status = status + 1;                                          // one more token done
    b->out = (b->out+1 < b->limit) ? b->out+1 : b->first;         // next extraction point
  }
  return status;   // number of items extracted
}

// b        : pointer to control structure for this circular buffer
// if buffer is empty, return number_of_available_slots (buffer capacity)
// if buffer is full, return 0
// otherwise return -(number_of_free_slots)
int32_t FiolBufStatus(int32_t target, int32_t dir){
  fiol *b;
  int status ;

  if(target > max_id) return -999999999;
  target = (target < 0) ?  my_buffer_id : target ;
  b = (dir == 0) ? fioltable[target].inbound : fioltable[target].outbound ;
//   printf("FiolBufStatus(%d): first, in, out, limit = %d %d %d %d\n",target,b->first, b->in, b->out, b->limit);
  if(b->in == b->out) {
    status = b->limit - b->first - 1 ;                //  buffer is empty, return buffer limit - 1
  }else{
    if(b->in > b->out) {                              // in - out tokens in buffer
      status = -(b->limit - b->first -1 - (b->in - b->out)); // number of free slots
    }else{                                            // out - in - 1 free slots
      status = -(b->out - b->in - 1) ;                // number of free slots
    }
  }
  return status;
}

#if defined(SELF_TEST)

static uint64_t rdtsc(void) {   // version rapide "out of order"
  uint32_t lo, hi;
  __asm__ volatile ("rdtsc"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
  return (uint64_t)lo | (((uint64_t)hi) << 32);
}

#define BUFSZ 1024*64
#define NW 15
#define NREP 10

uint32_t verify(uint32_t *v, uint32_t nv){
  int i;
  uint32_t e=0;
  for(i=0 ; i<nv ; i++) {
    if(v[i] != i) e++;
    if(e == 1) printf("first error at position %d\n",i);
  }
  return e;
}

main(int argc, char **argv){
  int32_t buf[NW*1024*1024];
  int32_t status, status0;
  int ierr, id;
  int localrank, localsize;
  int globalrank, globalsize;
  int peerrank, peersize;
  fiol *ptr = NULL;
  fiol *tst = NULL;
  struct shmid_ds shm_buf;
  size_t size;
  int sizei;
  int i, nrep;
  uint64_t t0 , t1, t2;
  unsigned char *tmp;
  MPI_Comm MY_World = MPI_COMM_NULL;
  MPI_Comm MY_Peers = MPI_COMM_NULL;
  MPI_Comm temp_comm;
  int myhost, myhost0, myhost1;

  ierr = MPI_Init( &argc, &argv );
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&globalrank);
  ierr = MPI_Comm_size(MPI_COMM_WORLD,&globalsize);

  myhost  = gethostid();
  myhost0 = myhost & 0x7FFFFFFF  ; // lower 31 bits
  myhost1 = (myhost >> 31) & 0x1 ; // upper bit
  ierr = MPI_Comm_split(MPI_COMM_WORLD,myhost0,globalrank,&temp_comm) ;  // split WORLD using the lower 31 bits of host id , weight=rank in base
  ierr = MPI_Comm_split(temp_comm     ,myhost1,globalrank,&MY_World) ;   // re split using the upper bit of host id , weight=rank in base
  ierr = MPI_Comm_rank(MY_World,&localrank);                             // rank of this PE on this SMP node
  ierr = MPI_Comm_size(MY_World,&localsize);                             // number of PEs on this SMP node
  ierr = MPI_Comm_split(MPI_COMM_WORLD,localrank,globalrank,&MY_Peers) ; // communicator with PES of same local rank on other SMP nodes
  ierr = MPI_Comm_rank(MY_Peers,&peerrank);
  ierr = MPI_Comm_size(MY_Peers,&peersize);
  printf("number of peers in MY_Peers = %d\n",peersize);

  size = (BUFSZ +4 ) * localsize * 2 * 4 ;
  sizei = size;

  if(localrank == 0){      // node SERVER PE
    id = shmget(IPC_PRIVATE,size,IPC_CREAT|S_IRUSR|S_IWUSR);  /* rank 0 allocates shared memory segment */
    ptr = shmat(id,NULL,0);

    status = FiolTableInit((unsigned char *)ptr, sizei, localsize, localrank);  // initialize as master
    tst = fioltable[1].inbound;
    // for this test we are using channel[0] inbound

    shmctl(id,IPC_RMID,&shm_buf);      /* mark segment for deletion preventively (only works on linux) */
    for(i=0 ; i<BUFSZ ; i++) ptr->data[i] = 0;  // force pages to be properly instanciated in memory

    for(i=0 ; i<sizeof(buf)/sizeof(int32_t) ; i++) buf[i] = i;

    t0 = rdtsc();
    for(i=0 ; i<NW ; i++) FiolBufInsert1(1, buf[i]);
    t1 = rdtsc();
    printf("single insert of %d words, time = %ld, per item = %ld\n",NW,t1-t0,(t1-t0)/NW);
    status = FiolBufStatus(1, 0) ;
    printf("0 buffer status is %d, expecting %d\n",status,-(BUFSZ-NW-1));
    printf("0 first, in, out, limit = %d %d %d %d\n\n",tst->first, tst->in, tst->out, tst->limit);

    ierr = MPI_Bcast(&id,1,MPI_INTEGER,0,MY_World);
    printf("0 shared memory segment ID = %d, ptr = %p\n\n",id,ptr);
    status = FiolBufStatus(1, 0) ;
    printf("0 shared buffer status is %d, expecting %d\n",status,-(BUFSZ-NW-1));
    printf("0 first, in, out, limit = %d %d %d %d\n\n",tst->first, tst->in, tst->out, tst->limit);

    ierr = MPI_Barrier(MY_World);

    t0 = rdtsc();
    status = FiolBufInsert(1, buf, NW*4);
    t1 = rdtsc();
    printf("0 status = %d, expected %d, time = %ld, per item = %ld\n",status,NW*4,t1-t0,(t1-t0)/(NW*4));
    printf("0 first, in, out, limit = %d %d %d %d\n\n",tst->first, tst->in, tst->out, tst->limit);
    status = 0;
    status0 = FiolBufInsert(1, buf, NW*1000);
    t0 = rdtsc();
    for(i=0 ; i<NREP ; i++) status += FiolBufInsert(1, buf, NW*1000);
    t1 = rdtsc();
    printf("0 status = %d, expected %d, time = %ld, per item = %ld\n",
	   status,NREP*NW*1000,t1-t0,(t1-t0)/(NREP*NW*1000));
    printf("0 first, in, out, limit = %d %d %d %d\n\n",tst->first, tst->in, tst->out, tst->limit);

    ierr = MPI_Barrier(MY_World);
    status = FiolBufStatus(-1, 0) ;
    printf("0 shared buffer status is %d, expecting %d\n",status,BUFSZ-1);
    printf("0 first, in, out, limit = %d %d %d %d\n\n",tst->first, tst->in, tst->out, tst->limit);
    printf("0 buf[0] = %d, buf[last] = %d, buf[last+1] = %d\n",buf[0],buf[NW*1000-1],buf[NW*1000]);
    printf("0 stalled_insert = %ld, stalled_extract = %ld, block_insert = %ld, block_extract = %ld\n",
	   stalled_insert,stalled_extract,block_insert,block_extract);

    // local/global "all reduce" test
    ierr = MPI_Barrier(MPI_COMM_WORLD);
    t0 = rdtsc();
    ierr = MPI_Barrier(MY_World);
    t2 = rdtsc()-t0;
    t0 = rdtsc();
    FiolBufInsert1(-1, localrank);                                   // insert into my own outbound buffer
    for(i=0 ; i<localsize ; i++) buf[i] = FiolBufExtract1(i);        // get from all PEs on node outbound buffer
    for(i=1 ; i<localsize ; i++) buf[0] = buf[0] + buf[i];           // local sum
    for(i=0 ; i<localsize ; i++) FiolBufInsert1(i,buf[0]);           // broadcast sum to all inbound buffers
    status = FiolBufExtract1(-1);                                    // get answer from my own inbound buffer
    t1 = rdtsc();
    printf("%5d answer1 = %d, time = %ld, barrier = %ld\n",localrank,status,t1-t0,t2);

    ierr = MPI_Barrier(MPI_COMM_WORLD);
    t0 = rdtsc();
    FiolBufInsert1(-1, localrank);                                   // insert into my own outbound buffer
    for(i=0 ; i<localsize ; i++) buf[i] = FiolBufExtract1(i);        // get from all PEs on node outbound buffer
    for(i=1 ; i<localsize ; i++) buf[0] = buf[0] + buf[i];           // local sum
    ierr = MPI_Allreduce(&buf[0],&buf[1],1,MPI_INTEGER,MPI_SUM,MY_Peers);
    buf[0] = buf[1];
    for(i=0 ; i<localsize ; i++) FiolBufInsert1(i,buf[0]);           // broadcast sum to all inbound buffers
    status = FiolBufExtract1(-1);                                    // get answer from my own inbound buffer
    t1 = rdtsc();
    printf("%5d answer2 = %d, time = %ld\n",localrank,status,t1-t0);
  }else if(localrank == 1){             // node CLIENT PE 1
    ierr = MPI_Bcast(&id,1,MPI_INTEGER,0,MY_World);
    ptr = shmat(id,NULL,0);

    status = FiolTableInit((unsigned char *)ptr, sizei, localsize, localrank);
    tst = fioltable[1].inbound;

    printf("%d shared memory segment ID = %d, ptr = %p\n\n",localrank,id,ptr);
    status = FiolBufStatus(-1, 0) ;
    printf("1 shared buffer status is %d, expecting %d\n",status,-(BUFSZ-NW-1));
    printf("1 first, in, out, limit = %d %d %d %d\n\n",tst->first, tst->in, tst->out, tst->limit);

    for(i=0 ; i<sizeof(buf)/sizeof(int32_t) ; i++) buf[i] = 999999;
    ierr = MPI_Barrier(MY_World);

    t0 = rdtsc();
    for(i=0 ; i<NW ; i++) buf[i] = FiolBufExtract1(-1);
    t1 = rdtsc();
    printf("1 single extract of %d words, time = %ld, per item = %ld\n",NW,t1-t0,(t1-t0)/NW);
    printf("time = %ld, per item = %ld\n",t1-t0,(t1-t0)/NW);
    printf("1 buf[0] = %d, buf[1] = %d, buf[last] = %d, buf[last+1] = %d\n",buf[0],buf[1],buf[NW-1],buf[NW]);

    t0 = rdtsc();
    status = FiolBufExtract(-1, buf, NW*4);
    t1 = rdtsc();
    printf("1 status = %d, expected %d, time = %ld, per item = %ld\n",status,NW*4,t1-t0,(t1-t0)/(NW*4));
    printf("1 buf[0] = %d, buf[1] = %d, buf[last] = %d, buf[last+1] = %d\n",buf[0],buf[1],buf[NW*4-1],buf[NW*4]);

    status = 0;
    status0 = FiolBufExtract(-1, buf, NW*1000);
    t0 = rdtsc();
    for(i=0 ; i<NREP ; i++) status += FiolBufExtract(-1, buf, NW*1000);
    t1 = rdtsc();
    ierr = verify(buf, status0);
    printf("1 status = %d, expected %d, time = %ld, errors = %d, per item = %ld\n",
	   status,NREP*NW*1000,t1-t0,ierr,(t1-t0)/(NREP*NW*1000));
    printf("1 buf[0] = %d, buf[1] = %d, buf[last] = %d, buf[last+1] = %d\n",buf[0],buf[1],buf[status0-1],buf[status0]);

    ierr = MPI_Barrier(MY_World);
    status = FiolBufStatus(-1, 0) ;
    printf("1 shared buffer status is %d, expecting %d\n",status,BUFSZ-1);
    printf("1 first, in, out, limit = %d %d %d %d\n\n",tst->first, tst->in, tst->out, tst->limit);
    printf("1 stalled_insert = %ld, stalled_extract = %ld, block_insert = %ld, block_extract = %ld\n",
	   stalled_insert,stalled_extract,block_insert,block_extract);
    // local/global "all reduce" test
    ierr = MPI_Barrier(MPI_COMM_WORLD);
    t0 = rdtsc();
    ierr = MPI_Barrier(MY_World);
    t2 = rdtsc()-t0;
    t0 = rdtsc();
    FiolBufInsert1(-1, localrank);                                   // insert into my own outbound buffer
    status = FiolBufExtract1(-1);                                    // get answer from my own inbound buffer
    t1 = rdtsc();
    printf("%5d answer1 = %d, time = %ld, barrier = %ld\n",localrank,status,t1-t0,t2);

    ierr = MPI_Barrier(MPI_COMM_WORLD);
    t0 = rdtsc();
    FiolBufInsert1(-1, localrank);                                   // insert into my own outbound buffer
    status = FiolBufExtract1(-1);                                    // get answer from my own inbound buffer
    t1 = rdtsc();
    printf("%5d answer2 = %d, time = %ld\n",localrank,status,t1-t0);
  }else{                    // node CLIENT PEs 2 -> localsize - 1
    ierr = MPI_Bcast(&id,1,MPI_INTEGER,0,MY_World);
    ptr = shmat(id,NULL,0);

    status = FiolTableInit((unsigned char *)ptr, sizei, localsize, localrank);

    ierr = MPI_Barrier(MY_World);  // all other ranks participate in barriers
    ierr = MPI_Barrier(MY_World);
    // local/global "all reduce" test
    ierr = MPI_Barrier(MPI_COMM_WORLD);
    t0 = rdtsc();
    ierr = MPI_Barrier(MY_World);
    t2 = rdtsc()-t0;
    t0 = rdtsc();
    FiolBufInsert1(-1, localrank);                                   // insert into my own outbound buffer
    status = FiolBufExtract1(-1);                                    // get answer from my own inbound buffer
    t1 = rdtsc();
    printf("%5d answer1 = %d, time = %ld, barrier = %ld\n",localrank,status,t1-t0,t2);

    ierr = MPI_Barrier(MPI_COMM_WORLD);
    t0 = rdtsc();
    FiolBufInsert1(-1, localrank);                                   // insert into my own outbound buffer
    status = FiolBufExtract1(-1);                                    // get answer from my own inbound buffer
    t1 = rdtsc();
    printf("%5d answer2 = %d, time = %ld\n",localrank,status,t1-t0);
  }

  ierr = shmdt(ptr);
  ierr = MPI_Finalize();
}
#endif
