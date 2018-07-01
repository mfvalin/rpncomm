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

// circular buffer management functions

/*
 *     circular data buffer layout
 * 
 *     first : points to first element of buffer            (never updated)
 *     limit : points one past the last element in buffer   (never updated)
 *     in    : points to insertion position                 (only updated by thread inserting data)
 *     out   : points to extraction position                (only updated by thread extracting data)
 * 
 *     this can then work without having to use locks when updating in and out
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
 */


// this struct contains the circular buffer indexes
typedef struct{   // first, in, out, limit are in "origin 0"
  uint32_t first;  // index in data of first point of circular buffer
  uint32_t in;     // index in data of insertion point into circular buffer
  uint32_t out;    // index in data of extraction point into circular buffer
  uint32_t limit;  // index in data of point after last point of circular buffer
  uint32_t data[1];
} fiol ;

#if defined(SELF_TEST)
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <sys/stat.h>

static int64_t stalled_insert  = 0;   // nb of times waiting on buffer full
static int64_t block_insert    = 0;   // number of block insertions
static int64_t stalled_extract = 0;   // nb of times waiting on buffer empty
static int64_t block_extract   = 0;   // number of block extractions

#endif

// initialize fiol struct
// b        : pointer to control structure for this circular buffer
// datasize : number of 4 byte items in this circular buffer
void FiolBufInit(fiol *b, int32_t datasize){
  b->first = 0;
  b->in    = 0;
  b->out   = 0;
  b->limit = datasize;
}

uint32_t RPN_COMM_fetch(uint32_t *p);
void RPN_COMM_store(uint32_t *p, uint32_t v);

#define FiolBufInsertNew FiolBufInsert

// older version, one by one insertion
static int32_t FiolBufInsertBy1(fiol volatile *b, void *userdata, int32_t nw){
  int status = 0;
  int nextin = b->in ;
  int stalled = 0;
  uint32_t limit = b->limit;;
  uint32_t first = b->first;;

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

// newer version, insertion by blocks when possible
// b        : pointer to control structure for this circular buffer
// userdata : data from user
// nw       : number of 4 byte items from user
// return   : number of items succesfully inserted into buffer
// Blocking version
int32_t FiolBufInsert(fiol volatile *b, void *userdata, int32_t nw){
  uint32_t status = 0;        // number of items processed
  uint32_t in    = b->in;     // in is not volatile for this routine as it is the only one updating it
  uint32_t volatile *out = &(b->out);
  uint32_t limit = b->limit;  // never changes, no point in keeping it volatile via *b
  uint32_t first = b->first;  // never changes, no point in keeping it volatile via *b
  uint32_t temp, temp2;
  int32_t caninsert;     // how many values we can insert in one pass 

#if defined(SELF_TEST)
  uint32_t timeout = 1000000;
  printf("mass insert of %d items\n",nw);
#endif
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

// Blocking insert of a single item
void FiolBufInsert1(fiol volatile *b, int32_t userdata){
  int nextin = b->in ;

  nextin = ((nextin + 1) < b->limit) ? nextin + 1 : b->first;   // next insertion point, possible wraparound
  while(nextin == b->out){                                      // circular buffer is full
#if defined(DIAG)
    stalled_insert++;
#endif
  }
  b->data[b->in] = userdata ;                                   // insert data
  b->in = nextin;                                               // update in after successful insertion
}

// Non Blocking version of FiolBufInsert
// 
int32_t FiolBufInsertNB(fiol volatile *b, void *userdata, int32_t nw){
  int status = 0;
  int nextin = b->in ;

  while(nw--){                                                    // not all data inserted
    nextin = ((nextin + 1) < b->limit) ? nextin + 1 : b->first;   // next insertion point
    if(nextin == b->out)return status;                            // circular buffer is full
    b->data[b->in] = ((int32_t *)userdata)[status];               // insert data
    status = status + 1;                                          // one more token done
    b->in = nextin;                                               // update in after successful insertion
  }
  return status;                                                  // return number of items inserted
}

static int32_t FiolBufExtractOrig(fiol volatile *b, void *userdata, int32_t nw){
  int status = 0;

  while(nw--){                                                    // not all data inserted
    while(b->in == b->out){                                       // circular buffer is empty
#if defined(DIAG)
      stalled_extract++;
#endif
    }
    ((int32_t *)userdata)[status] = b->data[b->out]; // extract data
    status = status + 1;                                          // one more token done
    b->out = (b->out+1 < b->limit) ? b->out+1 : b->first;         // next extraction point
  }
  return status;
}

// b        : pointer to control structure for this circular buffer
// userdata : user buffer that will receive data
// nw       : number of 4 byte items to user
// return   : number of items succesfully extracted from buffer
// Blocking version
uint32_t FiolBufExtract(volatile fiol *b, void *userdata, int32_t nw){
  uint32_t status = 0;
  uint32_t out   = b->out;    // not volatile, as this function is the only one updating it
  uint32_t first = b->first;  // never changes, no point in keeping it volatile via *b
  uint32_t limit = b->limit;  // never changes, no point in keeping it volatile via *b
  uint32_t canextract;        // number of items that can be extracted in one pass
  uint32_t temp;

#if defined(SELF_TEST)
  uint32_t timeout = 1000000;
  printf("mass extract of %d items\n",nw);
#endif
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

// Blocking single datun version
// returns extracted datum
int32_t FiolBufExtract1(fiol volatile *b){
  int temp = 0;

    while(b->in == b->out){                                       // circular buffer is empty
#if defined(DIAG)
      stalled_extract++;
#endif
    }
    temp = b->data[b->out];                                     // extract datum
    b->out = (b->out+1 < b->limit) ? b->out+1 : b->first;         // next extraction point
  return temp;  // return extracted item
}

// Non Blocking version of FiolBufExtract
int32_t FiolBufExtractNB(fiol volatile *b, void *userdata, int32_t nw){
  int status = 0;

  while(nw--){                                                    // not all data inserted
    if(b->in == b->out) return status;                            // circular buffer is empty
    ((int32_t *)userdata)[status] = b->data[b->out];              // extract data
    status = status + 1;                                          // one more token done
    b->out = (b->out+1 < b->limit) ? b->out+1 : b->first;         // next extraction point
  }
  return status;   // number of items extracted
}

// b        : pointer to control structure for this circular buffer
// if buffer is empty, return +(number of available slots) (buffer capacity)
// if buffer is full, return 0
// otherwise return -(number of free slots)
int32_t FiolBufStatus(fiol *b){
  int status ;
//   printf("first, in, out, limit = %d %d %d %d\n",b->first, b->in, b->out, b->limit);
  if(b->in == b->out) {
    status = b->limit - b->first - 1 ;                //  buffer is empty, return buffer capacity - 1
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

#define BUFSZ 1024*2
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
  int ierr, localrank, id;
  fiol *ptr = NULL;
  struct shmid_ds shm_buf;
  size_t size;
  int i, nrep;
  uint64_t t0 , t1;

  ierr = MPI_Init( &argc, &argv );
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&localrank);

  if(localrank == 0){
    size = 1234567 ;
    id = shmget(IPC_PRIVATE,size,IPC_CREAT|S_IRUSR|S_IWUSR);  /* rank 0 allocates shared memory segment */
    ptr = shmat(id,NULL,0);
    shmctl(id,IPC_RMID,&shm_buf);      /* mark segment for deletion preventively (only works on linux) */
    for(i=0 ; i<BUFSZ ; i++) ptr->data[i] = 0;  // force pages to be properly instanciated in memory
    FiolBufInit(ptr, BUFSZ);

    for(i=0 ; i<sizeof(buf)/sizeof(int32_t) ; i++) buf[i] = i;
    t0 = rdtsc();
    for(i=0 ; i<NW ; i++) FiolBufInsert1(ptr, buf[i]);
    t1 = rdtsc();
    printf("single insert of %d words, time = %ld, per item = %ld\n",NW,t1-t0,(t1-t0)/NW);
    status = FiolBufStatus(ptr) ;
    printf("0 buffer status is %d, expecting %d\n",status,-(BUFSZ-NW-1));
    printf("0 first, in, out, limit = %d %d %d %d\n\n",ptr->first, ptr->in, ptr->out, ptr->limit);

    ierr = MPI_Bcast(&id,1,MPI_INTEGER,0,MPI_COMM_WORLD);
    printf("0 shared memory segment ID = %d, ptr = %p\n\n",id,ptr);
    status = FiolBufStatus(ptr) ;
    printf("0 shared buffer status is %d, expecting %d\n",status,-(BUFSZ-NW-1));
    printf("0 first, in, out, limit = %d %d %d %d\n\n",ptr->first, ptr->in, ptr->out, ptr->limit);

    ierr = MPI_Barrier(MPI_COMM_WORLD);

    t0 = rdtsc();
    status = FiolBufInsert(ptr, buf, NW*4);
    t1 = rdtsc();
    printf("0 status = %d, expected %d, time = %ld\n",status,NW*4,t1-t0);
    printf("0 first, in, out, limit = %d %d %d %d\n\n",ptr->first, ptr->in, ptr->out, ptr->limit);
    status = 0;
    status0 = FiolBufInsert(ptr, buf, NW*1000);
    t0 = rdtsc();
    for(i=0 ; i<NREP ; i++) status += FiolBufInsert(ptr, buf, NW*1000);
    t1 = rdtsc();
    printf("0 status = %d, expected %d, time = %ld, per item = %ld\n",
	   status,NREP*NW*1000,t1-t0,(t1-t0)/(NREP*NW*1000));
    printf("0 first, in, out, limit = %d %d %d %d\n\n",ptr->first, ptr->in, ptr->out, ptr->limit);

    ierr = MPI_Barrier(MPI_COMM_WORLD);
    status = FiolBufStatus(ptr) ;
    printf("0 shared buffer status is %d, expecting %d\n",status,BUFSZ-1);
    printf("0 first, in, out, limit = %d %d %d %d\n\n",ptr->first, ptr->in, ptr->out, ptr->limit);
    printf("0 buf[0] = %d, buf[last] = %d, buf[last+1] = %d\n",buf[0],buf[NW*1000-1],buf[NW*1000]);
    printf("0 stalled_insert = %ld, stalled_extract = %ld, block_insert = %ld, block_extract = %ld\n",
	   stalled_insert,stalled_extract,block_insert,block_extract);

  }else if(localrank == 1){
    ierr = MPI_Bcast(&id,1,MPI_INTEGER,0,MPI_COMM_WORLD);
    ptr = shmat(id,NULL,0);
    printf("%d shared memory segment ID = %d, ptr = %p\n\n",localrank,id,ptr);
    status = FiolBufStatus(ptr) ;
    printf("1 shared buffer status is %d, expecting %d\n",status,-(BUFSZ-NW-1));
    printf("1 first, in, out, limit = %d %d %d %d\n\n",ptr->first, ptr->in, ptr->out, ptr->limit);

    for(i=0 ; i<sizeof(buf)/sizeof(int32_t) ; i++) buf[i] = 999999;
    ierr = MPI_Barrier(MPI_COMM_WORLD);

    t0 = rdtsc();
    for(i=0 ; i<NW ; i++) buf[i] = FiolBufExtract1(ptr);
    t1 = rdtsc();
    printf("1 single extract of %d words, time = %ld, per item = %ld\n",NW,t1-t0,(t1-t0)/NW);
//     printf("status = %d, expected %d, time = %ld\n",status,NW,t1-t0);
    printf("1 buf[0] = %d, buf[1] = %d, buf[last] = %d, buf[last+1] = %d\n",buf[0],buf[1],buf[NW-1],buf[NW]);

    t0 = rdtsc();
    status = FiolBufExtract(ptr, buf, NW*4);
    t1 = rdtsc();
    printf("1 status = %d, expected %d, time = %ld\n",status,NW*4,t1-t0);
    printf("1 buf[0] = %d, buf[1] = %d, buf[last] = %d, buf[last+1] = %d\n",buf[0],buf[1],buf[NW*4-1],buf[NW*4]);

    status = 0;
    status0 = FiolBufExtract(ptr, buf, NW*1000);
    t0 = rdtsc();
    for(i=0 ; i<NREP ; i++) status += FiolBufExtract(ptr, buf, NW*1000);
    t1 = rdtsc();
    ierr = verify(buf, status0);
    printf("1 status = %d, expected %d, time = %ld, errors = %d, per item = %ld\n",
	   status,NREP*NW*1000,t1-t0,ierr,(t1-t0)/(NREP*NW*1000));
    printf("1 buf[0] = %d, buf[1] = %d, buf[last] = %d, buf[last+1] = %d\n",buf[0],buf[1],buf[status0-1],buf[status0]);

    ierr = MPI_Barrier(MPI_COMM_WORLD);
    status = FiolBufStatus(ptr) ;
    printf("1 shared buffer status is %d, expecting %d\n",status,BUFSZ-1);
    printf("1 first, in, out, limit = %d %d %d %d\n\n",ptr->first, ptr->in, ptr->out, ptr->limit);
    printf("1 stalled_insert = %ld, stalled_extract = %ld, block_insert = %ld, block_extract = %ld\n",
	   stalled_insert,stalled_extract,block_insert,block_extract);
  }else{  // all other ranks only participate in barriers
    ierr = MPI_Barrier(MPI_COMM_WORLD);
    ierr = MPI_Barrier(MPI_COMM_WORLD);
  }

  ierr = shmdt(ptr);
  ierr = MPI_Finalize();
}
#endif
