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
  int32_t first;  // index in data of first point of circular buffer
  int32_t in;     // index in data of insertion point into circular buffer
  int32_t out;    // index in data of extraction point into circular buffer
  int32_t limit;  // index in data of point after last point of circular buffer
  int32_t data[1];
} fiol ;

static int64_t stalled_insert = 0;
static int64_t stalled_extract = 0;

// initialize fiol struct
// b        : pointer to control structure for this circular buffer
// datasize : number of 4 byte items in this circular buffer
void FiolBufInit(fiol *b, int32_t datasize){
  b->first = 0;
  b->in    = 0;
  b->out   = 0;
  b->limit = datasize;
}

// b        : pointer to control structure for this circular buffer
// userdata : data from user
// nw       : number of 4 byte items from user
// blocking : if 0, do not wait if buffer is full ; if non zero, wait
// return   : number of items succesfully inserted into buffer
int32_t FiolBufInsert(fiol volatile *b, void *userdata, int32_t nw, int32_t blocking){
  int status = 0;
  int nextin = b->in ;

  while(nw--){                                                    // not all data inserted
    nextin = ((nextin + 1) < b->limit) ? nextin + 1 : b->first;   // next insertion point
    while(nextin == b->out){                                      // circular buffer is full
      stalled_insert++;
      if(blocking == 0) return status;                            // non blocking return amount done
    }
    b->data[b->in] = ((int32_t *)userdata)[status]; // insert data
    status = status + 1;                                          // one more token done
    b->in = nextin;                                               // update in after successful insertion
  }
  return status;
}

// b        : pointer to control structure for this circular buffer
// userdata : data to user
// nw       : number of 4 byte items to user
// blocking : if 0, do not wait if buffer is empty ; if non zero, wait
// return   : number of items succesfully extracted from buffer
int32_t FiolBufExtract(fiol volatile *b, void *userdata, int32_t nw, int32_t blocking){
  int status = 0;

  while(nw--){                                                    // not all data inserted
    while(b->in == b->out){                                       // circular buffer is empty
      stalled_extract++;
      if(blocking == 0) return status;                            // non blocking return amount done
    }
    ((int32_t *)userdata)[status] = b->data[b->out]; // extract data
    status = status + 1;                                          // one more token done
    b->out = (b->out+1 < b->limit) ? b->out+1 : b->first;         // next extraction point
  }
  return status;
}

// b        : pointer to control structure for this circular buffer
// if buffer is empty, return +(number of available slots) (buffer capacity)
// if buffer is full, return 0
// otherwise return -(number of free slots)
int32_t FiolBufStatus(fiol *b){
  int status ;
//   printf("first, in, out, limit = %d %d %d %d\n",b->first, b->in, b->out, b->limit);
  if(b->in == b->out) {
    status = b->limit - b->first ;                 //  buffer is empty, return buffer capacity
  }else{
    if(b->in > b->out) {                           // in - out tokens in buffer
      status = -(b->limit - b->first - (b->in - b->out)); // number of free slots
    }else{                                         // out - in - 1 free slots
      status = -(b->out - b->in ) ;                // number of free slots
    }
  }
  return status;
}

#if defined(SELF_TEST)
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <sys/stat.h>

static uint64_t rdtsc(void) {   // version rapide "out of order"
  uint32_t lo, hi;
  __asm__ volatile ("rdtsc"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
  return (uint64_t)lo | (((uint64_t)hi) << 32);
}

#define BUFSZ 128
#define NW 15
main(int argc, char **argv){
  int32_t buf[NW*1024*1024];
  int32_t status;
  int ierr, localrank, id;
  fiol *ptr = NULL;
  struct shmid_ds shm_buf;
  size_t size;
  int i, nrep;
  uint64_t t0 , t1;

  ierr = MPI_Init( &argc, &argv );
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&localrank);

  if(localrank == 0){
    size = 123456 ;
    id = shmget(IPC_PRIVATE,size,IPC_CREAT|S_IRUSR|S_IWUSR);  /* rank 0 allocates shared memory segment */
    ptr = shmat(id,NULL,0);
    shmctl(id,IPC_RMID,&shm_buf);      /* mark segment for deletion preventively (only works on linux) */
    for(i=0 ; i<BUFSZ ; i++) ptr->data[i] = 0;  // force pages to be properly instanciated in memory
    FiolBufInit(ptr, BUFSZ);

    for(i=0 ; i<sizeof(buf)/sizeof(int32_t) ; i++) buf[i] = i;
    t0 = rdtsc();
    status = FiolBufInsert(ptr, buf, NW, 1);
    t1 = rdtsc();
    printf("status = %d, expected %d, time = %ld\n",status,NW,t1-t0);
    status = FiolBufStatus(ptr) ;
    printf("buffer status is %d, expecting %d\n",status,-(BUFSZ-NW));
    printf("first, in, out, limit = %d %d %d %d\n\n",ptr->first, ptr->in, ptr->out, ptr->limit);

    ierr = MPI_Bcast(&id,1,MPI_INTEGER,0,MPI_COMM_WORLD);
    printf("0 shared memory segment ID = %d, ptr = %p\n\n",id,ptr);
    status = FiolBufStatus(ptr) ;
    printf("shared buffer status is %d, expecting %d\n",status,-(BUFSZ-NW));
    printf("first, in, out, limit = %d %d %d %d\n\n",ptr->first, ptr->in, ptr->out, ptr->limit);

    ierr = MPI_Barrier(MPI_COMM_WORLD);

    t0 = rdtsc();
    status = FiolBufInsert(ptr, buf, NW*4, 1);
    t1 = rdtsc();
    printf("status = %d, expected %d, time = %ld\n",status,NW*4,t1-t0);
    t0 = rdtsc();
    status = FiolBufInsert(ptr, buf, NW*1000, 1);
    t1 = rdtsc();
    printf("status = %d, expected %d, time = %ld\n",status,NW*1000,t1-t0);

    ierr = MPI_Barrier(MPI_COMM_WORLD);
    status = FiolBufStatus(ptr) ;
    printf("shared buffer status is %d, expecting %d\n",status,BUFSZ);
    printf("first, in, out, limit = %d %d %d %d\n\n",ptr->first, ptr->in, ptr->out, ptr->limit);
    printf("buf[0] = %d, buf[last] = %d, buf[last+1] = %d\n",buf[0],buf[NW*1000-1],buf[NW*1000]);
    
  }else if(localrank == 1){
    ierr = MPI_Bcast(&id,1,MPI_INTEGER,0,MPI_COMM_WORLD);
    ptr = shmat(id,NULL,0);
    printf("%d shared memory segment ID = %d, ptr = %p\n\n",localrank,id,ptr);
    status = FiolBufStatus(ptr) ;
    printf("shared buffer status is %d, expecting %d\n",status,-(BUFSZ-NW));
    printf("first, in, out, limit = %d %d %d %d\n\n",ptr->first, ptr->in, ptr->out, ptr->limit);

    for(i=0 ; i<sizeof(buf)/sizeof(int32_t) ; i++) buf[i] = 999999;
    ierr = MPI_Barrier(MPI_COMM_WORLD);

    t0 = rdtsc();
    status = FiolBufExtract(ptr, buf, NW, 1);
    t1 = rdtsc();
    printf("status = %d, expected %d, time = %ld\n",status,NW,t1-t0);
    printf("buf[0] = %d, buf[1] = %d, buf[last] = %d, buf[last+1] = %d\n",buf[0],buf[1],buf[NW-1],buf[NW]);

    t0 = rdtsc();
    status = FiolBufExtract(ptr, buf, NW*4, 1);
    t1 = rdtsc();
    printf("status = %d, expected %d, time = %ld\n",status,NW*4,t1-t0);
    printf("buf[0] = %d, buf[1] = %d, buf[last] = %d, buf[last+1] = %d\n",buf[0],buf[1],buf[NW*4-1],buf[NW*4]);

    t0 = rdtsc();
    status = FiolBufExtract(ptr, buf, NW*1000, 1);
    t1 = rdtsc();
    printf("status = %d, expected %d, time = %ld\n",status,NW*1000,t1-t0);
    printf("buf[0] = %d, buf[1] = %d, buf[last] = %d, buf[last+1] = %d\n",buf[0],buf[1],buf[NW*1000-1],buf[NW*1000]);

    ierr = MPI_Barrier(MPI_COMM_WORLD);
    status = FiolBufStatus(ptr) ;
    printf("shared buffer status is %d, expecting %d\n",status,BUFSZ);
    printf("first, in, out, limit = %d %d %d %d\n\n",ptr->first, ptr->in, ptr->out, ptr->limit);
  }else{  // all other ranks only participate in barriers
    ierr = MPI_Barrier(MPI_COMM_WORLD);
    ierr = MPI_Barrier(MPI_COMM_WORLD);
  }
  printf("stalled_insert = %ld, stalled_extract=%ld\n",stalled_insert,stalled_extract);

  ierr = shmdt(ptr);
  ierr = MPI_Finalize();
}
#endif
