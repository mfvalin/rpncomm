#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>

#define N 128

// needs export  MPICH_RMA_OVER_DMAPP=1  on Cray XC systems
//  MUST USE DMAPP on Cray XC systems
//  cc -o one_sided_passive one_sided_passive.c  -Wl,--whole-archive,-ldmapp,--no-whole-archive

int main( int argc, char **argv )
{
    int buf[N];
    int rank, i ;
    volatile int *winbuf;
    MPI_Win window;
    MPI_Aint winsize;

    MPI_Init( &argc, &argv );
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    winsize = 4096;
    MPI_Win_allocate(winsize,4,MPI_INFO_NULL,MPI_COMM_WORLD,&winbuf,&window);
    for ( i=0; i< N ; i++ ) {
      buf[i]    =  i *1000;
      winbuf[i] = -i;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if(rank == 0) {
      MPI_Win_lock(MPI_LOCK_SHARED,1-rank,0,window);
      printf("writer locked rank %d\n",1-rank);; fflush(stdout);
      MPI_Put(buf,N,MPI_INT,1-rank,0,N,MPI_INT,window);
      printf("writer put to rank %d\n",1-rank);; fflush(stdout);
      MPI_Win_unlock(1-rank,window);
      printf("writer unlocked rank %d\n",1-rank);; fflush(stdout);
    }else{
      while(winbuf[1] < 0){
        printf("."); fflush(stdout);
        usleep(1000);
      }
    printf("\n");
    printf("reader locked rank %d\n",rank);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank ==1)for ( i=0; i< N ; i++ ) { printf("  %4d",winbuf[i]); } ; printf("\n");

    MPI_Win_free(&window);
    MPI_Finalize();
    return 0;
}

