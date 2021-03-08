#include <stdio.h>
#include <stdlib.h>
#include "mpi_stub.h"
int MPI_Comm_rank(int comm, int *rank){
   *rank = 0;
   return(MPI_SUCCESS);
}
int MPI_Comm_size(int comm, int *size){
   *size = 1;
   return(MPI_SUCCESS);
}
int MPI_Allgather(void *outx, int nout, int outtype, void *inx, int nin, int intype, int comm){
   int *out = outx;
   int *in = inx;
/*   if( nin != 1 || nout != 1 || outtype != MPI_INTEGER || intype != MPI_INTEGER ) {  */
   if( nin != 1 || nout != 1 ) {
/*      fprintf(stderr,"ERROR: MPI_Allgather, number of elements is not one or type is not MPI_INTEGER\n");  */
      fprintf(stderr,"ERROR: MPI_Allgather, number of elements is not one \n");
      exit(1);
   *out = *in;
   }
   return(MPI_SUCCESS);
}

