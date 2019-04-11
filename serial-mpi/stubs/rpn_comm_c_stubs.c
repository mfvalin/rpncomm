#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
int MPI_Comm_rank(MPI_Comm comm, int *rank){
   *rank = 0;
   return(MPI_SUCCESS);
}
int MPI_Comm_size(MPI_Comm comm, int *size){
   *size = 1;
   return(MPI_SUCCESS);
}
int MPI_Allgather(const void *outx, int nout, MPI_Datatype outtype,const void *inx, int nin, MPI_Datatype intype, MPI_Comm comm){
   int *out = (int *) outx;
   int *in = (int *) inx;
/*   if( nin != 1 || nout != 1 || outtype != MPI_INTEGER || intype != MPI_INTEGER ) {  */
   if( nin != 1 || nout != 1 ) {
/*      fprintf(stderr,"ERROR: MPI_Allgather, number of elements is not one or type is not MPI_INTEGER\n");  */
      fprintf(stderr,"ERROR: MPI_Allgather, number of elements is not one \n");
      exit(1);
   *out = *in;
   }
   return(MPI_SUCCESS);
}
int MPI_Barrier(MPI_Comm comm){
  return(MPI_SUCCESS);
}
int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm){
  return(MPI_SUCCESS);
}
int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm){
  fprintf(stderr,"ERROR: MPI_Allreduce not supported in serial mode\n");
  exit(1);
}

int MPI_Comm_c2f(int i){
 return i;
}

int MPI_Comm_f2c(int i){
 return i;
}

int MPI_File_c2f(int i){
 return i;
}

int MPI_File_f2c(int i){
 return i;
}

int MPI_Group_c2f(int i){
 return i;
}

int MPI_Group_f2c(int i){
 return i;
}

int MPI_Info_c2f(int i){
 return i;
}

int MPI_Info_f2c(int i){
 return i;
}

int MPI_Op_c2f(int i){
 return i;
}

int MPI_Op_f2c(int i){
 return i;
}

int MPI_Request_c2f(int i){
 return i;
}

int MPI_Request_f2c(int i){
 return i;
}

int MPI_Type_c2f(int i){
 return i;
}

int MPI_Type_f2c(int i){
 return i;
}

int MPI_Win_c2f(int i){
 return i;
}

int MPI_Win_f2c(int i){
 return i;
}

