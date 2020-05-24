#include <mpi.h>

static int available = MPI_THREAD_SINGLE;  // assume worst case

int MPI_Init_thread(int *argc, char ***argv, int asked, int *provided){
  int status;
  int required = asked;

  if(required ==  MPI_THREAD_SINGLE || required == MPI_THREAD_FUNNELED) required = MPI_THREAD_SERIALIZED;
  status = PMPI_Init_thread(argc, argv, required, provided);
  available = *provided;
  return status;
}
int MPI_Init(int *argc, char ***argv){
  int status;

  status = PMPI_Init_thread(argc, argv, MPI_THREAD_SERIALIZED, &available);
  return status;
}
void rpn_comm_mpi_c_thread_provided(int *provided){
  *provided = available;
}
