!===================================================================
module ftn_mpi_init
  implicit none
  include 'mpif.h'
  integer, save :: provided = MPI_THREAD_SINGLE  ! save provided threaded version, assume worst case
end module ftn_mpi_init
!===================================================================
subroutine MPI_Init(ierr)  ! override library version
  use ftn_mpi_init         ! to force usage of threaded version
  implicit none
  integer, intent(OUT) :: ierr
  call PMPI_Init_thread(MPI_THREAD_FUNNELED, provided, ierr)
  return
end subroutine MPI_Init
!===================================================================
subroutine MPI_Init_thread(required, available, ierr)
  use ftn_mpi_init
  implicit none
  integer, intent(OUT) :: ierr, available
  integer, intent(IN) :: required
  integer :: asked
  asked = required
  if(asked == MPI_THREAD_SINGLE .or. asked == MPI_THREAD_FUNNELED) asked = MPI_THREAD_SERIALIZED
  call PMPI_Init_thread(asked, provided, ierr)
  available = provided
  return
end subroutine MPI_Init_thread
!===================================================================
subroutine rpn_comm_mpi_thread_provided(what)  ! make provided version available
  use ftn_mpi_init                             ! information available
  implicit none
  integer, intent(OUT) :: what
  what = provided
end subroutine rpn_comm_mpi_thread_provided
!===================================================================
