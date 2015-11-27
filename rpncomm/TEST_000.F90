subroutine rpn_comm_test_000
  implicit none
  integer ierr, myrank,totpes
  include 'mpif.h'
  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD,totpes,ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,myrank,ierr)
  print *,'I am PE ',myrank+1,' of ',totpes
  call test_null
  call mpi_finalize(ierr)
  stop
end
subroutine test_null
  use ISO_C_BINDING
  implicit none
  include 'RPN_COMM.inc'
  type(rpncomm_context) :: a
  type(rpncomm_ptr) :: b
  type(rpncomm_window) :: c
  type(rpncomm_file) :: d
  type(rpncomm_request) :: e
  type(rpncomm_datatype) :: f
  type(rpncomm_communicator) :: g
  type(rpncomm_group) :: h
  a%p = C_NULL_PTR
  b%p = C_NULL_PTR
  c%p = C_NULL_PTR
  d%p = C_NULL_PTR
  e%p = C_NULL_PTR
  f%p = C_NULL_PTR
  g%p = C_NULL_PTR
  h%p = C_NULL_PTR
  print *,RPN_COMM_is_null(a)
  print *,RPN_COMM_is_null(b)
  print *,RPN_COMM_is_null(c)
  print *,RPN_COMM_is_null(d)
  print *,RPN_COMM_is_null(e)
  print *,RPN_COMM_is_null(f)
  print *,RPN_COMM_is_null(g)
  print *,RPN_COMM_is_null(h)
return
end
	
