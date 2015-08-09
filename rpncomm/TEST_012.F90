#if defined(TEST)
program my_test
call mpi_init(ierr)
call rpn_comm_test_012
call mpi_finalize(ierr)
stop
end
#endif
subroutine rpn_comm_test_012
  use ISO_C_BINDING
  implicit none
  include 'mpif.h'
  include 'RPN_COMM_types.inc'
  include 'RPN_COMM_ftoc.inc'

  type(rpncomm_datatype) :: my_type
  type(rpncomm_group) :: my_group
  type(rpncomm_window) :: my_win
  type(rpncomm_communicator) :: my_com
  type(rpncomm_file) :: my_file
  type(rpncomm_request) :: my_req
  type(rpncomm_info) :: my_info
  type(rpncomm_operator) :: my_op
  integer(C_INT) :: result

  my_type = NULL_rpncomm_datatype
  result = RPN_COMM_c2f(my_type)
  print *,'type result =',result,MPI_DATATYPE_NULL

  my_group = NULL_rpncomm_group
  result = RPN_COMM_c2f(my_group)
  print *,'group result =',result

  my_win = NULL_rpncomm_window
  result = RPN_COMM_c2f(my_win)
  print *,'window result =',result

  my_com = NULL_rpncomm_communicator
  result = RPN_COMM_c2f(my_com)
  print *,'com result =',result

  my_file = NULL_rpncomm_file
  result = RPN_COMM_c2f(my_file)
  print *,'file result =',result

  my_req = NULL_rpncomm_request
  result = RPN_COMM_c2f(my_req)
  print *,'request result =',result

  my_info = NULL_rpncomm_info
  result = RPN_COMM_c2f(my_info)
  print *,'info result =',result

  my_op = NULL_rpncomm_operator
  result = RPN_COMM_c2f(my_op)
  print *,'operator result =',result

  return
end subroutine rpn_comm_test_012