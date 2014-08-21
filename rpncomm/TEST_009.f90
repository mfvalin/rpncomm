subroutine rpn_comm_test_009
! call RPN_COMM_haloflip_test(halox, haloy, ni, nj, nk)
  call mpi_init(ierr)
  call RPN_COMM_haloflip_test(1,     2,      5,  6,  1)
  call RPN_COMM_haloflip_test(3,     3,     75, 27, 80)
  call mpi_finalize(ierr)
  stop
end
