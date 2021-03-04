#if defined(TEST)
program my_test
call mpi_init(ierr)
call rpn_comm_test_014
call mpi_finalize(ierr)
stop
end
#endif
subroutine rpn_comm_test_014
  use ISO_C_BINDING
  implicit none
  include 'mpif.h'
  include 'RPN_COMM.inc'
  integer :: ierr, narg, arglen, rank, pop
  external :: Userinit
  integer :: Pelocal,Petotal,Pex,Pey,MultiGrids,Grids,Io
  character(len=128) :: AppID, string

  call MPI_Init(ierr)
  call MPI_Comm_Rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_Size(MPI_COMM_WORLD, pop, ierr)
  do narg = 1, 6
    if(narg == 1) then
      call get_command_argument(1 , string, arglen, ierr)
      if(ierr .ne. 0) goto 999
      read(string,*,err=999) Pex
!       print *,'pex =',pex
    endif
    if(narg == 2) then
      call get_command_argument(2 , string, arglen, ierr)
      if(ierr .ne. 0) goto 999
      read(string,*,err=999) Pey
!       print *,'pey =',pey
    endif
    if(narg == 3) then
      call get_command_argument(3 , string, arglen, ierr)
      if(ierr .ne. 0) goto 999
      read(string,*,err=999) Multigrids
!       print *,'Multigrids =',Multigrids
    endif
    if(narg == 4) then
      call get_command_argument(4 , string, arglen, ierr)
      if(ierr .ne. 0) goto 999
      read(string,*,err=999) Grids
!       print *,'Grids =',Grids
    endif
    if(narg == 5) then
      call get_command_argument(5 , AppID, arglen, ierr)
      if(ierr .ne. 0) goto 999
!       print *,'AppID = "',trim(AppID)//'"'
    endif
    if(narg == 6) then
      call get_command_argument(6 , string, arglen, ierr)
      if(ierr .ne. 0) goto 999
      read(string,*,err=999) Io
!       print *,'Io =',Io
    endif
  enddo

  if(rank >= pop/2) AppID = "<N02>"
  ierr = RPN_COMM_init_all_levels(Userinit,Pelocal,Petotal,Pex,Pey,Multigrids,Grids,AppID,Io)

777 continue
  call RPN_COMM_finalize(ierr)
  stop
999 continue
  print *,"ERROR in arguments",ierr
  go to 777
  return
end subroutine rpn_comm_test_014

subroutine Userinit
  return
end subroutine Userinit\
