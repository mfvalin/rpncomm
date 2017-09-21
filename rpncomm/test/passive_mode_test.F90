program passive_one_sided_test
! needs export  MPICH_RMA_OVER_DMAPP=1  on Cray XC systems
! MUST USE DMAPP on Cray XC systems
! ftn  -o passive_mode_test_cray passive_mode_test.F90 ifetch.o  -Wl,--whole-archive,-ldmapp,--no-whole-archive
  use ISO_C_BINDING
  implicit none
  include 'mpif.h'
  integer :: nprocess, myrank, ierr, window, i, temp
  type(C_PTR) :: base_addr
  integer, dimension(:), pointer :: winarray
  INTEGER(KIND=MPI_ADDRESS_KIND) TARGET_DISP
  real*8 :: t1, t2
  integer, external :: ifetch

  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD,nprocess,ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,myrank,ierr)
  print *,"I am process",myrank+1,' of',nprocess

  call mpi_win_allocate(1024, 1, MPI_INFO_NULL, MPI_COMM_WORLD, base_addr, window, ierr)
  call c_f_pointer(base_addr,winarray,[256])
  winarray(:) = -1

  call mpi_barrier(MPI_COMM_WORLD,ierr)
  if(myrank == 0) then
    do i = 1, nprocess-1
      temp = i
      call mpi_win_lock(MPI_LOCK_SHARED,temp,0,window,ierr)
      TARGET_DISP = 0
      call mpi_put(temp,1,MPI_INTEGER,temp,TARGET_DISP,1,MPI_INTEGER,window,ierr)
      call mpi_win_unlock(i,window,ierr)
    enddo
  else
    t1 = mpi_wtime()
!     print *,'in rank',myrank,'entering while'
    do while(ifetch(winarray(1)) <= 0)
    enddo
    t2 = mpi_wtime()
    t2 = t2-t1
    print *,'in rank',myrank,' winarray(1)=',winarray(1),' t=',t2*1000000.
  endif
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  call mpi_win_free(window,ierr)

  call mpi_finalize(ierr)
end program
