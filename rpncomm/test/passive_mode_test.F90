program passive_one_sided_test
! needs export  MPICH_RMA_OVER_DMAPP=1 at runtime on Cray XC systems
! MUST USE DMAPP on Cray XC systems
! ftn  -o passive_mode_test_cray passive_mode_test.F90 ifetch.o  -Wl,--whole-archive,-ldmapp,--no-whole-archive
  use ISO_C_BINDING
  implicit none
  include 'mpif.h'
  integer :: nprocess, myrank, ierr, window, i, temp, j
  type(C_PTR) :: base_addr
  integer, dimension(:), pointer :: winarray
  INTEGER(KIND=MPI_ADDRESS_KIND) TARGET_DISP, WIN_SIZE
  real*8 :: t1, t2, times(128)
  integer, external :: ifetch

  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD,nprocess,ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,myrank,ierr)

  WIN_SIZE = 1024
  call mpi_win_allocate(WIN_SIZE, 1, MPI_INFO_NULL, MPI_COMM_WORLD, base_addr, window, ierr)
  call c_f_pointer(base_addr,winarray,[256])

  call mpi_barrier(MPI_COMM_WORLD,ierr)
  print *,"I am process",myrank+1,' of',nprocess
  do j = 1,2
  winarray(:) = -1
  if(myrank == 0) then
    do i = 1, nprocess-1
      temp = i
      TARGET_DISP = 0
      t1 = mpi_wtime()
      call mpi_win_lock(MPI_LOCK_SHARED,i,0,window,ierr)
      call mpi_put(temp,1,MPI_INTEGER,i,TARGET_DISP,1,MPI_INTEGER,window,ierr)
      call mpi_win_unlock(i,window,ierr)
      t2 = mpi_wtime()
      times(i) = t2 - t1
    enddo
    print *," server pass",j
  else
    t1 = mpi_wtime()
    do while(ifetch(winarray(1)) <= 0)
    enddo
    t2 = mpi_wtime()
    times(j) = t2-t1
  endif
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  enddo

  if(myrank == 0) then
    print 101,(times(1:nprocess-1)*1000000)
  else
    print *,'in rank',myrank,' winarray(1)=',winarray(1),' t=',int(times(1:2)*1000000.)
  endif
101 format(20F6.0)

  call mpi_win_free(window,ierr)

  call mpi_finalize(ierr)
end program
