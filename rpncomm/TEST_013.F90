!=======================================================================
module test_013
  use ISO_C_BINDING
#include <time_trace.hf>
  integer, save :: pe_mex, pe_mey
  type(time_context), save :: trace
end module test_013
subroutine rpn_comm_test_013
  use test_013
  implicit none
  include 'mpif.h'
  integer :: Pex, Pey, Pelocal, Petotal, halox, haloy, step, is, nbeads, nbent, nused, nprc
  external :: TestUserInit
  integer :: lni, lnj, gnk, gni, gnj, nk, lnk, lniymax, lniy
  character(len=128) :: RPN_COMM_TEST_CFG
  integer :: i, j, ier, errors, iter, tag, nsteps, iter2
  integer, external :: xch_halo_test, transpose_test
  integer(C_INT), dimension(:), pointer   :: local
  integer(C_INT), dimension(:,:), allocatable, target :: global
  integer(C_INT), dimension(:,:), allocatable :: blob
  type(C_PTR) :: local_c
  namelist /TEST_CFG/ lni, lnj, iter, nsteps, nk, iter2
  logical :: file_exists

  Pex = 0
  Pey = 0
  lni = 40     ! default configuration values
  lnj = 38
  iter = 5
  iter2 = 5
  nsteps = 1
  call RPN_COMM_init(TestUserInit,Pelocal,Petotal,Pex,Pey)
  pe_mex = mod(Pelocal,Pex)
  pe_mey = Pelocal / Pex
!   print *,'PE',Pelocal+1,' of',Petotal
!   nk = Pex * 2
  nk = max(80,Pex)     ! default configuration value
!   do while(nk < 60)
!     nk = nk + pex
!   enddo
  call get_environment_variable("RPN_COMM_TEST_CFG",RPN_COMM_TEST_CFG,i,ier)
  if(ier == 0) then  ! environment variable found, get dimensions
    read(RPN_COMM_TEST_CFG,*)lni, lnj, nk
  endif
  INQUIRE(FILE="RPN_COMM_TEST_CFG.txt", EXIST=file_exists)   ! check for explicit config file
  if(file_exists) then                                       ! file exists, read configuration overrides from namelist
    open(unit=55,file="RPN_COMM_TEST_CFG.txt",FORM="FORMATTED")
    read(unit=55,NML=TEST_CFG)
    if(Pelocal == 0) write(6,NML=TEST_CFG)
    close(unit=55)
  endif
  gni = lni * Pex
  gnj = lnj * pey
  gnk = max(nk,Pex)
  lnk = (nk + Pex - 1) / Pex        ! gnk distributed over Pex
  gnk = lnk * Pex                   ! make sure gnk is a multiple of Pex
  nk = gnk
  lniymax = (gni + Pey - 1) / Pey   ! gni distributed over Pey
  lniy = lniymax
  if(pe_mey == Pey - 1)  lniy = gni - lniymax*(Pey - 1)
  if(Pelocal == 0) then
    print 100,' Pelocal,Petotal,Pex,Pey =',Pelocal,Petotal,Pex,Pey
    print 100,' lni, lnj, lnk, lniy =', lni, lnj, lnk, lniy
    print 100,' gni, gnj, gnk =', gni, gnj, gnk
!     print 100,' mex, mey =',pe_mex, pe_mey
  endif
  if(gni - lniymax*(Pey - 1) < 0 ) then   ! OOPS cannot distribute properly over Pey
    if(Pelocal == 0) print 100,'ERROR: improper value for gni, lniymax, remainder =',lniymax,gni - lniymax*(Pex - 1)
    goto 777
  endif

  call time_trace_init(trace)

  do step = 1, nsteps
    call time_trace_step(trace, step)
    is = step
    errors = 0
    do halox = 1, 9, 2
      tag = 100 * halox
      haloy = halox
      if(iter > 0) then
        if(iter == 1) then  ! only one iteration
          errors = xch_halo_test(lni, lnj, nk, halox, haloy, gni, gnj,  0, tag) ! allocate, then deallocate
        else                ! more than one iteration
          errors = xch_halo_test(lni, lnj, nk, halox, haloy, gni, gnj,  2, tag) ! allocate, do not deallocate
        endif
      endif
      do i = 2, iter-1
        if(iter > 2) errors = errors + xch_halo_test(lni, lnj, nk, halox, haloy, gni, gnj, -2, tag)   ! neither allocate nor deallocate
      enddo
      if(iter > 1) then
        errors = errors + xch_halo_test(lni, lnj, nk, halox, haloy, gni, gnj, -1, tag) ! do not allocate, deallocate
      endif
    enddo
    tag = 1000
    if(iter2 > 0) then
      if(iter2 == 1) then  ! only one iteration
        errors = errors + transpose_test(lni, lnj, gnk, halox, haloy, gni, gnj, Pex, Pey, 0, tag) ! allocate, then deallocate
      else
        errors = errors + transpose_test(lni, lnj, gnk, halox, haloy, gni, gnj, Pex, Pey, 2, tag) ! allocate, do not deallocate
      endif
    endif
    do i = 2, iter2-1
      errors = errors + transpose_test(lni, lnj, gnk, halox, haloy, gni, gnj, Pex, Pey, -2, tag)   ! neither allocate nor deallocate
    enddo
    if(iter2 > 1) then
      errors = errors + transpose_test(lni, lnj, gnk, halox, haloy, gni, gnj, Pex, Pey, -1, tag) ! do not allocate, deallocate
    endif
    print 100,'step, iterations, errors =', step, iter,iter2, errors
!     if(Pelocal == 0) then
!     endif
  enddo

  call time_trace_dump_text(trace, 'time_list', Pelocal)

  local_c = time_trace_get_buffer_data(trace, nbeads, nbent, 1)
  call c_f_pointer(local_c, local, [nbent])  ! C pointer to Fortran pointer
  allocate(global(nbent,Petotal))
  call MPI_gather(local, nbent, MPI_INTEGER, global, nbent, MPI_INTEGER, 0, MPI_COMM_WORLD,ier)

  goto 777
  if(Pelocal == 0)then
    allocate(blob(2+2*Petotal,nbent/2))
    blob = 0
    nused = time_trace_expand(C_LOC(global(1,1)), nbent, blob, blob(3,1), 2+2*Petotal, nbent/2, 0)
    do i = 2, Petotal
      nused = time_trace_expand(C_LOC(global(1,i)), nbent, blob, blob(1+2*i,1), 2+2*Petotal, nbent/2, 1)
    enddo
    print *,'nused, size =',nused, nbent/2
    nprc = min(Petotal,5)
    do i = 1, abs(nused)
      if(blob(2,i) == -1) then   ! step flag, combine 2 unsigned 32 bit integers into a 64 bit integer
        print 102,blob(1,i), blob(2,i), &
                  ((i8_from_2_i4(blob(2*j+1,i), blob(2*j+2,i)), 0) , j=1,nprc)
      else
        print 101,blob(1:2+nprc*2,i)
      endif
    enddo
  endif

777 continue
  call RPN_COMM_finalize(ier)
100 format(A50,10I7)
101 format(12I10)
102 format(2I10,5(I18,I2))
end subroutine rpn_comm_test_013
!=======================================================================
subroutine TestUserInit(NX,NY) ! try to get NX,NY from file TEST.cfg if it exists
!=======================================================================
  implicit none
  integer, intent(OUT) :: nx, ny
  integer :: i, ier
  character(len=128) :: RPN_COMM_TEST_SHAPE
  call get_environment_variable("RPN_COMM_TEST_SHAPE",RPN_COMM_TEST_SHAPE,i,ier)
  if(ier == 0) then
    read(RPN_COMM_TEST_SHAPE,*)NX,NY
  else
    nx = 0
    ny = 0
    print *,'ERROR: environment variable RPN_COMM_TEST_SHAPE not found'
    print *,'       RPN_COMM_TEST_SHAPE="NPEX NPEY"'
  endif
  return
end subroutine TestUserInit
!=======================================================================
integer function xch_halo_test(lni, lnj, nk, halox, haloy, gni, gnj, mode, tag)
!=======================================================================
  use ISO_C_BINDING
  use test_013
  implicit none
  include 'mpif.h'
  external MPI_Barrier
  include 'RPN_COMM.inc'
  !
  integer, pointer, dimension(:,:,:), static :: localarray
  integer, intent(IN) :: lni, lnj, nk, halox, haloy, gni, gnj, mode
  integer, intent(INOUT) :: tag

  integer :: i, j, k, iter, ierr
  integer :: lminx, lmaxx, lminy, lmaxy
  integer :: minx1, maxx1, miny1, maxy1
  integer :: minx, maxx, miny, maxy
  integer :: npol_row, errors, expected, ii, jj
  logical :: periodx, periody
  !
  xch_halo_test=-1
  periodx = .false.
  periody = .false.
  npol_row = -1
  lminx = lni * pe_mex + 1
  lmaxx = lminx + lni -1
  lminy = lnj * pe_mey + 1
  lmaxy = lminy + lnj -1
  !
  minx = lminx-halox
  maxx = lmaxx+halox
  miny = lminy-haloy
  maxy = lmaxy+haloy

!   print 100,' tag =',tag
  tag = tag + 10
!   call mpi_barrier(MPI_COMM_WORLD,ierr)
!   call time_trace_barr(trace, tag, .true., MPI_COMM_WORLD, MPI_barrier)
  if(mode >= 0) then
    allocate(localarray(minx:maxx,miny:maxy,nk))
!     print 100,' allocated localarray',minx, maxx, miny, maxy, nk
  endif
  !
  localarray = 99999
  do k = 1,nk
  do j = lminy,lmaxy
  do i = lminx,lmaxx
    localarray(i,j,k) = mod(k,100) + 100*mod(j,1000) + 100000*mod(i,1000)
  enddo
  enddo
  enddo
!   do j = maxy, miny, -1
!     print 100,' avant ',localarray(minx:maxx,j,1)
!   enddo
!   print *,''
!   call mpi_barrier(MPI_COMM_WORLD,ierr)
  call time_trace_barr(trace, tag+2, .true., MPI_COMM_WORLD, MPI_barrier)
! SUBROUTINE RPN_COMM_xch_halo(g,minx,maxx,miny,maxy,ni,nj,nk,halox,haloy,periodx,periody,gni,npol_row)
  call RPN_COMM_xch_halo(localarray, &
                minx - lminx + 1, maxx - lminx + 1, &
                miny - lminy + 1, maxy - lminy + 1, &
                lni,lnj,nk,halox,haloy, &
                periodx,periody,gni,npol_row)

!   do j = maxy, miny, -1
!     print 100,' apres ',localarray(minx:maxx,j,1)
!   enddo
!   call mpi_barrier(MPI_COMM_WORLD,ierr)
  call time_trace_barr(trace, tag+4, .true., MPI_COMM_WORLD, MPI_barrier)
  errors=0
  miny1 = max(1,miny)
  maxy1 = min(gnj,maxy)
  minx1 = max(1,minx)
  maxx1 = min(gni,maxx)
!   print 100,' minx1, maxx1, miny1, maxy1 =',minx1, maxx1, miny1, maxy1
  do k = 1, nk
  do j = miny1, maxy1
  do i = minx1, maxx1
    expected = mod(k,100) + 100*mod(j,1000) + 100000*mod(i,1000)
    if(localarray(i,j,k) /= expected) errors=errors+1
  enddo
  enddo
  enddo

  if(abs(mode) <= 1) then
    deallocate(localarray)
!     print 100,' deallocated localarray'
  endif
!   call mpi_barrier(MPI_COMM_WORLD,ierr)
!   call time_trace_barr(trace, tag+6, .true., MPI_COMM_WORLD, MPI_barrier)

  xch_halo_test=errors
  return
100 format(A40,8I7)
end function xch_halo_test
!=======================================================================
!=======================================================================
integer function transpose_test(lni, lnj, gnk, halox, haloy, gni, gnj, Pex, Pey, mode, tag)
!=======================================================================
  use ISO_C_BINDING
  use test_013
  implicit none
  include 'mpif.h'
  external MPI_Barrier
  include 'RPN_COMM.inc'

  integer, intent(IN) :: lni, lnj, gnk, halox, haloy, gni, gnj, Pex, Pey, mode
  integer, intent(INOUT) :: tag
  integer :: errors, lnk, ierr, my_row, my_col, lniy
  integer, pointer, dimension(:,:,:),   static :: sx
  integer, pointer, dimension(:,:,:,:), static :: tx, sy, ty
  !
  errors = 0
  lnk = gnk / Pex
  lniy = (gni + Pey -1) / Pey
  if(mode >= 0) then
    allocate(sx(lni,lnj,gnk))
    allocate(tx(lni,lnj,lnk,Pex))
    allocate(sy(lniy,lnj,lnk,Pey))
    allocate(ty(lniy,lnj,lnk,Pey))
!     print *,'allocated with lni,lnj,gnk,Pex =',lni,lnj,gnk,Pex
  endif
  my_row = RPN_COMM_comm(RPN_COMM_EW)
  my_col = RPN_COMM_comm(RPN_COMM_NS)
  sx = 0
  tx = 0
  sy = 0
  ty = 0
  tag = tag + 10
  call time_trace_barr(trace, tag+2, .true., MPI_COMM_WORLD, MPI_barrier)
  call MPI_ALLTOALL(sx, lni*lnj*lnk, MPI_INTEGER, tx, lni*lnj*lnk, MPI_INTEGER, my_row, ierr)
  call time_trace_barr(trace, tag+3, .true., MPI_COMM_WORLD, MPI_barrier)
  call MPI_ALLTOALL(sy, lniy*lnj*lnk, MPI_INTEGER, ty, lniy*lnj*lnk, MPI_INTEGER, my_col, ierr)
  call time_trace_barr(trace, tag+4, .true., MPI_COMM_WORLD, MPI_barrier)
  call MPI_ALLTOALL(ty, lniy*lnj*lnk, MPI_INTEGER, sy, lniy*lnj*lnk, MPI_INTEGER, my_col, ierr)
  call time_trace_barr(trace, tag+5, .true., MPI_COMM_WORLD, MPI_barrier)
  call MPI_ALLTOALL(tx, lni*lnj*lnk, MPI_INTEGER, sx, lni*lnj*lnk, MPI_INTEGER, my_row, ierr)
  call time_trace_barr(trace, tag+6, .true., MPI_COMM_WORLD, MPI_barrier)
  if(abs(mode) <= 1) then
    deallocate(ty)
    deallocate(sy)
    deallocate(tx)
    deallocate(sx)
!     print *,'deallocated'
  endif
  transpose_test = errors
end function transpose_test
