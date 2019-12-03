!=======================================================================
module test_013
  integer, save :: pe_mex, pe_mey
end module test_013
subroutine rpn_comm_test_013
  use test_013
  implicit none
  integer :: Pex, Pey, Pelocal, Petotal, halox, haloy
  external :: TestUserInit
  integer :: lni, lnj, lnk, gni, gnj, nk
  character(len=128) :: RPN_COMM_TEST_CFG
  integer :: i, ier, errors, iter, tag
  integer, external :: xch_halo_test

  Pex = 0
  Pey = 0
  lni = 40
  lnj = 38
  halox = 3
  haloy = 3
  iter = 8
  call RPN_COMM_init(TestUserInit,Pelocal,Petotal,Pex,Pey)
  pe_mex = mod(Pelocal,Pex)
  pe_mey = Pelocal / Pex
!   print *,'PE',Pelocal+1,' of',Petotal
  nk = Pex * 2
!   do while(nk < 60)
!     nk = nk + pex
!   enddo
  call get_environment_variable("RPN_COMM_TEST_CFG",RPN_COMM_TEST_CFG,i,ier)
  if(ier == 0) then
    read(RPN_COMM_TEST_CFG,*)lni, lnj, nk, halox, haloy
  endif
  gni = lni * Pex
  gnj = lnj * pey
  lnk = nk
  if(Pelocal == 0) then
    print 100,' Pelocal,Petotal,Pex,Pey =',Pelocal,Petotal,Pex,Pey
    print 100,' lni, lnj, lnk =', lni, lnj, lnk
    print 100,' gni, gnj, nk =', gni, gnj, nk
!     print 100,' mex, mey =',pe_mex, pe_mey
  endif  

  tag = 100
  do i = 1, iter
    if(i == 1) errors = xch_halo_test(lni, lnj, nk, halox, haloy, gni, gnj,  2, tag) ! allocate, do not deallocate
    errors = errors + xch_halo_test(lni, lnj, nk, halox, haloy, gni, gnj, -2, tag)   ! neither
    if(i == iter) errors = errors + xch_halo_test(lni, lnj, nk, halox, haloy, gni, gnj, -1, tag) ! do not allocate, deallocate
  enddo
  print 100,' iterations, halox, haloy, pe, errors =',  iter+2, halox, haloy, Pelocal, errors

  halox = 9
  haloy = 9

  do i = 1, iter
    if(i == 1) errors = xch_halo_test(lni, lnj, nk, halox, haloy, gni, gnj,  2, tag) ! allocate, do not deallocate
    errors = errors + xch_halo_test(lni, lnj, nk, halox, haloy, gni, gnj, -2, tag)   ! neither
    if(i == iter) errors = errors + xch_halo_test(lni, lnj, nk, halox, haloy, gni, gnj, -1, tag) ! do not allocate, deallocate
  enddo
  print 100,' iterations, halox, haloy, pe, errors =',  iter+2, halox, haloy, Pelocal, errors

  call RPN_COMM_finalize(ier)
100 format(A40,8I7)
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

  print 100,' tag =',tag
  tag = tag + 10
  call mpi_barrier(MPI_COMM_WORLD,ierr)
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
  call mpi_barrier(MPI_COMM_WORLD,ierr)
! SUBROUTINE RPN_COMM_xch_halo(g,minx,maxx,miny,maxy,ni,nj,nk,halox,haloy,periodx,periody,gni,npol_row)
  call RPN_COMM_xch_halo(localarray, &
                minx - lminx + 1, maxx - lminx + 1, &
                miny - lminy + 1, maxy - lminy + 1, &
                lni,lnj,nk,halox,haloy, &
                periodx,periody,gni,npol_row)

!   do j = maxy, miny, -1
!     print 100,' apres ',localarray(minx:maxx,j,1)
!   enddo
  call mpi_barrier(MPI_COMM_WORLD,ierr)
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
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  xch_halo_test=errors
  return
100 format(A40,8I7)
end function xch_halo_test
!=======================================================================
