!/! RPN_COMM - Library of useful routines for C and FORTRAN programming
! ! Copyright (C) 1975-2017  Division de Recherche en Prevision Numerique
! !                          Environnement Canada
! !
! ! This library is free software; you can redistribute it and/or
! ! modify it under the terms of the GNU Lesser General Public
! ! License as published by the Free Software Foundation,
! ! version 2.1 of the License.
! !
! ! This library is distributed in the hope that it will be useful,
! ! but WITHOUT ANY WARRANTY; without even the implied warranty of
! ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! ! Lesser General Public License for more details.
! !
! ! You should have received a copy of the GNU Lesser General Public
! ! License along with this library; if not, write to the
! ! Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! ! Boston, MA 02111-1307, USA.
! !/
!
! map of regions in array f 
! A, B, C, D external corner regions to get from neighbors
! 1, 2, 3, 4, 5, 6, 7, 8 regions from which neighbors get their external corners
! array f has a useful size lni x lnj, with haloes along x (hx) and y (hy)
!
!           +---+---+                      +---+---+----(lnj+hy)
!           | A | 1 |                      | 2 | B |
!           +---+---+----------------------+---+---+----(lnj)
!           | 8 |                             /| 3 |
!           +---+                            / +---+----(lnj-hy+1)
!               |                    (lni.lnj) |
!               |                              |
!               |                              |
!               |           array f            |
!               |                              |
!               |                              |
!           +---+ (1,1)                        +---+----(hy)
!           | 7 |/                             | 4 |
!           +---+---+----------------------+---+---+----(1)
!           | D | 6 |                      | 5 | C |
!           +---+---+----------------------+---+---+----(1-hy)
!           |   |   |                      |   |   |
!           |   |  (hx)             (lni-hx+1) |   |
!           |  (1)                           (lni) |
!         (1-hx)                                (lni+hx)
!

#define Z1 f(1       :hx      ,lnj+1   :lnj+hy,:)
#define Z2 f(lni-hx+1:lni     ,lnj+1   :lnj+hy,:)
#define Z3 f(lni+1   :lni+hx  ,lnj-hy+1:lnj   ,:)
#define Z4 f(lni+1   :lni+hx  ,1       :hy    ,:)
#define Z5 f(lni-hx+1:lni     ,1-hy    :0     ,:)
#define Z6 f(1       :hx      ,1-hy    :0     ,:)
#define Z7 f(1-hx    :0       ,1       :hy    ,:)
#define Z8 f(1-hx    :0       ,lnj-hy+1:lnj   ,:)

#define A f(1-hx    :0       ,lnj+1   :lnj+hy,:)
#define B f(lni+1   :lni+hx  ,lnj+1   :lnj+hy,:)
#define C f(lni+1   :lni+hx  ,1-hy    :0     ,:)
#define D f(1-hx    :0       ,1-hy    :0     ,:)

#define TAG 0

#if defined(SELF_TEST)

#define rpn_comm_test rpn_comm

module rpn_comm_test
#include <mpif.h>
  integer, dimension(:,:), pointer :: pe_id
  integer :: pe_mex, pe_mey, pe_grid, pe_nx, pe_ny
  logical :: bnd_north, bnd_south, bnd_east, bnd_west
end module rpn_comm
program test
  use rpn_comm
  implicit none
  integer :: NX=25
  integer :: NY=15
  integer :: NK=80
  integer :: HX=2
  integer :: HY=2
  integer :: ierror, pe_me, siz, i, j, k, errors, lni, lnj
  integer, dimension(:,:,:), pointer :: f
  character(len=128) :: argument
  real*8 :: t1, t2, t3
!   logical interior
  
  call MPI_init(ierror)
  argument = "0"
  CALL GET_COMMAND_ARGUMENT(1 , argument)
  read(argument,*) pe_nx

  pe_grid = MPI_COMM_WORLD
  call MPI_comm_size(pe_grid,siz,ierror)
  call MPI_comm_rank(pe_grid,pe_me,ierror)
!   print *,'rank',pe_me+1,' of',siz
  call MPI_barrier(MPI_COMM_WORLD,ierror)

  allocate(f(1-HX:NX+HX,1-HY:NY+HY,NK))
  pe_ny = siz / pe_nx
  if(pe_nx * pe_ny .ne. siz) then
    if(pe_me == 0) print *,'pe_nx * pe_ny .ne. siz',pe_nx,pe_ny,siz
    goto 1
  endif
  if(pe_me == 0) print *,'pe_nx , pe_ny , siz',pe_nx,pe_ny,siz
  pe_mex = mod(pe_me,pe_nx)
  pe_mey = pe_me / pe_nx

  bnd_north = (pe_mey == (pe_ny - 1))
  bnd_south = (pe_mey == 0)
  bnd_east  = (pe_mex == (pe_nx - 1))
  bnd_west  = (pe_mex == 0)
!   interior = .not. (bnd_west .or. bnd_east .or. bnd_south .or. bnd_north)
  allocate(pe_id(-1:pe_nx,-1:pe_ny))
  pe_id = -1
  lni = nx
  lnj = ny

  k = 0
  do j = 0,pe_ny-1
  do i = 0,pe_nx-1
    pe_id(i,j) = k
    k = k + 1
  enddo
  enddo

  if(pe_me == 0 .and. pe_nx <= 10 .and. pe_ny <= 10) then
    do j = pe_ny,-1, -1
      print 101,pe_id(:,j)
    enddo
  endif
  call MPI_barrier(MPI_COMM_WORLD,ierror)

  do k = 1,nk
    do j = 1-hy , ny+hy
    do i = 1-hx , nx+hx
      f(i,j,k) = (i + nx * pe_mex + 10) * 1000 + (j + ny * pe_mey + 10)
    enddo
    enddo
  enddo
  if(bnd_north .and. .not.  bnd_west) A = 999999
  if(bnd_west  .and. .not. bnd_north) A = 999999
  if(bnd_north .and. .not.  bnd_east) B = 999999
  if(bnd_east  .and. .not. bnd_north) B = 999999
  if(bnd_south .and. .not.  bnd_east) C = 999999
  if(bnd_east  .and. .not. bnd_south) C = 999999
  if(bnd_west  .and. .not. bnd_south) D = 999999
  if(bnd_south .and. .not.  bnd_west) D = 999999

  if(nx <= 5 .and. ny <= 3 .and. pe_nx*pe_ny <= 6) then
    do j = ny+hy, 1-hy, -1
      print 101,f(:,j,1)
    enddo
  endif
  call MPI_barrier(MPI_COMM_WORLD,ierror)

  call RPN_COMM_propagate_boundary_circular(f,1-hx,nx+hx,1-hy,ny+hy,nx,ny,nk,hx,hy)  ! to prime the MPI pump
  t1 = MPI_wtime()
  call RPN_COMM_propagate_boundary_circular(f,1-hx,nx+hx,1-hy,ny+hy,nx,ny,nk,hx,hy)
!   call RPN_COMM_propagate_boundary(f,1-hx,nx+hx,1-hy,ny+hy,nx,ny,nk,hx,hy)
  t2 = MPI_wtime()
  call MPI_barrier(MPI_COMM_WORLD,ierror)
  t3 = MPI_wtime()
!   call RPN_COMM_propagate_boundary_circular(f,1-hx,nx+hx,1-hy,ny+hy,nx,ny,nk,hx,hy)
  errors = 0
  do k = 1 , nk
  do j = 1-hy , ny+hy
  do i = 1-hx , nx+hx
    if(f(i,j,k) .ne. (i + nx * pe_mex + 10) * 1000 + (j + ny * pe_mey + 10)) errors = errors + 1
  enddo
  enddo
  enddo
  print *,'pe =',pe_me,'time =',nint((t3-t1)*1000000),nint((t2-t1)*1000000),' errors =',errors
  if(nx <= 5 .and. ny <= 3 .and. pe_nx*pe_ny <= 6) then
    print *," "
    do j = ny+hy, 1-hy, -1
      print 101,f(:,j,1)
    enddo
  endif
1 continue
  call MPI_finalize(ierror)
101 format(15I8.6)
end program
#endif

subroutine RPN_COMM_propagate_boundary(f,minx,maxx,miny,maxy,lni,lnj,nk,hx,hy)
  use rpn_comm
  implicit none
  integer, intent(IN) :: minx,maxx,miny,maxy,lni,lnj,nk,hx,hy
  integer, dimension(minx:maxx,miny:maxy,nk), intent(INOUT) :: f

  logical :: north, south, east, west
  integer :: northpe, southpe, eastpe, westpe
  integer :: ierror
  integer, dimension(MPI_STATUS_SIZE) :: statn, stats, state, statw ! status for wait
  integer :: reqn1, reqs1, reqe1, reqw1                             ! irecv requests
  integer :: reqn2, reqs2, reqe2, reqw2                             ! irecv requests
  integer, dimension(hx,hy,nk) :: ta, tb, tc, td                    ! recv buffers
  integer, dimension(hx,hy,nk) :: temp                              ! send buffer
  integer :: npts

  north=(bnd_north)
  northpe=pe_id(pe_mex,pe_mey+1)
  south=(bnd_south)
  southpe=pe_id(pe_mex,pe_mey-1)
  east=(bnd_east)
  eastpe=pe_id(pe_mex+1,pe_mey)
  west=(bnd_west)
  westpe=pe_id(pe_mex-1,pe_mey)
  npts = hx * hy * nk
! print *,pe_mex,pe_mey,north,south,east,west,npts
! print *,northpe,southpe,eastpe,westpe
  if( .not. (north .or. south .or. east .or. west) ) return ! not on boundary, nothing to do
  if(pe_nx == 1 .and. pe_ny ==1) return ! single tile, nothing to do
! post appropriate nonblocking receives
  if(north) then
    if(.not. east) then  ! get B from east
      call  MPI_irecv(tb, npts, MPI_INTEGER, eastpe, TAG, pe_grid, reqe1, ierror)
    endif
    if(.not. west) then  ! get A from west
      call  MPI_irecv(ta, npts, MPI_INTEGER, westpe, TAG, pe_grid, reqw1, ierror)
    endif
  endif
  if(south) then
    if(.not. east) then  ! get C from east
      call  MPI_irecv(tc, npts, MPI_INTEGER, eastpe, TAG, pe_grid, reqe2, ierror)
    endif
    if(.not. west) then  ! get D from west
      call  MPI_irecv(td, npts, MPI_INTEGER, westpe, TAG, pe_grid, reqw2, ierror)
    endif
  endif
  if(east)  then
    if(.not. north) then  ! get B from north
      call  MPI_irecv(tb, npts, MPI_INTEGER, northpe, TAG, pe_grid, reqn1, ierror)
    endif
    if(.not. south) then  ! get C from south
      call  MPI_irecv(tc, npts, MPI_INTEGER, southpe, TAG, pe_grid, reqs1, ierror)
    endif
  endif
  if(west)  then
    if(.not. north) then  ! get A from north
      call  MPI_irecv(ta, npts, MPI_INTEGER, northpe, TAG, pe_grid, reqn2, ierror)
    endif
    if(.not. south) then  ! get D from south
      call  MPI_irecv(td, npts, MPI_INTEGER, southpe, TAG, pe_grid, reqs2, ierror)
    endif
  endif
! extract and blocking send to appropriate destination
temp = 0
  if(north) then
    if(.not. east) then  ! send 2 to east
      temp = Z2
      call mpi_send(temp, npts, MPI_INTEGER, eastpe, TAG, pe_grid, ierror)
    endif
    if(.not. west) then  ! send 1 to west
      temp = Z1
      call mpi_send(temp, npts, MPI_INTEGER, westpe, TAG, pe_grid, ierror)
    endif
  endif
  if(south) then
    if(.not. east) then  ! send 5 to east
      temp = Z5
      call mpi_send(temp, npts, MPI_INTEGER, eastpe, TAG, pe_grid, ierror)
    endif
    if(.not. west) then  ! send 6 to west
      temp = Z6
      call mpi_send(temp, npts, MPI_INTEGER, westpe, TAG, pe_grid, ierror)
    endif
  endif
  if(east)  then
    if(.not. north) then  ! send 3 to north
      temp = Z3
      call mpi_send(temp, npts, MPI_INTEGER, northpe, TAG, pe_grid, ierror)
    endif
    if(.not. south) then  ! send 4 to south
      temp = Z4
      call mpi_send(temp, npts, MPI_INTEGER, southpe, TAG, pe_grid, ierror)
    endif
  endif
  if(west)  then
    if(.not. north) then  ! send 8 to north
      temp = Z8
      call mpi_send(temp, npts, MPI_INTEGER, northpe, TAG, pe_grid, ierror)
    endif
    if(.not. south) then  ! send 7 to south
      temp = Z7
      call mpi_send(temp, npts, MPI_INTEGER, southpe, TAG, pe_grid, ierror)
    endif
  endif
! wait for non blocking receives and plug back
  if(north) then
    if(.not. east) then  ! wait for B from east
      CALL MPI_WAIT(reqe1, state, ierror)
      B = tb
    endif
    if(.not. west) then  ! get A from west
      CALL MPI_WAIT(reqw1, statw, ierror)
      A = ta
    endif
  endif
  if(south) then
    if(.not. east) then  ! wait for C from east
      CALL MPI_WAIT(reqe2, state, ierror)
      C = tc
    endif
    if(.not. west) then  ! wait for D from west
      CALL MPI_WAIT(reqw2, statw, ierror)
      D = td
    endif
  endif
  if(east)  then
    if(.not. north) then  ! wait for B from north
      CALL MPI_WAIT(reqn1, state, ierror)
      B = tb
    endif
    if(.not. south) then  ! wait for C from south
      CALL MPI_WAIT(reqs1, state, ierror)
      C = tc
    endif
  endif
  if(west)  then
    if(.not. north) then  ! wait for A from north
      CALL MPI_WAIT(reqn2, state, ierror)
      A = ta
    endif
    if(.not. south) then  ! wait for D from south
      CALL MPI_WAIT(reqs2, state, ierror)
      D = td
    endif
  endif

  return
end subroutine

subroutine RPN_COMM_propagate_boundary_circular(f,minx,maxx,miny,maxy,lni,lnj,nk,hx,hy)
  use rpn_comm
  implicit none
  integer, intent(IN) :: minx,maxx,miny,maxy,nk   ! dimensions of array f
  integer, dimension(minx:maxx,miny:maxy,nk), intent(INOUT) :: f
  integer, intent(IN) :: lni, lnj                 ! number of private points in array f
  integer, intent(IN) :: hx, hy                   ! useful horizontal halo around f

  integer :: northpe, southpe, eastpe, westpe, upstream, downstream
  integer :: ierror
  integer, dimension(MPI_STATUS_SIZE) :: status ! status for sendrecv
  integer, dimension(hx,hy,nk) :: ta, tb, tc, td                    ! recv buffers
  integer, dimension(hx,hy,nk) :: temp                              ! send buffer
  integer :: npts

#define TAG 0

  if( .not. (bnd_north .or. bnd_south .or. bnd_east .or. bnd_west) ) return ! not on boundary, nothing to do
  if(pe_nx == 1 .and. pe_ny ==1) return                     ! one tile only, nothing to do

  if(pe_nx ==1 .or. pe_ny == 1) then  ! pe_nx x 1 or 1 x pe_ny configuration, use alternate version
    call RPN_COMM_propagate_boundary(f,minx,maxx,miny,maxy,lni,lnj,nk,hx,hy)
    return
  endif

  northpe=pe_id(pe_mex,pe_mey+1)
  southpe=pe_id(pe_mex,pe_mey-1)
  eastpe=pe_id(pe_mex+1,pe_mey)
  westpe=pe_id(pe_mex-1,pe_mey)
  npts = hx * hy * nk           ! number of points to send/receive (zones A, B, C, D, Z1 .. Z8)

! clockwise move, followed by counterclockwise move (upstream/downstream defined by clockwise move)
! get from upstream send downstream , then get from downstream send upstream
  if(bnd_north .and. bnd_east) then       ! get A send 4 : get C send 1
    upstream   = westpe
    downstream = southpe
    call MPI_sendrecv(Z4, npts, MPI_INTEGER, downstream, TAG, ta, npts, MPI_INTEGER, upstream  , TAG, pe_grid, status, ierror)
    A = ta
    call MPI_sendrecv(Z1, npts, MPI_INTEGER, upstream  , TAG, tc, npts, MPI_INTEGER, downstream, TAG, pe_grid, status, ierror)
    C = tc
  else if(bnd_north .and. bnd_west) then  ! get D send 2 : get B send 7
    upstream   = southpe
    downstream = eastpe
    call MPI_sendrecv(Z2, npts, MPI_INTEGER, downstream, TAG, td, npts, MPI_INTEGER, upstream  , TAG, pe_grid, status, ierror)
    D = td
    call MPI_sendrecv(Z7, npts, MPI_INTEGER, upstream  , TAG, tb, npts, MPI_INTEGER, downstream, TAG, pe_grid, status, ierror)
    B = tb
  else if(bnd_south .and. bnd_east) then  ! get B send 6 : get D send 3
    upstream   = northpe
    downstream = westpe
    call MPI_sendrecv(Z6, npts, MPI_INTEGER, downstream, TAG, tb, npts, MPI_INTEGER, upstream  , TAG, pe_grid, status, ierror)
    B = tb
    call MPI_sendrecv(Z3, npts, MPI_INTEGER, upstream  , TAG, td, npts, MPI_INTEGER, downstream, TAG, pe_grid, status, ierror)
    D = td
  else if(bnd_south .and. bnd_west) then  ! get C send 8 : get A send 5
    upstream   = eastpe
    downstream = northpe
    call MPI_sendrecv(Z8, npts, MPI_INTEGER, downstream, TAG, tc, npts, MPI_INTEGER, upstream  , TAG, pe_grid, status, ierror)
    C = tc
    call MPI_sendrecv(Z5, npts, MPI_INTEGER, upstream  , TAG, ta, npts, MPI_INTEGER, downstream, TAG, pe_grid, status, ierror)
    A = ta
  else if(bnd_north) then             ! get A send 2 : get B send 1
    upstream   = westpe
    downstream = eastpe
    call MPI_sendrecv(Z2, npts, MPI_INTEGER, downstream, TAG, ta, npts, MPI_INTEGER, upstream  , TAG, pe_grid, status, ierror)
    A = ta
    call MPI_sendrecv(Z1, npts, MPI_INTEGER, upstream  , TAG, tb, npts, MPI_INTEGER, downstream, TAG, pe_grid, status, ierror)
    B = tb
  else if(bnd_south) then             ! get C send 6 : get D send 5
    upstream   = eastpe
    downstream = westpe
    call MPI_sendrecv(Z6, npts, MPI_INTEGER, downstream, TAG, tc, npts, MPI_INTEGER, upstream  , TAG, pe_grid, status, ierror)
    C = tc
    call MPI_sendrecv(Z5, npts, MPI_INTEGER, upstream  , TAG, td, npts, MPI_INTEGER, downstream, TAG, pe_grid, status, ierror)
    D = td
  else if(bnd_east)  then             ! get B send 4 : get C send 3
    upstream   = northpe
    downstream = southpe
    call MPI_sendrecv(Z4, npts, MPI_INTEGER, downstream, TAG, tb, npts, MPI_INTEGER, upstream  , TAG, pe_grid, status, ierror)
    B = tb
    call MPI_sendrecv(Z3, npts, MPI_INTEGER, upstream  , TAG, tc, npts, MPI_INTEGER, downstream, TAG, pe_grid, status, ierror)
    C = tc
  else if(bnd_west)  then             ! get D send 8 : get A send 7
    upstream   = southpe
    downstream = northpe
    call MPI_sendrecv(Z8, npts, MPI_INTEGER, downstream, TAG, td, npts, MPI_INTEGER, upstream  , TAG, pe_grid, status, ierror)
    D = td
    call MPI_sendrecv(Z7, npts, MPI_INTEGER, upstream  , TAG, ta, npts, MPI_INTEGER, downstream, TAG, pe_grid, status, ierror)
    A = ta
  endif

  return
end subroutine
