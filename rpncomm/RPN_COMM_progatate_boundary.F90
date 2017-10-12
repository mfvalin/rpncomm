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
! map of regions in array f 
! A, B, C, D external corner regions to get from neighbors
! 1, 2, 3, 4, 5, 6, 7, 8 regions from which neighbors get their external corners
!           +---+---+                      +---+---+----(lnj+hy)
!           | A | 1 |                      | 2 | B |
!           +---+---+----------------------+---+---+----(lnj)
!           | 8 |                              | 3 |
!           +---+                              +---+----(lnj-hy+1)
!               |                              |
!               |                              |
!               |                              |
!               |                              |
!               |                              |
!               |                              |
!           +---+                              +---+----(hy)
!           | 7 |                              | 4 |
!           +---+---+----------------------+---+---+----(1)
!           | D | 6 |                      | 5 | C |
!           +---+---+----------------------+---+---+----(1-hy)
!           |   |   |                      |   |   |
!           |   |  (hx)             (lni-hx+1) |   |
!           |  (1)                           (lni) |
!         (1-hx)                                (lni+hx)
subroutine RPN_COMM_progatate_boundary(f,minx,maxx,miny,maxy,lni,lnj,nk,hx,hy)
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

  north=(bnd_north)
  northpe=pe_id(pe_mex,pe_mey+1)
  south=(bnd_south)
  southpe=pe_id(pe_mex,pe_mey-1)
  east=(bnd_east)
  eastpe=pe_id(pe_mex+1,pe_mey)
  west=(bnd_west)
  westpe=pe_id(pe_mex-1,pe_mey)
  npts = hx * hy * nk

! post appropriate nonblocking receives
  if(north) then
    if(.not. east) then  ! get B from east
      call irecv(tb, npts, MPI_INTEGER, eastpe, TAG, pe_grid, reqe1, ierror)
    endif
    if(.not. west) then  ! get A from west
      call irecv(ta, npts, MPI_INTEGER, westpe, TAG, pe_grid, reqw1, ierror)
    endif
  endif
  if(south) then
    if(.not. east) then  ! get C from east
      call irecv(tc, npts, MPI_INTEGER, eastpe, TAG, pe_grid, reqe2, ierror)
    endif
    if(.not. west) then  ! get D from west
      call irecv(td, npts, MPI_INTEGER, westpe, TAG, pe_grid, reqw2, ierror)
    endif
  endif
  if(east)  then
    if(.not. north) then  ! get B from north
      call irecv(tb, npts, MPI_INTEGER, northpe, TAG, pe_grid, reqn1, ierror)
    endif
    if(.not. south) then  ! get C from south
      call irecv(tc, npts, MPI_INTEGER, southpe, TAG, pe_grid, reqs1, ierror)
    endif
  endif
  if(west)  then
    if(.not. north) then  ! get A from north
      call irecv(ta, npts, MPI_INTEGER, northpe, TAG, pe_grid, reqn2, ierror)
    endif
    if(.not. south) then  ! get D from south
      call irecv(td, npts, MPI_INTEGER, southpe, TAG, pe_grid, reqs2, ierror)
    endif
  endif
! extract and blocking send to appropriate destination
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
      D = tc
    endif
  endif
  if(west)  then
    if(.not. north) then  ! wait for B from north
      CALL MPI_WAIT(reqn1, state, ierror)
      B = tb
    endif
    if(.not. south) then  ! wait for C from south
      CALL MPI_WAIT(reqs1, state, ierror)
      C = tc
    endif
  endif
  if(east)  then
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
