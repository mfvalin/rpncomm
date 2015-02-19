!/! RMNLIB - Library of useful routines for C and FORTRAN programming
! ! Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
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
! distribute gnk 2D arrays, destination is gnk of lnk 2D plane of local 3D array
! using IO PE set setno
! liste_i contains the list of 2D planes to be inserted into the 3D destination array
! liste_o(k) is set to .true. when 2D plane k has been received
! start_x(i) contains the first global x index of local array in grid PE (i-1,any)
! count_x(i) contains the number of points along x of local array in grid PE (i-1,any)
! nx dimension of start_x and count_x (ERROR if not equal to pe_nx)
! start_y(j) contains the first global y index of local array in grid PE (any,j-1)
! count_y(j) contains the number of points along y of local array in grid (any,j-1)
! ny dimension of start_y and count_y (ERROR if not equal to pe_ny)
! mini,maxi,minj,maxj are used to determine the halo width in the local array
!                     this is used to save a halo exchange after distribution
!
! start_x and start_y are in origin(1)
!
subroutine RPN_COMM_shuf_dist(setno,  &
                              global,gni,gnj,gnk,  &
                              local,mini,maxi,minj,maxj,lnk,  &
                              liste_i,liste_o,  &
                              start_x,count_x,nx,start_y,count_y,ny,  &
                              periodx,periody,status)
  use RPN_COMM_io_pe_tables
  implicit none
  integer, intent(IN) :: setno,gni,gnj,gnk,mini,maxi,minj,maxj,lnk,nx,ny
  integer, intent(IN), dimension(gni,gnj,gnk) :: global
  integer, intent(OUT), dimension(mini:maxi,minj:maxj,lnk) :: local
  integer, intent(IN), dimension(nx)    :: start_x,count_x
  integer, intent(IN), dimension(ny)    :: start_y,count_y
  integer, intent(IN), dimension(gnk)   :: liste_i
  integer, intent(OUT), dimension(lnk)  :: liste_o
  logical, intent(IN) :: periodx,periody
  integer, intent(OUT) :: status
  integer :: i, k, low, high

  status = RPN_COMM_ERROR
  if(setno < 1 .or. setno > iosets) return    ! setno out of bounds
  if(io_set(setno)%ioset .ne. setno) return   ! set no longer valid
  if(pe_nx .ne. nx .or. pe_ny .ne. ny) return ! wrong number of PEs in grid
  if(periody) return                          ! not supported for the time being
!
! distribute one 2D plane at a time, one group of IO PEs at a time
! in a group of IO PEs, no column has more than 1 IO PE, neither has any row
!
  do k= 1 , gnk                         ! loop over 2D planes to distribute
    do i= 1 , io_set(setno)%ngroups     ! loop over groups in this set of IO PEs
      low = 1 + (i-1) * i               ! index of first PE in goup
      high = min( io_set(setno)%npe , low+io_set(setno)%groupsize-1 )  ! index of last PE in group
      call RPN_COMM_shuf_dist_1(setno, &
                                global(1,1,k),gni,gnj,liste_i(k),  &
                                local,mini,maxi,minj,maxj,lnk,  &
                                liste_o, io_set(setno)%x(low:high),  io_set(setno)%y(low:high), (high-low+1), &
                                start_x,count_x,start_y,count_y,  &
                                periodx,periody,status)
      if(status == RPN_COMM_ERROR) return
      liste_o(liste_i(k)) = .true.
    enddo
  enddo
  contains
!
! distribute one 2D array at a time
! no column has 2 IO PES, neither has any row
!
  subroutine RPN_COMM_shuf_dist_1(setno,  &
                                  global,gni,gnj,gk,  &
                                  local,mini,maxi,minj,maxj,lnk,  &
                                  liste_o, pe_x, pe_y, npes, &
                                  start_x,count_x,start_y,count_y,  &
                                  periodx,periody,status)
    use rpn_comm
    use RPN_COMM_io_pe_tables
    implicit none
    integer, intent(IN) :: setno,gni,gnj,gk,mini,maxi,minj,maxj,lnk
    integer, intent(IN), dimension(gni,gnj) :: global
    integer, intent(OUT), dimension(mini:maxi,minj:maxj,lnk) :: local
    integer, intent(IN), dimension(0:pe_nx-1) :: start_x,count_x
    integer, intent(IN), dimension(0:pe_ny-1) :: start_y,count_y
    integer, intent(OUT), dimension(lnk) :: liste_o
    integer, intent(IN), dimension(npes) :: pe_x, pe_y
    integer, intent(IN) :: npes
    logical, intent(IN) :: periodx,periody
    integer, intent(OUT) :: status
    logical :: on_column
    integer :: root, ybase, kcol, ierr, halox, haloy, lni, lnj
    integer, dimension(0:pe_nx-1) :: listofk
    integer :: i, j
    integer, dimension(:,:), pointer :: fullrow
    integer, dimension(:,:,:), pointer :: local_1
    integer, dimension(0:pe_nx-1) :: cx, dx
    integer, dimension(0:pe_ny-1) :: cy, dy

    liste_o(gk) = .false.
    on_column = .false.
    root = -1 
    kcol = -1
    haloy = minj + 1
    halox = mini + 1
    lni = count_x(pe_mex)
    lnj = count_y(pe_mey)
    do i = 1 , npes
      if(pe_mex == pe_x(i)) then
        on_column = .true.
        root = pe_y(i)
        if(root == pe_mey) kcol = gk             ! level of 2D plane for this column
        exit
      endif
    enddo
!
!   first and last PE along column get count + haloy, others get count + 2*haloy
!
    cy = count_y + haloy
    cy(1:pe_ny-2) = cy(1:pe_ny-2) + haloy       ! count, adjusted for halo (periody assumed false)
    cy = cy * gni                               ! and multiply by row length
!
!   all except first PE along column have their starting point bumped down by one haloy width
!
    dy = start_y - 1
    dy(1:pe_ny-1) = dy(1:pe_ny-1) - haloy       ! displacement, adjusted for halo (periody assumed false)
    dy = dy * gni                               ! and multiply by row length
    ybase = minj
    if(pe_mey == 0) ybase = 1                    ! south PE gets no south halo (periody assumed false)
!
    if(on_column)then  ! i am on a column where there is a member of the set
      allocate(fullrow(gni,minj:maxj))
      call mpi_bcast(kcol,1,MPI_INTEGER,root,pe_mycol,ierr)
      call mpi_scatterv(global,cy,dy,MPI_INTEGER,   &
                        fullrow(1,ybase),cy(pe_mey),MPI_INTEGER, &
                        root,pe_mycol,ierr)
    else
      allocate(fullrow(1,1))
    endif
!
!   fullrow now contains what will be redistributed along x, haloy is accounted for
!   we may now process halo along x periodicity condition
!
    if(periodx .and. halox > 0)then
      do j = minj , maxj
        fullrow(mini:0,j)          = fullrow(lni-halox+1:lni,j)   ! east halo
        fullrow(lni+1:lni+halox,j) = fullrow(1:halox,j)           ! west halo
      enddo
    endif
!
!   get list of which PE has which gk along row
!   if column had no IO PE, -1 gets transmitted
!   if IO PE on column had nothing to contribute, 0 should get transmitted
!
    call mpi_allgather(kcol,1,MPI_INTEGER,listofk,1,MPI_INTEGER,pe_myrow,ierr)
!
!   we can now compute sizes and displacements for the final alltoallv
!
!   final alltoallv
!
    local = 0
    status = RPN_COMM_OK
  end subroutine RPN_COMM_shuf_dist_1
end subroutine RPN_COMM_shuf_dist
!
subroutine RPN_COMM_shuf_coll(setno,status)
  use rpn_comm
  use RPN_COMM_io_pe_tables
  implicit none
  integer, intent(IN) :: setno
  integer, intent(OUT) :: status
!
  status = RPN_COMM_ERROR
  if(setno < 1 .or. setno > iosets) return    ! setno out of bounds
  if(io_set(setno)%ioset .ne. setno) return   ! set no longer valid
!
!
! collect one 2D plane at a time
! no column has 2 IO PES, neither has any row
!
  contains
!
  subroutine RPN_COMM_shuf_coll_1(setno,status)
    use rpn_comm
    use RPN_COMM_io_pe_tables
    implicit none
    integer, intent(IN) :: setno
    integer, intent(OUT) :: status
!
    status = RPN_COMM_OK
  end subroutine RPN_COMM_shuf_coll_1
end subroutine RPN_COMM_shuf_coll
!
