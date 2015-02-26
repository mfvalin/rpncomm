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
! important notes:
!    the following parameters MUST de the same on ALL PEs
!      setno
!      mini,maxi,minj,maxj,lnk  (dimensions of the "local" array)
!      gni,gnj,gnk   (dimensions of the "global" array)
!      nx,ny,start_x,count_x,start_y,count_y
!      periodx,periody
!
!    even if not used/necessary on a given PE
!      the "global" array must exist ( (1,1,1) dimension array is OK )
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
! assumption: no column has 2 IO PES, neither has any row
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
    logical :: on_column, error, eff_periodx
    integer :: root, ybase, kcol, ierr, halox, haloy, lni, lnj, size2d
    integer, dimension(0:pe_nx-1) :: listofk
    integer :: i, j, k, i0, in, ioff, lni0
    integer, dimension(:,:), pointer :: fullrow
    integer, dimension(:,:,:), pointer :: local_1
    integer, dimension(0:pe_nx-1) :: cxs, dxs, cxr, dxr
    integer, dimension(0:pe_ny-1) :: cy, dy

    on_column = .false.    ! precondition for failure in case a hasty exit is needed
    status = RPN_COMM_ERROR
    nullify(fullrow)
    nullify(local_1)
!
    do i = 2 , npes
      if(pe_x(i) <= pe_x(i-1) .or. pe_y(i) <= pe_y(i-1)) return  ! pe_x and pe_y must be monotonous and increasing
    enddo
    root = -1 
    kcol = -1
    haloy = 1 - minj    ! halos are implicitly specified by lower bound of x and y dimensions
    halox = 1 - mini
    size2d = (maxj-minj+1) * (maxi-mini+1)
    eff_periodx = periodx .and. (halox > 0)
    lni = count_x(pe_mex)
    lnj = count_y(pe_mey)
    if(maxi < lni+halox .or. maxj < lnj+haloy) return  ! upper boun cannot accomodate halo
    do i = 1 , npes
      if(pe_mex == pe_x(i)) then
        on_column = .true.                       ! there is an IO PE on the column
        root = pe_y(i)                           ! y coordinate of IO PE, will be the root for scatterv
        if(root == pe_mey) kcol = gk             ! level of 2D plane that this PE will distribute (it is the root)
        exit
      endif
    enddo
    if(on_column) call mpi_bcast(kcol,1,MPI_INTEGER,root,pe_mycol,ierr)    ! send gk to all PEs on the column
!
!   get list of which PE has which gk along row
!   if column had no IO PE, -1 gets transmitted
!   if IO PE on column had nothing to contribute, 0 should get transmitted
!
    listofk = 0
    call mpi_allgather(kcol,1,MPI_INTEGER,listofk,1,MPI_INTEGER,pe_myrow,ierr)
    if(maxval(listofk) <= 0) then   ! no contribution from any IO PE
      status = RPN_COMM_OK
      return
    endif
    do i = 0 , pe_nx-1
       if(listofk(i) > 0) liste_o(listofk(i)) = .false.
    enddo
!
!   first and last PE on column get count + haloy rows, others get count + 2*haloy rows
!
    cy = count_y + haloy
    cy(1:pe_ny-2) = cy(1:pe_ny-2) + haloy       ! count, adjusted for halo (periody assumed false)
    cy = cy * gni                               ! and multiply by row length
!
!   all except first PE on column have their starting point bumped down by one haloy width
!
    dy = start_y - 1
    dy(1:pe_ny-1) = dy(1:pe_ny-1) - haloy       ! displacement, adjusted for halo (periody assumed false)
    dy = dy * gni                               ! and multiply by row length
    ybase = minj
    if(pe_mey == 0) ybase = 1                   ! south PE gets no south halo (periody assumed false)
!
    if(on_column .and. kcol > 0)then    ! this PE is on a column where a member of the IO PE group has something to send
      allocate(fullrow(gni,minj:maxj))
      allocate(local_1(mini:maxi,minj:maxj,pe_nx))
      call mpi_scatterv(global,cy,dy,MPI_INTEGER,   &          ! 
                        fullrow(1,ybase),cy(pe_mey),MPI_INTEGER, &
                        root,pe_mycol,ierr)
!
!     fullrow now contains what will be redistributed along x, haloy is accounted for
!     we may now process reshaping and halo along x periodicity condition
!
      lni0 = count_x(0)
      do k = 1 , pe_nx
        i0 = mini
        in = maxi
        if(k == 1) i0 = 1
        if(k == pe_nx) in = count_x(pe_nx-1)
        ioff = start_x(k-1) - 1
        do j = minj , maxj
          local_1(  :  ,j,k) = 0
          local_1(i0:in,j,k) = fullrow(i0+ioff:in+ioff,j)
          if(eff_periodx .and. (k == 1)) &
                 local_1(    1-halox:0    ,j,k) = fullrow(gni-halox+1:gni,j)   ! west halo from global east
          if(eff_periodx .and. (k == pe_nx)) &
                 local_1(lni0+1:lni0+halox,j,k) = fullrow(1:halox,j)           ! east halo from global west
        enddo
      enddo
      cxs = size2d  ! send sizes and displacements for the final alltoallv
      dxs(0) = 0
      do i = 1 , pe_nx-1
        dxs(i) = dxs(i-1) + cxs(i-1)
      enddo
    else
      allocate(fullrow(1,1))
      allocate(local_1(1,1,pe_nx))
      cxs = 0                      ! nothing to send from here, size = displacement = 0
      dxs = 0
    endif
!
!   we can now compute receive sizes and displacements for the final alltoallv
!
    do i = 0 , pe_nx-1
      if(listofk(i) > 0) then
        cxr(i) = size2d
        dxr(i) = (listofk(i)-1) * size2d
      else
        cxr(i) = 0
        dxr(i) = 0
      endif
    enddo
!   final alltoallv
!
    local = 0
!
    if(associated(fullrow)) deallocate(fullrow)
    if(associated(local_1)) deallocate(local_1)
    status = RPN_COMM_OK   ! success at last !
    do i = 0 , pe_nx-1
       if(listofk(i) > 0) liste_o(listofk(i)) = .true.
    enddo
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
