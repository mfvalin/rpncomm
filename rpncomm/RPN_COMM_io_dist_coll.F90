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
subroutine RPN_COMM_io_dist_coll_test(nparams,params)
  use rpn_comm
  implicit none
#include <RPN_COMM_interfaces_int.inc>
  integer, intent(IN) :: nparams
  integer, intent(IN), dimension(nparams) :: params
  logical :: periodx, periody
  integer :: setno, me_io, n_io
  integer :: i, j, k
  integer, parameter :: gni=8
  integer, parameter :: gnj=5
  integer, parameter :: gnk=3
  integer, parameter :: lnk=5
  integer, dimension(gnk) :: liste_k
  logical, dimension(lnk) :: liste_o
  integer, dimension(gni,gnj,gnk) :: global
  integer, dimension(:,:,:), allocatable :: local
  integer :: lni, lnj
  integer :: mini,maxi,minj,maxj,status
  integer, dimension(pe_nx) :: start_x,count_x
  integer, dimension(pe_ny) :: start_y,count_y
  integer :: i0,in,j0,jn,nerrors,nvalid,expected
!
  periodx = .false.
  periody = .false.
  liste_k = 0
  liste_o = .false.
!
  lni = (gni+pe_nx-1)/pe_nx
  mini = 1-1
  maxi = lni+1
  count_x = lni
  start_x(1) = 1
  do i = 2,pe_nx
    start_x(i) = start_x(i-1) + count_x(i-1)
  enddo
  count_x(pe_nx) = gni + 1 - start_x(pe_nx)
  lni = count_x(pe_mex+1)
  i0 = mini
  if(pe_mex == 0 .and. (.not. periodx)) i0 = 1
  in = maxi
  if(pe_mex == pe_nx-1 .and. (.not. periodx)) in = count_x(pe_nx)
  print *,"start_x =",start_x
  print *,"count_x =",count_x
!
  lnj = (gnj+pe_ny-1)/pe_ny
  minj = 1-1
  maxj = lnj+1
  count_y = lnj
  start_y(1) = 1
  do j = 2,pe_ny
    start_y(j) = start_y(j-1) + count_y(j-1)
  enddo
  count_y(pe_ny) = gnj + 1 - start_y(pe_ny)
  lnj = count_y(pe_mey+1)
  j0 = minj
  if(pe_mey == 0 .and. (.not. periody)) j0 = 1
  jn = maxj
  if(pe_mey == pe_ny-1 .and. (.not. periody)) jn = count_y(pe_ny)
  print *,"start_y =",start_y
  print *,"count_y =",count_y
!
  allocate(local(mini:maxi,minj:maxj,lnk))
  local = 99999
  global = 88888
! create IO PE set
  setno = RPN_COMM_create_io_set(min(pe_nx,pe_ny),0)
!  print *,'params=',params
  print *,'IO PE set created :',setno
  me_io = RPN_COMM_is_io_pe(setno)
  n_io = RPN_COMM_io_pe_size(setno)
  if(me_io .ne. -1) then
    print *,"I am a busy IO pe!",me_io+1,' of',n_io
    do k=1,1               !me_io+1
      liste_k(k)=4 - me_io   !  *n_io
      do j = 1,gnj
      do i = 1,gni
        global(i,j,k) = liste_k(k) + j*10 + i*1000
      enddo
      enddo
    enddo
    print *,"level list =",liste_k
    do k= gnk,1,-1
      if(liste_k(k) <= 0) cycle
      print *,"===== level ==",liste_k(k),"  ====="
      do j=gnj,1,-1
        print 100,j,global(:,j,k)
      enddo
    enddo
  else
    print *,"I am a relaxed  NON-IO pe !"
  endif
print *,'lni,lnj,mini,maxi,minj,maxj',lni,lnj,mini,maxi,minj,maxj
!return
  call RPN_COMM_shuf_dist(setno,  &
                          global,gni,gnj,gnk,  &
                          local,mini,maxi,minj,maxj,lnk,  &
                          liste_k,liste_o,  &
                          start_x,count_x,pe_nx,start_y,count_y,pe_ny,  &
                          periodx,periody,status)
  print *,'liste_o apres=',liste_o
  nerrors = 0
  nvalid = 0
  do k = lnk,1,-1
    if(liste_o(k)) then
      do j = j0,jn
      do i = i0,in
        nvalid = nvalid + 1
        expected = k + (start_y(pe_mey+1)+j-1)*10 + (start_x(pe_mex+1)+i-1)*1000
        if(expected .ne. local(i,j,k)) then
          print *,'i,j,k,expected,local(i,j,k)',i,j,k,expected,local(i,j,k)
          nerrors = nerrors + 1
          if(nerrors>3)goto 666
        endif
      enddo
      enddo
    else
      print *,'no data at level',k
    endif
  enddo
  print *,"nerrors, nvalid=",nerrors,nvalid
  goto 777
666 continue
#if  ! defined(DEPRECATED)
  do k = lnk,1,-1
    if(liste_o(k)) then
      print *,"===== level",k,"  ====="
      do j = maxj,minj,-1
        print 100,j,mini,maxi,local(i0:in,j,k)
      enddo
    else
      print *,'no data at level',k
    endif
  enddo
#endif
!
100 format(I3,20I6.5)
777 continue
  return
end subroutine RPN_COMM_io_dist_coll_test
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
!      nx,ny,start_x,count_x,start_y,count_y  (start_x, start_y : ORIGIN 1)
!      periodx,periody
!
!    even if not used/necessary on a given PE
!      the "global" array must exist ( (1,1,1) array is OK on non IO PEs )
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
  logical, intent(OUT), dimension(lnk)  :: liste_o
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
!print *,"DEBUG: k,i,status,low,high,liste_i(k)=",k,i,status,low,high,liste_i(k)
      if(status == RPN_COMM_ERROR) return
!      liste_o(liste_i(k)) = .true.
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
    integer, intent(IN) :: setno
    integer, intent(IN) :: gni,gnj,gk
    integer, intent(IN), dimension(gni,gnj) :: global
    integer, intent(IN) :: mini,maxi,minj,maxj,lnk
    integer, intent(OUT), dimension(mini:maxi,minj:maxj,lnk) :: local
    logical, intent(OUT), dimension(lnk) :: liste_o
    integer, intent(IN) :: npes
    integer, intent(IN), dimension(npes) :: pe_x, pe_y
    integer, intent(IN), dimension(0:pe_nx-1) :: start_x,count_x
    integer, intent(IN), dimension(0:pe_ny-1) :: start_y,count_y
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

print *,'IN shuffle'
    on_column = .false.    ! precondition for failure in case a hasty exit is needed
    status = RPN_COMM_ERROR
    nullify(fullrow)
    nullify(local_1)
!
    do i = 2 , npes
      if(pe_x(i) <= pe_x(i-1) .or. pe_y(i) <= pe_y(i-1)) then
        print *,"ERROR: PE list is not monotonous and increasing"
        return  ! pe_x and pe_y must be monotonous and increasing
      endif
    enddo
    root = -1 
    kcol = -1
    haloy = 1 - minj    ! halos are implicitly specified by lower bound of x and y dimensions
    halox = 1 - mini
    size2d = (maxj-minj+1) * (maxi-mini+1)    ! size of a 2D local array
    eff_periodx = periodx .and. (halox > 0)   ! global along x perodioc adjustment needed
    lni = count_x(pe_mex)                     ! useful number of points on local tile
    lnj = count_y(pe_mey)
    if(maxi < lni+halox .or. maxj < lnj+haloy) then
      print *,"ERROR: upper bound of array too small to accomodate halo"
      print *,"mini,maxi,lni,minj,maxj,lnj",mini,maxi,lni,minj,maxj,lnj
      return  ! OOPS, upper bound cannot accomodate halo
    endif
    do i = 1 , npes
      if(pe_mex == pe_x(i)) then
        on_column = .true.                       ! there is an IO PE on the column (and only one)
        root = pe_y(i)                           ! y coordinate of IO PE, will be the root for scatterv
        if(root == pe_mey) kcol = gk             ! level of 2D plane that this PE will distribute (it is the root)
        exit
      endif
    enddo
    if(on_column) call mpi_bcast(kcol,1,MPI_INTEGER,root,pe_mycol,ierr)    ! send gk to all PEs on the column
!
!   get list of which PE has which gk along row
!   if column had no IO PE, -1 gets transmitted
!   if IO PE on column has nothing to contribute, 0 should get transmitted
!
    listofk = 0
    call mpi_allgather(kcol,1,MPI_INTEGER,listofk,1,MPI_INTEGER,pe_myrow,ierr)
!print *,"DEBUG: listofk=",listofk
    if(maxval(listofk) <= 0) then   ! no contribution from any IO PE, job id done for this round
      print *,"WARNING: no work to do"
      status = RPN_COMM_OK
      return
    endif
!   maybe we should check instead if liste_o(listofk(i)) is already .true. which would point to a possible duplicate
    do i = 0 , pe_nx-1
       if(listofk(i) > 0) liste_o(listofk(i)) = .false.
    enddo
!
!   first and last PE on column get count + haloy rows, others get count + 2*haloy rows
!   periody condition to be dealt with later
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
      if(minj < 1)    fullrow(:,minj:0) = 0
      if(maxj > lnj)  fullrow(:,lnj+1:maxj) = 0
!print *,'IN shuffle before scatter'
      call mpi_scatterv(global,cy,dy,MPI_INTEGER,   &          ! 
                        fullrow(1,ybase),cy(pe_mey),MPI_INTEGER, &
                        root,pe_mycol,ierr)
!print *,'IN shuffle after scatter'
!print *,"DEBUG: ======= fullrow for level",kcol
!do j=maxj,minj,-1
!  print 100,j,fullrow(:,j)
!enddo
!
!     fullrow now contains what will be redistributed along x, haloy is accounted for
!     we may now process reshaping and halo along x periodicity condition
!
      allocate(local_1(mini:maxi,minj:maxj,pe_nx))    ! reshape for distribution along x
      lni0 = count_x(0)
      do k = 1 , pe_nx
        i0 = mini
        in = maxi
        if(k == 1) i0 = 1                        ! lower bound along x of west PE
        if(k == pe_nx) in = count_x(pe_nx-1)     ! upper boung along x of east PE
        ioff = start_x(k-1) - 1                  ! offset along x in global space for PE no k-1 along x
        do j = minj , maxj
          local_1(  :  ,j,k) = 0
          local_1(i0:in,j,k) = fullrow(i0+ioff:in+ioff,j)
          if(eff_periodx .and. (k == 1)) &
                 local_1(    1-halox:0    ,j,k) = fullrow(gni-halox+1:gni,j)   ! west halo from global east
          if(eff_periodx .and. (k == pe_nx)) &
                 local_1(lni0+1:lni0+halox,j,k) = fullrow(1:halox,j)           ! east halo from global west
        enddo
      enddo
                    ! send sizes and displacements for the final alltoallv
      cxs = size2d  ! a full 2D local slice will be sent to all PEs on row
      dxs(0) = 0
      do i = 1 , pe_nx-1
        dxs(i) = dxs(i-1) + cxs(i-1)
      enddo
!print *,"DEBUG: ======= local_1 for level",kcol
!do k = pe_nx,1,-1
!  print *,"======= PE", k
!  do j = maxj,minj,-1
!    print 100,j,local_1(:,j,k)
!  enddo
!enddo
    else                           ! we are not on a column where there is an IO PE
      allocate(fullrow(1,1))
      allocate(local_1(1,1,1))
      cxs = 0                      ! nothing to send from here, size = displacement = 0
      dxs = 0
    endif
!
!   receive sizes and displacements for the final alltoallv
!
    do i = 0 , pe_nx-1
      if(listofk(i) > 0) then    ! something (2D slice) will be received from column i
        cxr(i) = size2d
        dxr(i) = (listofk(i)-1) * size2d   ! offset into local is listofk(i)-1 2D planes
      else                       ! nothing to receive from column i
        cxr(i) = 0
        dxr(i) = 0
      endif
    enddo
!
!print *,"DEBUG: kcol,cxs,dxs,11111,cxr,dxr,11111,listofk"
!print 100,kcol,cxs,dxs,11111,cxr,dxr,11111,listofk
!do k=lnk,1,-1
!  print *,'=== lv=',k
!  do j=maxj,minj,-1
!    print 100,j,local(:,j,k)
!  enddo
!enddo
    call mpi_alltoallv(local_1, cxs, dxs, MPI_INTEGER,  &
                       local,   cxr, dxr, MPI_INTEGER,  &
                       pe_myrow, ierr)
!
!do k=lnk,1,-1
!  print *,'=== lv=',k
!  do j=maxj,minj,-1
!    print 100,j,local(:,j,k)
!  enddo
!enddo
!print *,"DEBUG: exiting dist_1"
    if(associated(fullrow)) deallocate(fullrow)
    if(associated(local_1)) deallocate(local_1)
    status = RPN_COMM_OK   ! success at last !
    do i = 0 , pe_nx-1     ! mark 2D array at position listofk(i) as received
       if(listofk(i) > 0) liste_o(listofk(i)) = .true.
    enddo
  100 format(I3,20I6.5)
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
