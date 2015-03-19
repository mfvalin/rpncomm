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
  integer, parameter :: gni=10
  integer, parameter :: gnj=7
  integer, parameter :: gnk=2
  integer, parameter :: lnk=6
  integer, parameter :: iope_extra=3
  integer, dimension(gnk) :: liste_k, liste_k2
  logical, dimension(lnk) :: liste_o
  integer, dimension(gni,gnj,gnk) :: global,global2
  integer, dimension(:,:,:), allocatable :: local
  integer :: lni, lnj
  integer :: mini,maxi,minj,maxj,status
  integer, dimension(pe_nx) :: start_x,count_x
  integer, dimension(pe_ny) :: start_y,count_y
  integer :: i0,in,j0,jn,nerrors,nvalid,expected,effective_lnk
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
  if(pe_me == 0) then
    print *,"start_x =",start_x
    print *,"count_x =",count_x
  endif
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
  if(pe_me == 0) then
    print *,"start_y =",start_y
    print *,"count_y =",count_y
  endif
!
  allocate(local(mini:maxi,minj:maxj,lnk))
  local = 99999
  global = 88888
! create IO PE set
  setno = RPN_COMM_create_io_set( min( min(pe_nx,pe_ny)+iope_extra , lnk) ,0)  ! make sure not to overflow lnk
!  print *,'params=',params
  print *,'IO PE set created :',setno
  me_io = RPN_COMM_is_io_pe(setno)
  n_io = RPN_COMM_io_pe_size(setno)
  if(me_io .ne. -1) then
    print *,"I am a busy IO pe!",me_io+1,' of',n_io
    do k=1,1               !me_io+1
      liste_k(k) = lnk - me_io   !  *n_io
!      liste_k(k) = 1 + me_io   !  *n_io
      do j = 1,gnj
      do i = 1,gni
        global(i,j,k) = liste_k(k) + j*10 + i*1000
      enddo
      enddo
    enddo
    print *,"level list =",liste_k
    do k= gnk,1,-1
      if(liste_k(k) <= 0) cycle
      print *,"===== source level ==",liste_k(k),"  ====="
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
  print *,"level ,nerrors, nvalid=",k,nerrors,nvalid
  goto 777
666 continue
#if  ! defined(DEPRECATED)
  do k = lnk,1,-1
    if(liste_o(k)) then
      print *,"===== level",k," local ====="
      do j = maxj,minj,-1
!        print 100,j,local(i0:in,j,k)
        print 100,j,local(:,j,k)
      enddo
    else
      print *,'no data at level',k
    endif
  enddo
#endif
!
100 format(I3,20I6.5)
777 continue
!
! collect test follows distribute test
!
  global2 = 99999
  liste_k2 = -99999
  effective_lnk = 1
  do k = 1 , lnk
    if(liste_o(k)) effective_lnk = k  ! highest level available
  enddo
!  print *,"====== calling  shuf_coll ======",effective_lnk
  call RPN_COMM_shuf_coll(setno,  &
                          global2,gni,gnj,gnk,  &
                          local,mini,maxi,minj,maxj,effective_lnk,  &
                          liste_k2,  &
                          start_x,count_x,pe_nx,start_y,count_y,pe_ny,  &
                          status)
! global2 should be identical to global once k has been adjusted
! expected k + global - mod(global,10)
  if(me_io .ne. -1) then ! I am an IO PE
    print *,"====== after shuf_coll ======"
    print *,"DEBUG: pe_me, liste_k2", pe_me,liste_k2
    nerrors = 0
    nvalid = 0
    do k = 1,gnk
      if(liste_k2(k) <= 0) cycle
      do j = 1,gnj
      do i = 1,gni
        expected = global(i,j,k) - mod(global(i,j,k),10) + liste_k2(k)
        nvalid = nvalid + 1
        if(expected .ne. global2(i,j,k)) nerrors = nerrors + 1
      enddo
      enddo
    enddo
    print *,"number of errors =",nerrors
    print *,"number of points =",nvalid
    do k = 1, gnk
      if(liste_k2(k) > 0) then
        print *,"===== k, collected level",k,liste_k2(k),"  ====="
        do j = gnj , 1 , -1
          print 100,j,global2(:,j,k)
        enddo
      endif
    enddo
  endif
  return
end subroutine RPN_COMM_io_dist_coll_test
!
! distribute gnk 2D arrays, destination is gnk of lnk 2D plane of local 3D array
! using IO PE set setno
! liste_i contains the list of 2D planes to be inserted into the 3D destination array
! liste_o(k) is set to .true. when 2D plane k has been received
! start_x(i) contains the first global x index of local array in grid PE (i-1,any)
! count_x(i) contains the number of points along x of local array in grid PE (i-1,any)
! nx is the dimension of start_x and count_x (ERROR if not equal to pe_nx)
! start_y(j) contains the first global y index of local array in grid PE (any,j-1)
! count_y(j) contains the number of points along y of local array in grid (any,j-1)
! ny is the dimension of start_y and count_y (ERROR if not equal to pe_ny)
! mini,maxi,minj,maxj are used to determine the halo width in the local array
!                     this is used to save a halo exchange after distribution
!
! start_x, start_y, count_x, count_y  are in origin(1)
!
! important notes:
!    the following parameters MUST de the same on ALL PEs
!      setno
!      mini,maxi,minj,maxj,lnk  (dimensions of the "local" array)
!      gni,gnj,gnk   (dimensions of the "global" array)
!      nx,ny,start_x,count_x,start_y,count_y  (start_x, start_y : ORIGIN 1)
!      periodx,periody
!      liste_o
!
!    even if not used/necessary on a given PE
!      the "global" array must exist ( (1,1,1) array is OK on non IO PEs )
!      array liste_i must exist  ( a dimension(1) array is OK on non IO PEs )
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
  integer :: i, k, low, high, listeik, n

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
      low = 1 + (i-1) * io_set(setno)%groupsize      ! index of first PE in goup
      high = min( io_set(setno)%npe , low+io_set(setno)%groupsize-1 )  ! index of last PE in group
!print *,"DEBUG: group=",i," of",io_set(setno)%ngroups
!print *,"DEBUG: , PEs x",io_set(setno)%x(low:high)
!print *,"DEBUG: , PEs y",io_set(setno)%y(low:high)
      listeik = 0
      do n = low, high  ! is this PE part of this group ?
        if(pe_mex == io_set(setno)%x(n) .and. pe_mey == io_set(setno)%y(n)) listeik = liste_i(k)
      enddo
!      print *,"DEBUG: shuf_dist, k,i=",k,i
      call RPN_COMM_shuf_dist_1(setno, &
                                global(1,1,k),gni,gnj,listeik,  &
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

!print *,'IN shuffle'
    on_column = .false.    ! precondition for failure in case a hasty exit is needed
    status = RPN_COMM_ERROR
    nullify(fullrow)
    nullify(local_1)
!
!    do i = 2 , npes   ! assumption to be revised, may not stay true with PE dispersion
!      if(pe_x(i) <= pe_x(i-1) .or. pe_y(i) <= pe_y(i-1)) then
!        print *,"WARNING: PE list is not monotonous and increasing"
!        return  ! pe_x and pe_y must be monotonous and increasing
!      endif
!    enddo
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
if(pe_me==0) print *,"DEBUG: kcol,listofk", kcol,listofk
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
!====================================================================================
!
! collect nk 2D array sections into IO PEs from the grid PEs using IO PE set setno
!
! nk 2D array sections coming in, reassembled on the IO PEs of set setno
! each IO PE ends up with a small number (nk / size of IO PE set) of full arrays
!
! liste_o(n) is set to k when 2D array level k has been received
! start_x(i) contains the first global x index of local array in grid PE (i-1,any)
! count_x(i) contains the number of points along x of local array in grid PE (i-1,any)
! nx is the dimension of start_x and count_x (ERROR if not equal to pe_nx)
! start_y(j) contains the first global y index of local array in grid PE (any,j-1)
! count_y(j) contains the number of points along y of local array in grid (any,j-1)
! ny is the dimension of start_y and count_y (ERROR if not equal to pe_ny)
! mini,maxi,minj,maxj are used to determine the halo width in the local array
!                     this is used to save a halo exchange after distribution
!
! start_x, start_y, count_x, count_y  are in origin(1)
!
! important notes:
!    the following parameters MUST de the same on ALL PEs
!      setno
!      mini,maxi,minj,maxj,nk  (dimensions of the "local" array)
!      gni,gnj,gnk   (dimensions of the "global" array)
!      nx,ny,start_x,count_x,start_y,count_y  (start_x, start_y : ORIGIN 1)
!
!    even if not used/necessary on a given PE
!      the "global" array must exist ( (1,1,1) array is OK on non IO PEs )
!      array liste_o must exist  ( a dimension(1) array is OK on non IO PEs )
!
!====================================================================================
subroutine RPN_COMM_shuf_coll(setno,  &
                              global,gni,gnj,gnk,  &
                              local,mini,maxi,minj,maxj,nk,  &
                              liste_o,  &
                              start_x,count_x,nx,start_y,count_y,ny,  &
                              status)
  use rpn_comm
  use RPN_COMM_io_pe_tables
  implicit none
  integer, intent(IN) :: setno,gni,gnj,gnk,mini,maxi,minj,maxj,nk,nx,ny
  integer, intent(OUT), dimension(gni,gnj,gnk) :: global
  integer, intent(IN), dimension(mini:maxi,minj:maxj,nk) :: local
  integer, intent(IN), dimension(nx)    :: start_x,count_x
  integer, intent(IN), dimension(ny)    :: start_y,count_y
  integer, intent(OUT), dimension(gnk)  :: liste_o
  integer, intent(OUT) :: status
  integer :: iset, npass, setsize, igroup, groupsize
  integer k0, k1, kn, low, high
!
  status = RPN_COMM_ERROR
  if(setno < 1 .or. setno > iosets) return    ! setno out of bounds
  if(io_set(setno)%ioset .ne. setno) return   ! set no longer valid
  if(pe_nx .ne. nx .or. pe_ny .ne. ny) return ! wrong number of PEs in grid
!
  if(io_set(setno)%me >= 0) liste_o = 0      ! this PE is a member of the IO set
  setsize = io_set(setno)%npe
  groupsize = io_set(setno)%groupsize
  npass = (nk+setsize-1) / setsize
  if(gnk < npass) then   ! OUCH, cannot collect
    print *,"ERROR: cannot collect. setsize, groupsize,npass,gnk,nk",setsize, groupsize,npass,gnk,nk
    return
  endif
  k0 = 1
!
! loop over the number of passes necessary to distribute nk over setsize IO PEs
! each IO PE will receive (npass) or (npass - 1) full arrays
!
  do iset = 1 , npass   ! process up to setsize levels per iteration
    k1 = k0             ! base level for the first group of this pass
!
!   loop over the groups in the IO PE set
!   the active routine needs the "no column has 2 IO PES, neither has any row" condition
!   IO PEs low -> high in set will potentially receive a full 2D array
!   from levels k1 -> kn   ( levels k1 -> min(nk , k1+groupsize-1) )
!   if more IO PEs in group than levels to distribute some IO PEs will receive nothing
!
    do igroup = 1 , io_set(setno)%ngroups  ! loop over groups in IO PE set
      kn = min( nk , k1+groupsize-1)
      low = 1 + (igroup-1) * io_set(setno)%groupsize      ! index of first IO PE in goup
      high = min( io_set(setno)%npe , low+groupsize-1 )   ! index of last PE in group
#if defined(FULL_DEBUG)
      print *,"DEBUG: shuf_coll, iset, igroup =",iset, igroup
      print *,"DEBUG: shuf_coll, low, high =",low, high
      print *,"DEBUG: shuf_coll, nk, k1, kn =",nk, k1, kn
#endif
      call RPN_COMM_shuf_coll_1(setno,  &
                                global(1,1,iset),gni,gnj,  &
                                local,mini,maxi,minj,maxj,nk,k1,kn,  &
                                liste_o(iset), io_set(setno)%x(low:high), io_set(setno)%y(low:high), (high-low+1), &
                                start_x,count_x,start_y,count_y,  &
                                status)
      k1 = k1 + groupsize  ! base level for next group
    enddo
    k0 = k0 + setsize      ! base level for the first group of next pass
  enddo
!
! collect one 2D plane at a time
! no column has 2 IO PES, neither has any row
!
  contains
!
  subroutine RPN_COMM_shuf_coll_1(setno,  &
                                  global, gni, gnj, &
                                  local, mini, maxi, minj, maxj, nk, k1, kn, &
                                  levnk, pe_x, pe_y, npes, &
                                  start_x, count_x, start_y, count_y, &
                                  status)
!
! collect one 2D array at a time from a set of 2D array sections
! if kn-k1+1 smaller than groupsize, some IO PEs will not receive anything
! the Mth PE will receive level k1+M-1
! assumption: no column has 2 IO PES, neither has any row
!
!
    use rpn_comm
    use RPN_COMM_io_pe_tables
    implicit none
    integer, intent(IN) :: setno,gni,gnj,mini,maxi,minj,maxj,nk,k1,kn
    integer, intent(OUT), dimension(gni,gnj) :: global
    integer, intent(IN), dimension(mini:maxi,minj:maxj,nk) :: local
    integer, intent(IN) :: npes
    integer, intent(IN), dimension(npes)  :: pe_x, pe_y
    integer, intent(IN), dimension(0:pe_nx-1)    :: start_x,count_x
    integer, intent(IN), dimension(0:pe_ny-1)    :: start_y,count_y
    integer, intent(OUT)  :: levnk
    integer, intent(OUT) :: status
    integer :: i, j, kexpected
    integer, dimension(0:pe_nx-1) :: cxs, dxs, cxr, dxr
    integer, dimension(0:pe_ny-1) :: cy, dy
    integer, dimension(:,:,:), allocatable :: local_1
    integer, dimension(:,:), allocatable :: local_2
    logical :: io_on_column
    integer :: blocki, blockj, k, nlev, klev, ierr, column_root, dimenj
!
    status = RPN_COMM_ERROR
    nlev = kn - k1 + 1             ! number of levels to distribute
    if(nlev > npes) then           ! nlev cannot be larger than npes
      ! add error message
      return
    endif
    kexpected = 0
    do i = 1 , npes  ! if this PE is part of the group, flag level
      if( pe_mex == pe_x(i) .and. pe_mey == pe_y(i) ) then
        kexpected = k1 + i - 1   ! expected level on this PE
        if(kexpected <= kn)   levnk = -kexpected    ! if successful, levnk will be +kexpected
      endif
    enddo
    io_on_column = .false.
    column_root = -1
    klev = 0
    do i = 1 , nlev
      if(pe_x(i) == pe_mex) then
        io_on_column = .true.         ! there a useful IO PE on this column
        klev = k1 + i -1              ! level klev will be collected on this column
        column_root = pe_y(i)         ! by PE having pe_mey == column_root 
      endif
    enddo
#if defined(FULL_DEBUG)
      print *,"DEBUG: shuf_coll1, io_on_column, nk, k1, kn =", io_on_column, nk, k1, kn
      print *,"DEBUG: shuf_coll1, kexpected, klev, column_root =",kexpected, klev, column_root
#endif
!
!   PASS 1 , row alltoallv to get a local (gni,lnj) array with the proper level
!            on PEs in the same column as IO PEs
!            PEs where pe_mex = pe_x(l) receive level (k1+l-1) into (mini:maxi,lnj,pe_nx) array
!            PEs on row will send local(:,1:lnj,k1+l-1) to pe_x(l)
!
    blocki = maxi - mini + 1       ! number of points along x in array local
    blockj = count_y(pe_mey)       ! number of useful points along y in row pe_mey
    dimenj = maxj - minj + 1
!
    if(io_on_column) then               ! this PE will be collection one level for this row
      allocate( local_1(mini:maxi,blockj,0:pe_nx-1) )
      local_1 = 77777
      allocate( local_2(gni,blockj) )
    else
      allocate( local_1(1,1,1) )
      allocate( local_2(1,1) )
    endif
!
    cxs = 0
    dxs = 0
    do i = 1 , nlev
      j = pe_x(i)                ! PE j will receive level k1 + (i-1)
      cxs(j) =  blocki * blockj
      dxs(j) = (blocki * dimenj) * (k1 + i - 2)  ! offset of level  k1 + (i-1)
      dxs(j) = dxs(j) + blocki * (1-minj)        ! correct for halo along y
    enddo
!
    cxr = 0
    dxr = 0
    if(klev > 0) then          ! this PE collects level klev from all PEs in row
      cxr = blocki * blockj    ! x y 2D block size
      do i = 0 , pe_nx -1
        dxr(i) = (blocki * blockj) * i
      enddo
    endif
!
#if defined(FULL_DEBUG)
    print *,"DEBUG: ========= local ============="
    do k = kn , k1 , -1
      do j = maxj , minj , -1
        print 101,k,j,local(:,j,k)
      enddo
    enddo
    print *,"DEBUG: blocki, blockj",blocki, blockj
    print *,"DEBUG: cxs",cxs
    print *,"DEBUG: dxs",dxs
    print *,"DEBUG: cxr",cxr
    print *,"DEBUG: dxr",dxr
#endif
    call mpi_alltoallv(local,   cxs, dxs, MPI_INTEGER,  &   ! send from local, size cxs, displacement dxs
                       local_1, cxr, dxr, MPI_INTEGER,  &   ! receive into local_1, 
                       pe_myrow, ierr)
!
!   REORG    move from (mini:maxi,lnj,pe_nx) to (gni,lnj)
    if(io_on_column) then
#if defined(FULL_DEBUG)
      print *,"DEBUG: ======== local_1 ========="
#endif
      do k = 0 , pe_nx-1
#if defined(FULL_DEBUG)
        do j = blockj , 1 , -1
          print 101,k,j,local_1(:,j,k)
        enddo
#endif
        do j = 1 , blockj
          local_2( start_x(k):start_x(k)+count_x(k)-1 , j ) = local_1( 1:count_x(k) , j , k )
        enddo
      enddo
#if defined(FULL_DEBUG)
      print *,"DEBUG: ======== local_2 ========="
      do j = blockj , 1 , -1
        print 100,j,local_2(:,j)
      enddo
#endif
    endif
    deallocate( local_1 )   ! local_1 no longer useful
!
!   PASS 2 , on columns where there is an IO PE, gatherv on IO PE into (gni,gnj) array
!
    if(io_on_column) then   ! column root will collect level  klev from local_2 into global
       cy = 0
       dy = 0
       if(column_root == pe_mey) then  ! this PE is the gather root
         cy = gni * count_y
         dy = gni * (start_y -1)
       endif
#if defined(FULL_DEBUG)
       print *,"DEBUG: ======== before gatherv on",column_root," ========="
       print *,"DEBUG: cy",cy
       print *,"DEBUG: dy",dy
#endif
       call mpi_gatherv(local_2, gni*blockj, MPI_INTEGER, &
                        global , cy,  dy,  MPI_INTEGER, &
                        column_root, pe_mycol, ierr)
       if(kexpected .ne. 0) levnk = kexpected
    endif
    deallocate( local_2 )
    status = RPN_COMM_OK
100 format(I3,20I6.5)
101 format(2I3,20I6.5)
  end subroutine RPN_COMM_shuf_coll_1
end subroutine RPN_COMM_shuf_coll
!
