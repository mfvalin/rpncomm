subroutine RPN_COMM_haloflip(g,minx,maxx,miny,maxy,ni,nj,nk,halox,haloy,gni)
  use rpn_comm
  implicit none
  include 'RPN_COMM_interfaces.inc'
  integer, intent(IN) :: minx,maxx,miny,maxy,ni,nj,nk,halox,haloy,gni
  integer, intent(INOUT) :: g(minx:maxx,miny:maxy,nk)

  integer :: i,j,k,ierr,gmin,gmax,i0,in
  integer :: j_src,j_dest
  integer, dimension(minx:maxx,haloy,nk,0:1) :: local
  integer, dimension(0:pe_nx) :: count, depl
  integer :: start, finish, start2, finish2, low, high, ihave
  integer :: send, recv, send_to, send_tag, recv_from, recv_tag, nwds
  integer :: ipe
  integer :: gis, gie, l_offset, g_offset
  integer, dimension(MPI_STATUS_SIZE) :: status

  if(.not. bnd_south .and. .not. bnd_north) return
  if(bnd_south) then
    j_src = 1
    j_dest = 1-haloy
  endif
  if(bnd_north) then
    j_src = nj - haloy + 1
    j_dest = nj + 1
  endif
  ierr = RPN_COMM_limit(pe_mex,pe_nx,1,gni,gmin,gmax,count,depl)
  do k=1,nk
  do j=1,haloy
  do i=minx,maxx
    local(i,j,k,0) = g(i,j_src+haloy-j,k)
  enddo
  enddo
  enddo
!  do k=nk,1,-1
!  do j=haloy,1,-1
!    print 1,local(:,j,k,0)
!  enddo
!  enddo
!  print *,"-------------------------------------------------------"

  start = gmin + minx -1 + gni/2
  finish = start + (maxx - minx)
!  print *,'start=',start,' finish=',finish
  if(start > gni) then
    start = start - gni
    finish = finish - gni
  endif
  if(finish > gni) then
    start2 = 1
    finish2 = finish - gni
    finish = gni
  else
    start2 = 0
    finish2 = 0
  endif
  send = 0;
  recv = 1;
  ihave = pe_mex
  nwds = nk * haloy * (maxx-minx+1)
  g_offset = minx - start
!  print *,'start=',start,' finish=',finish
!  print *,'start2=',start2,' finish2=',finish2
  do ipe = 1, pe_nx
    low = depl(ihave) + 1
    high = low + count(ihave) -1
    l_offset = 1 - low
    gis = max(low,start)
    gie = min(high,finish)
    print *,'low=',low,' high=',high
    if (gis <= gie) then
      print *,'%gis=',gis,gis+g_offset,' gie=',gie,gie+g_offset
      do k=1,nk
      do j=1,haloy
      do i=gis,gie
        g(i+g_offset,j_dest+j-1,k) = local(i+l_offset,j,k,send)
      enddo
      enddo
      enddo
    endif
    if(low <= finish2) then
!      print *,'low=',low,' finish2=',finish2,' high=',high
      gis = low
      gie = min(high,finish2)
      print *,'+gis=',gis,gis+minx+finish-start,' +gie=',gie,gie+minx+finish-start
      do k=1,nk
      do j=1,haloy
      do i=gis,gie
        g(i+minx+finish-start,j_dest+j-1,k) = local(i+l_offset,j,k,send)
      enddo
      enddo
      enddo
    endif
    send_to = pe_id(pe_mex+1,pe_mey)      ! send to east neighbor
    send_tag = pe_me
    recv_from = pe_id(pe_mex-1,pe_mey)     ! recv from west neighbor
    recv_tag = recv_from
!    print *,pe_me,' sending to ',send_to,' receiving from',recv_from
!    print 1,local(:,1,1,send)
    call mpi_sendrecv(local(minx,1,1,send) , nwds, MPI_INTEGER, send_to,   send_tag,   &
                      local(minx,1,1,recv) , nwds, MPI_INTEGER, recv_from, recv_tag,   &
                      pe_defcomm,status,ierr)
    send = 1 - send
    recv = 1 - recv
    ihave = ihave - 1
    if(ihave<0) ihave = ihave + pe_nx
  enddo

  i0 = 1
  in = count(pe_mex)
!  do k=1,nk
!  do j=1,haloy
!  do i=i0,in
!    g(i,j_dest+j-1,k) = local(i,j,k,0)
!  enddo
!  enddo
!  enddo
1 format(20I7.6)
  return
end subroutine RPN_COMM_haloflip
subroutine RPN_COMM_haloflip_test
  use rpn_comm
  implicit none
  include 'RPN_COMM_interfaces.inc'
  integer :: ierr
  integer, parameter :: halox=1, haloy=2, ni=5, nj=7, nk=1
  integer, parameter :: minx=1-halox-1
  integer, parameter :: maxx=ni+halox+1
  integer, parameter :: miny=1-haloy-1
  integer, parameter :: maxy=nj+haloy+1
  integer, dimension(:,:,:), allocatable :: g
  integer :: gni, i, j, k, ipe, gmin, gmax
  integer, dimension(0:100) :: count, depl

  call mpi_init(ierr)
  pe_defcomm = MPI_COMM_WORLD
  call mpi_comm_rank(MPI_COMM_WORLD,pe_mex,ierr)
  pe_me = pe_mex
  pe_mey = 0
  call mpi_comm_size(MPI_COMM_WORLD,pe_nx,ierr)
  pe_ny = 1
  allocate(pe_id(-1:pe_nx,0:0))
  do i = 0,pe_nx-1
    pe_id(i,0) = i
  enddo
  pe_id(-1,0) = pe_nx-1
  pe_id(pe_nx,0) = 0
  allocate(g(minx:maxx,miny:maxy,nk))
  gni = ni*pe_nx
  gni = gni -  mod(gni,2)   ! must be even
  ierr = RPN_COMM_limit(pe_mex,pe_nx,1,gni,gmin,gmax,count,depl)
  if(pe_mex == 0) then
    write(6,1)'gni   =',gni
    write(6,1)'count =',count(0:pe_nx-1)
    write(6,1)'depl  =',depl(0:pe_nx-1)
    call flush(6)
  endif
  g = 0
  do k=1,nk
  do j=1,nj
  do i=1,count(pe_mex)
    g(i,j,k) = 10000*(i+pe_mex*ni) + 100*j + k
  enddo
  enddo
  enddo
  bnd_south = .false.
  bnd_north = .true.
  call RPN_COMM_haloflip(g,minx,maxx,miny,maxy,ni,nj,nk,halox,haloy,gni)
  goto 9
  bnd_south = .true.
  bnd_north = .false.
!  do k=nk,1,-1
!  do j=maxy,miny,-1
!    print 1,g(:,j,k)
!  enddo
!  enddo
!  print *,"========================================================="
  call RPN_COMM_haloflip(g,minx,maxx,miny,maxy,ni,nj,nk,halox,haloy,gni)
9 continue
  do ipe = 0,pe_nx
    if(ipe == pe_mex) then
      write(6,*)'============',pe_mex,'============'
      do k=nk,1,-1
      do j=maxy,miny,-1
        write(6,1) "",g(:,j,k)
      enddo
      enddo
    endif
    call flush(6)
    call mpi_barrier(MPI_COMM_WORLD,ierr)
  enddo
  call mpi_finalize(ierr)
1 format(A,20I7.6)
  return
end subroutine RPN_COMM_haloflip_test
program my_test
  call RPN_COMM_haloflip_test
  stop
end
