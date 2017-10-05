subroutine RPN_COMM_progatate_boundary(f,minx,maxx,miny,maxy,lni,lnj,nk,hx,hy)
  use rpn_comm
  implicit none
  integer, intent(IN) :: minx,maxx,miny,maxy,lni,lnj,nk,hx,hy
  integer, dimension(minx:maxx,miny:maxy,nk), intent(INOUT) :: f

  logical :: north, south
  integer :: northpe, southpe
  logical :: east, west
  integer :: eastpe, westpe
  integer, dimension(MPI_STATUS_SIZE) :: status1, status2     ! irecv statuses
  integer, dimension(hx,hy,nk) :: put1, put2, get1, get2
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

  north=(bnd_north)
  northpe=pe_id(pe_mex,pe_mey+1)
  south=(bnd_south)
  southpe=pe_id(pe_mex,pe_mey-1)
  east=(bnd_east)
  eastpe=pe_id(pe_mex+1,pe_mey)
  west=(bnd_west)
  westpe=pe_id(pe_mex-1,pe_mey)

  if     (north .and. west) then
    put1 = Z2
    put2 = Z7
  else if(north .and. east) then
    put1 = Z1
    put2 = Z4
  else if(south .and. west) then
    put1 = Z5
    put2 = Z8
  else if(south .and. east) then
    put1 = Z3
    put2 = Z6
  else if(north) then
    put1 = Z1
    put2 = Z2
  else if(south) then
    put1 = Z6
    put2 = Z5
  else if(west)  then
    put1 = Z7
    put2 = Z8
  else if(east)  then
    put1 = Z3
    put2 = Z4
  else
  endif

  if     (north .and. west) then
  else if(north .and. east) then
  else if(south .and. west) then
  else if(south .and. east) then
  else if(north) then
  else if(south) then
  else if(west)  then
  else if(east)  then
  else
  endif

  if     (north .and. west) then
    B = get1
    D = get2
  else if(north .and. east) then
    A = get1
    C = get2
  else if(south .and. west) then
    A = get1
    C = get2
  else if(south .and. east) then
    B = get1
    D = get2
  else if(north) then
    A = get1
    B = get2
  else if(south) then
    C = get1
    D = get2
  else if(west)  then
    A = get1
    D = get2
  else if(east)  then
    B = get1
    C = get2
  else
  endif


end subroutine
