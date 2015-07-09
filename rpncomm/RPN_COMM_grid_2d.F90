module rpn_comm_grids
  type :: rpn_comm_dist_grid
    integer :: grid_id
    integer :: gni, gnj
    integer :: mini, maxi, lni
    integer :: minj, maxj, lnj
    integer, dimension(:), pointer :: start_i, count_i
    integer, dimension(:), pointer :: start_j, count_j
  end type

  integer, parameter :: MAX_GRIDS = 64
  type(rpn_comm_dist_grid), dimension(MAX_GRIDS), save :: gt

  integer :: used_grids = 0

contains
  function init_grid_slot(ix, nx, ny) result(indx)   ! initialize a (maybe new) slot in grid table
    implicit none
    integer, intent(IN) :: ix, nx, ny
    integer :: indx

    if(ix > MAX_GRIDS .or. ix <= 0) then
      indx = -1
      return   ! error, bad slot number
    endif

    gt(ix)%grid_id = -1
    gt(ix)%gni     = 0
    gt(ix)%gni     = 0
    gt(ix)%lni     = 0
    gt(ix)%mini    = 0
    gt(ix)%maxi    = 0
    gt(ix)%lnj     = 0
    gt(ix)%minj    = 0
    gt(ix)%maxj    = 0

    if( .not. associated(gt(ix)%start_i) ) then
      allocate(gt(ix)%start_i(nx), gt(ix)%count_i(nx))
      allocate(gt(ix)%start_j(ny), gt(ix)%count_j(ny))
    endif

    gt(ix)%start_i = 0
    gt(ix)%count_i = 0

    gt(ix)%start_j = 0
    gt(ix)%count_j = 0

    indx = ix

    return
  end function init_grid_slot

  function init_new_grid(nx, ny) result(indx)   ! initialize a (maybe new) slot in grid table
    implicit none
    integer, intent(IN) :: nx, ny
    integer :: indx

    integer :: i, ix

    if(used_grids == 0) then
      do i = 1 , MAX_GRIDS
        nullify(gt(i)%start_i)
        nullify(gt(i)%count_i)
        nullify(gt(i)%start_j)
        nullify(gt(i)%count_j)
        gt(i)%grid_id = -1
      enddo
    endif

    indx = -1
    do i = 1 , MAX_GRIDS
      if(gt(i)%grid_id <= 0) then  ! unused grid
        ix = i
        exit
      endif
    enddo
    if(ix == -1) return  ! no free slot found, table is full

    indx = init_grid_slot(ix, nx, ny)  ! fill slot with zeroes

    return
  end function init_new_grid

  function find_grid(id) result(indx)   ! find index in table gt associated to grid_id id
    implicit none
    integer, intent(IN) :: id
    integer :: indx

    integer :: i

    indx = -1
    i = and(id,Z'000000FF')
    if(i <= 0 .or. i > MAX_GRIDS) return  ! invalid index, not a grid_id
    if(gt(i)%grid_id /= id)       return  ! wrong grid id, index is not coherent
    indx = i

  end function find_grid

end module rpn_comm_grids

function rpn_comm_create_2dgrid(gni,gnj,mini,maxi,minj,maxj) result (grid_id)  !InTf!
  use rpn_comm
  use rpn_comm_grids
  implicit none
  integer, intent(IN) :: gni, gnj, mini, maxi, minj, maxj              !InTf!
  integer :: grid_id                                                   !InTf!

  integer :: ix, i, j, lni, lnj

  grid_id = RPN_COMM_ERROR
  ix = init_new_grid(pe_nx,pe_ny)
  if(ix <= 0) return                ! id <= 0 is an error

  lni = (gni + pe_nx -1) / pe_nx    ! max lni
  lnj = (gnj + pe_ny -1) / pe_ny    ! max lnj
  if( 1-mini > maxi-lni) return     ! halo x problem (halo x is assumed to be 1-mini)
  if( 1-minj > maxj-lnj) return     ! halo y problem (halo y is assumed to be 1-minj)

  gt(ix)%gni     = gni
  gt(ix)%gnj     = gnj
  gt(ix)%mini    = mini
  gt(ix)%maxi    = maxi
  gt(ix)%minj    = minj
  gt(ix)%maxj    = maxj

  gt(ix)%count_i = lni
  gt(ix)%count_i(pe_nx) = gni - lni * (pe_nx - 1)
  gt(ix)%start_i(1) = 0
  do i = 2 , pe_nx
    gt(ix)%start_i(i) = gt(ix)%start_i(i-1) + gt(ix)%count_i(i-1)
  enddo
  gt(ix)%lni     = gt(ix)%count_i(pe_mex + 1)    ! adjust local lni

  gt(ix)%count_j = lnj
  gt(ix)%count_j(pe_ny) = gnj - lnj * (pe_ny - 1)
  gt(ix)%start_j(1) = 0
  do j = 2 , pe_ny
    gt(ix)%start_j(j) = gt(ix)%start_j(j-1) + gt(ix)%count_j(j-1)
  enddo
  gt(ix)%lnj     = gt(ix)%count_j(pe_mey + 1)    ! adjust local lnj

  gt(ix)%grid_id = or(ix , Z'CAFE0000')
  grid_id = gt(ix)%grid_id

end function rpn_comm_create_2dgrid                                    !InTf!

function rpn_comm_get_2dgrid(grid_id,gni,gnj,mini,maxi,minj,maxj,starti,counti,startj,countj) result (status)  !InTf!
  use rpn_comm
  use rpn_comm_grids
  implicit none
  integer, intent(IN) :: grid_id                                                   !InTf!
  integer, intent(OUT) :: gni,gnj,mini,maxi,minj,maxj                              !InTf!
  integer, intent(OUT), dimension(:) :: starti,counti,startj,countj                !InTf!
  integer :: status                                                                !InTf!

  integer :: indx

  status = RPN_COMM_ERROR
  indx = find_grid(grid_id)
  if(indx <= 0 .or. indx > MAX_GRIDS) return

  if(size(starti) .ne. size(gt(indx)%start_i)) return
  if(size(counti) .ne. size(gt(indx)%count_i)) return
  if(size(startj) .ne. size(gt(indx)%start_j)) return
  if(size(countj) .ne. size(gt(indx)%count_j)) return

  gni    = gt(indx)%gni
  gnj    = gt(indx)%gnj
  mini   = gt(indx)%mini
  maxi   = gt(indx)%maxi
  minj   = gt(indx)%minj
  maxj   = gt(indx)%maxj
  starti = gt(indx)%start_i
  counti = gt(indx)%count_i
  startj = gt(indx)%start_j
  countj = gt(indx)%count_j

  status = RPN_COMM_OK
end function rpn_comm_get_2dgrid                                                   !InTf!
