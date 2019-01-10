module rpn_comm_ezwin_mod
  use ISO_C_BINDING
  implicit none
  include 'RPN_COMM.inc'
  type :: rpn_comm_ezwin
    type(C_PTR) :: p
    integer     :: ix
    integer     :: win
    integer     :: sz
  end type
  integer, parameter :: RPN_COMM_MAGIC2=Z'CAFEFADE'
  integer, parameter :: MAX_EZWINDOWS=64
  type(rpn_comm_ezwin), dimension(MAX_EZWINDOWS), save :: ezwtab
  integer, save :: ent=0
  contains
  function is_invalid(mywin) result(status)
    implicit none
    type(rpncomm_window), intent(IN) :: mywin
    logical :: status

    status = .true.
    if(.not. C_ASSOCIATED(mywin%p)) return
    if( mywin%t1 .ne. ieor(mywin%t2,RPN_COMM_MAGIC2)) return
    status = .false.
  end function is_invalid
end module
subroutine RPN_COMM_ezwin_create(sz,comm,mywin,ierr)
  use rpn_comm_ezwin_mod
  implicit none
  include 'mpif.h'
  integer, intent(IN) :: sz,comm
  integer, intent(OUT) :: ierr
  type(rpncomm_window), intent(OUT) :: mywin
  integer(KIND=MPI_ADDRESS_KIND) :: sz8

  mywin%p  = C_NULL_PTR    ! precondition to error condition
  mywin%t1 = 0
  mywin%t2 = 0

  if(ent >= MAX_EZWINDOWS) return    ! table is full
  ent = ent + 1    ! bump table pointer

  sz8 = sz
  call MPI_WIN_ALLOCATE(sz8, 4, MPI_INFO_NULL, comm, ezwtab(ent)%p,ierr)
  if(ierr .ne. MPI_SUCCESS) return    !! error creating window

  mywin%p = ezwtab(ent)%p
  mywin%t1 = ieor(ent,RPN_COMM_MAGIC2)
  mywin%t2 = ent
end subroutine RPN_COMM_ezwin_create

function RPN_COMM_ezwin_id(mywin) result(winid)
  use rpn_comm_ezwin_mod
  implicit none
  include 'mpif.h'
  type(rpncomm_window), intent(IN) :: mywin
  integer :: winid

  winid = MPI_WIN_NULL
  if(is_invalid(mywin)) return
  winid = ezwtab(mywin%t2)%win
end function RPN_COMM_ezwin_id

function RPN_COMM_ezwin_ptr(mywin) result(winptr)
  use rpn_comm_ezwin_mod
  implicit none
  include 'mpif.h'
  type(rpncomm_window), intent(IN) :: mywin
  type(C_PTR) :: winptr

  winptr = ezwtab(mywin%t2)%p
end function RPN_COMM_ezwin_ptr

function RPN_COMM_ezwin_size(mywin) result(sz)
  use rpn_comm_ezwin_mod
  implicit none
  include 'mpif.h'
  type(rpncomm_window), intent(IN) :: mywin
  integer :: sz

  sz = -1
  if(is_invalid(mywin)) return
  sz = ezwtab(mywin%t2)%sz
end function RPN_COMM_ezwin_size

subroutine RPN_COMM_ezwin_fetch_add(mywin, add, fetch, rank, offset, ierr)
  use rpn_comm_ezwin_mod
  implicit none
  include 'mpif.h'
  type(rpncomm_window), intent(IN) :: mywin
  integer, intent(IN) :: add, rank, offset
  integer, intent(OUT) :: fetch, ierr
  integer(KIND=MPI_ADDRESS_KIND) target_disp
  integer :: assert

  ierr = ieor(MPI_SUCCESS,RPN_COMM_MAGIC2)   ! fudge an error code
  if(is_invalid(mywin)) return
  target_disp = offset
  assert = 0
  call MPI_Win_lock(MPI_LOCK_EXCLUSIVE, rank, assert, ezwtab(mywin%t2)%win, ierr)
  call MPI_FETCH_AND_OP(add, fetch, MPI_INTEGER, rank, target_disp, MPI_SUM, ezwtab(mywin%t2)%win, ierr)
  call MPI_Win_unlock(rank, ezwtab(mywin%t2)%win, ierr)
end subroutine RPN_COMM_ezwin_fetch_add

subroutine RPN_COMM_ezwin_get(mywin,z,nw,rank,offset,ierr)
  use rpn_comm_ezwin_mod
  implicit none
  include 'mpif.h'
  type(rpncomm_window), intent(IN) :: mywin
  integer, dimension(*), intent(OUT) :: z
  integer, intent(IN) :: nw, rank, offset
  integer, intent(OUT) :: ierr
  integer(KIND=MPI_ADDRESS_KIND) target_disp
  integer :: assert

  ierr = ieor(MPI_SUCCESS,RPN_COMM_MAGIC2)   ! fudge an error code
  if(is_invalid(mywin)) return
  target_disp = offset
  assert = 0
  call MPI_Win_lock(MPI_LOCK_SHARED, rank, assert, ezwtab(mywin%t2)%win, ierr)
  call MPI_get(z, nw, MPI_INTEGER, rank, target_disp, nw, MPI_INTEGER, ezwtab(mywin%t2)%win, ierr)
  call MPI_Win_unlock(rank, ezwtab(mywin%t2)%win, ierr)
end subroutine RPN_COMM_ezwin_get

subroutine RPN_COMM_ezwin_put(mywin,z,nw,rank,offset,ierr)
  use rpn_comm_ezwin_mod
  implicit none
  include 'mpif.h'
  type(rpncomm_window), intent(IN) :: mywin
  integer, dimension(*), intent(IN) :: z
  integer, intent(IN) :: nw, rank, offset
  integer, intent(OUT) :: ierr
  integer(KIND=MPI_ADDRESS_KIND) target_disp
  integer :: assert

  ierr = ieor(MPI_SUCCESS,RPN_COMM_MAGIC2)   ! fudge an error code
  if(is_invalid(mywin)) return
  target_disp = offset
  assert = 0
  call MPI_Win_lock(MPI_LOCK_SHARED, rank, assert, ezwtab(mywin%t2)%win, ierr)
  call MPI_get(z, nw, MPI_INTEGER, rank, target_disp, nw, MPI_INTEGER, ezwtab(mywin%t2)%win, ierr)
  call MPI_Win_unlock(rank, ezwtab(mywin%t2)%win, ierr)
end subroutine RPN_COMM_ezwin_put
