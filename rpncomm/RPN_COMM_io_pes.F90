!
module RPN_COMM_io_pe_tables
  use rpn_comm
  type :: RPN_COMM_io_set
    integer, dimension(:), pointer :: x
    integer, dimension(:), pointer :: y
    integer :: ioset
    integer :: comm
    integer :: me
    integer :: npe
  end type
  integer, save :: iosets=0
  integer, PARAMETER :: MAXIOSETS=8
  type(RPN_COMM_io_set), dimension(MAXIOSETS), save :: io_set
contains
!
  function is_io_pe(setno) result(ordinal)  ! is this pe part of IO PE set setno. if yes return rank in set, else return -1
    implicit none
    integer, intent(IN) :: setno
    integer :: ordinal
    integer :: i

    ordinal = -1
    do i = 1 , io_set(setno)%npe
      if(pe_mex == io_set(setno)%x(i) .and. pe_mey == io_set(setno)%y(i)) then  ! I am an IO pe
        ordinal =  io_set(setno)%me
        return
      endif
    enddo
    
  end function is_io_pe
!
  function create_ioset(npes,method) result(setno)
    implicit none
    integer, intent(IN) :: npes, method
    integer :: setno
    integer :: ierr, my_color, i
!
    setno = -1
    if(iosets >= MAXIOSETS) return   ! OOPS, table is full
    iosets = iosets + 1
!
    allocate(io_set(iosets)%x(npes))
    io_set(iosets)%x = -1
    allocate(io_set(iosets)%y(npes))
    io_set(iosets)%y = -1
    io_set(iosets)%ioset = iosets
    io_set(iosets)%npe = npes
    io_set(iosets)%comm = MPI_COMM_NULL
    io_set(iosets)%me = -1
    my_color = MPI_UNDEFINED
    do i = 1 , npes
      if(pe_mex == io_set(iosets)%x(i) .and. pe_mey == io_set(iosets)%y(i)) then  ! I am an IO pe
        my_color = 1
        exit
      endif
    enddo
    call MPI_COMM_SPLIT(pe_medomm,my_color,pe_me,io_set(iosets)%comm,ierr)   ! communicator for this set
    if(io_set(iosets)%comm .ne. MPI_UNDEFINED) then                          ! get rank in communicator if part of set
      call mpi_comm_rank(io_set(iosets)%comm,io_set(iosets)%me,ierr)
    endif
    
  end function create_ioset
!
  subroutine make_io_pe_list(listx,listy,npe,method)
  implicit none
  integer, intent(IN) :: npe, method
  integer, intent(OUT), dimension(npe) :: listx, listy

  listx = -1
  listy = -1
  end subroutine make_io_pe_list
!
end module RPN_COMM_io_pe_tables
!
subroutine RPN_COMM_set_io_pes(setno,npe,method)
  use RPN_COMM_io_pe_tables
  implicit none
  integer, intent(OUT) :: setno
  integer, intent(IN) :: npe, method

  setno = create_ioset(npe,method)
  return
end subroutine RPN_COMM_set_io_pes
!
function RPN_COMM_is_io_pe(setno) result(ordinal)
  use RPN_COMM_io_pe_tables
  implicit none
  integer, intent(IN) :: setno
  integer :: ordinal
  ordinal = is_io_pe(setno)
end function RPN_COMM_is_io_pe
