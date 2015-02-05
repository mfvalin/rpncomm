!
module RPN_COMM_io_pe_tables
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
  function create_ioset(npes,method) result(setno)
    use rpn_comm
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
      if(pe_mex == io_set(iosets)%x(i) .and. pe_mey == io_set(iosets)%y(i)) then  ! I am an IO pe ?
        my_color = 1
        exit
      endif
    enddo
    call MPI_COMM_SPLIT(pe_medomm,my_color,pe_me,io_set(iosets)%comm,ierr)   ! communicator for this set
    
  end function create_ioset
end module RPN_COMM_io_pe_tables

!
subroutine RPN_COMM_io_pes
  use rpn_comm
  use RPN_COMM_io_pe_tables
  implicit none
  return
end subroutine RPN_COMM_io_pes
