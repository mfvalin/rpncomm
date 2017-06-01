module mgi_channel_names
  use ISO_C_BINDING
  implicit none
  save
  integer, parameter :: MAX_NAMES = 8
  character(C_CHAR), dimension(32,MAX_NAMES) :: names
  character(C_CHAR), dimension(MAX_NAMES) :: modes
  integer(C_INT) :: nch = 0
end module

subroutine register_mgi_channel_name(name, mode)
  use mgi_channel_names
  implicit none
  character(len=*), intent(IN) :: name
  character(len=1), intent(IN) :: mode

  nch = nch + 1
  names(:,nch) = transfer( trim(name)//achar(0) , names(:,nch) )
  modes(nch) = mode
  return
end subroutine register_mgi_channel_name

subroutine register_mgi_channels(cpl,comm)
  use mgi_channel_names
  implicit none
  integer, intent(IN) :: cpl,comm
  interface
    function mgi_create(name, mode, cpl) result(status) bind(C,name='MPI_mgi_create')
       import :: C_CHAR, C_INT
       character(C_CHAR), dimension(*), intent(IN) :: name
       character(C_CHAR), intent(IN), value :: mode
       integer(C_INT), intent(IN), value :: cpl
       integer(C_INT) :: status
    end function mgi_create
    function mgi_create_begin(cpl, comm) result(status) bind(C,name='MPI_mgi_create_begin')
       import :: C_CHAR, C_INT
       integer(C_INT), intent(IN), value :: cpl
       integer(C_INT), intent(IN), value :: comm
       integer(C_INT) :: status
    end function mgi_create_begin
    function mgi_create_end(cpl, comm) result(status) bind(C,name='MPI_mgi_create_end')
       import :: C_CHAR, C_INT
       integer(C_INT), intent(IN), value :: cpl
       integer(C_INT), intent(IN), value :: comm
       integer(C_INT) :: status
    end function mgi_create_end
  end interface
  integer :: i, status

  status = mgi_create_begin(cpl,comm)
  do i=1,nch
    status = mgi_create(names(:,i),modes(i),cpl)
  enddo
  status = mgi_create_end(cpl,comm)

  return
end subroutine register_mgi_channels


