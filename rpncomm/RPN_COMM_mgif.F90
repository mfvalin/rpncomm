module mpi_mgi_internals
  use ISO_C_BINDING
  implicit none
  save
  integer, parameter :: MAX_NAMES = 8
  character(C_CHAR), dimension(32,MAX_NAMES) :: names
  character(C_CHAR), dimension(MAX_NAMES) :: modes
  integer(C_INT) :: nch = 0
  interface
    !int MPI_mgi_create(const char *alias, char mode, int cpl);
    function mgi_create(name, mode, cpl) result(status) bind(C,name='MPI_mgi_create')
       import :: C_CHAR, C_INT
       character(C_CHAR), dimension(*), intent(IN) :: name
       character(C_CHAR), intent(IN), value :: mode
       integer(C_INT), intent(IN), value :: cpl
       integer(C_INT) :: status
    end function mgi_create
    !int MPI_mgi_create_begin(int cpl, MPI_Comm f_comm);
    function mgi_create_begin(cpl, comm) result(status) bind(C,name='MPI_mgi_create_begin_f')
       import :: C_CHAR, C_INT
       integer(C_INT), intent(IN), value :: cpl
       integer(C_INT), intent(IN), value :: comm
       integer(C_INT) :: status
    end function mgi_create_begin
    !int MPI_mgi_create_end(int cpl, MPI_Comm f_comm);
    function mgi_create_end(cpl, comm) result(status) bind(C,name='MPI_mgi_create_end_f')
       import :: C_CHAR, C_INT
       integer(C_INT), intent(IN), value :: cpl
       integer(C_INT), intent(IN), value :: comm
       integer(C_INT) :: status
    end function mgi_create_end
    ! int MPI_mgi_init(const char *alias);
    function MPI_mgi_init(alias) result(status) bind(C,name='MPI_mgi_init')
      import :: C_CHAR, C_INT
      character(C_CHAR), dimension(*), intent(IN) :: alias
      integer(C_INT) :: status
    end function MPI_mgi_init
    ! int MPI_mgi_open(int mpi_channel, unsigned char *mode);
    function MPI_mgi_open(mpi_channel, mode) result(status) bind(C,name='MPI_mgi_open')
      import :: C_CHAR, C_INT
      integer(C_INT), intent(IN), value :: mpi_channel
      character(C_CHAR), intent(IN), value :: mode
      integer(C_INT) :: status
    end function MPI_mgi_open
    ! int MPI_mgi_read(int channel, unsigned char *data, int nelm, unsigned char *dtyp);
    function MPI_mgi_read(channel, data, nelm, dtyp) result(status) bind(C,name='MPI_mgi_read')
      import :: C_CHAR, C_INT
      integer(C_INT), intent(IN), value :: channel, nelm
      integer(C_INT), dimension(*), intent(OUT) :: data
      character(C_CHAR), intent(IN), value :: dtyp
      integer(C_INT) :: status
    end function MPI_mgi_read
    ! int MPI_mgi_write(int channel, unsigned char *data, int nelm, unsigned char *dtyp);
    function MPI_mgi_write(channel, data, nelm, dtyp) result(status) bind(C,name='MPI_mgi_write')
      import :: C_CHAR, C_INT
      integer(C_INT), intent(IN), value :: channel, nelm
      integer(C_INT), dimension(*), intent(IN) :: data
      character(C_CHAR), intent(IN), value :: dtyp
      integer(C_INT) :: status
    end function MPI_mgi_write
    ! int MPI_mgi_close(int channel);
    function MPI_mgi_close(channel) result(status) bind(C,name='MPI_mgi_close')
      import :: C_INT
      integer(C_INT), intent(IN), value :: channel
      integer(C_INT) :: status
    end function MPI_mgi_close
  end interface
end module

subroutine register_mgi_channel_name(name, mode)
  use mpi_mgi_internals
  implicit none
  character(len=*), intent(IN) :: name
  character(len=1), intent(IN) :: mode

  nch = nch + 1
  names(:,nch) = transfer( trim(name)//achar(0) , names(:,nch) )
  modes(nch) = mode
  return
end subroutine register_mgi_channel_name

subroutine register_mgi_channels(cpl,comm)
  use mpi_mgi_internals
  implicit none
  integer, intent(IN) :: cpl,comm
  integer :: i, status

  status = mgi_create_begin(cpl,comm)
  do i=1,nch
    status = mgi_create(names(:,i),modes(i),cpl)
  enddo
  status = mgi_create_end(cpl,comm)

  return
end subroutine register_mgi_channels

function mgi_init(alias) result(status)
  use mpi_mgi_internals
  implicit none
  character(len=*), intent(IN) :: alias
  integer :: status
  character(C_CHAR), dimension(MAX_NAMES) :: temp

  temp = transfer(trim(alias)//achar(0),temp)
  status = MPI_mgi_init(temp)
end function mgi_init

function mgi_open(mpi_channel, mode) result(status)
  use mpi_mgi_internals
  implicit none
  integer, intent(IN) :: mpi_channel
  character(len=1), intent(IN) :: mode
  integer :: status

  status = MPI_mgi_open(mpi_channel, mode)
  return
end function mgi_open

function mgi_read(channel, data, nelm, dtyp)  result(status)
  use mpi_mgi_internals
  implicit none
  integer, intent(IN) :: channel, nelm
  integer, dimension(*), intent(OUT) :: data
  character(len=1), intent(IN) :: dtyp
  integer :: status

  status = MPI_mgi_read(channel, data, nelm, dtyp)
  return
end function mgi_read

function mgi_write(channel, data, nelm, dtyp)  result(status)
  use mpi_mgi_internals
  implicit none
  integer, intent(IN) :: channel, nelm
  integer, dimension(*), intent(IN) :: data
  character(len=1), intent(IN) :: dtyp
  integer :: status

  status = MPI_mgi_write(channel, data, nelm, dtyp)
  return
end function mgi_write

function mgi_close(channel) result(status)
  use mpi_mgi_internals
  implicit none
  integer, intent(IN) :: channel
  integer :: status

  status = MPI_mgi_close(channel)
  return
end function mgi_close

function mgi_term() result(status)
  use ISO_C_BINDING
  implicit none
  integer(C_INT) :: status

  status = 0
  return
end function mgi_term

