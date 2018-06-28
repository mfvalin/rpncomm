module rpn_comm_server_mod                     ! setup data used by rpn_comm_init and server code
  use ISO_C_BINDING 
  use rpn_comm_fiol_buf_mod
  implicit none
  integer, save :: internal_shared_tag = -1                   ! shared memory tag for this region
  type(C_PTR), save :: internal_shared_mem = C_NULL_PTR       ! shared memory pointer on SMP node for internal use (WORLD)

  integer, save :: server_shared_tag = -1                     ! shared memory tag for this region
  type(C_PTR), save :: server_shared_mem = C_NULL_PTR         ! shared memory pointer on SMP node for server PEs (GRID)

  type(C_FUNPTR), save :: fnserv = C_NULL_FUNPTR ! user supplied function to call if this process is a server
  PROCEDURE ( ), POINTER :: pserv => NULL()      ! Fortran equivalent of above
  integer, dimension(:), pointer, save :: imem => NULL()      ! shared memory area (Fortran version of server_shared_mem)

  integer, save :: nservers = 0                  ! number of server PEs in group
  integer, save :: ncompute = 0                  ! number of compute PEs in group (npegroup-nservers)
  integer, save :: npegroup = 0                  ! number of PEs in PE group (including the servers)
  integer, save :: mbytes = 1                    ! number of MegaBytes of shared memory for servers

  integer, parameter :: COMMAND_BUFFER_SIZE = 2048
  integer, parameter :: MAX_COMMANDS        = 128
  integer, parameter :: CMD_NONE            = 0
  integer, parameter :: CMD_QUEUED          = 1
  integer, parameter :: CMD_ACTIVE          = 2
  integer, parameter :: CMD_DONE            = 4

  type, bind(C) :: cmd_channel
    integer(C_INT) :: id                          ! server ID ( 0 <= ID <= nservers-1 )
    integer(C_INT) :: tag                         ! command tag (used for status)
    type(fiol_buf) :: command                     ! commands to server PE
    type(fiol_buf) :: reply                       ! replies from server PE
    integer(C_INT8_T), dimension(0:MAX_COMMANDS-1)    :: status ! command status (none/in/processing/done)
  end type

  ! in == out                  : buffer is empty
  ! in + 1 modulo limit == out : buffer is full
!   type :: command_channel
!     integer :: id                                 ! server ID ( 0 <= ID <= nservers-1 )
!     integer :: tag                                ! command tag (used for status)
!     integer :: first, limit                       ! first, limit : index limits for buf 
!     integer :: in, out                            ! in : insertion point : out : extraction point (commands to server)
!     integer :: inr, outr                          ! inr : insertion point : outr : extraction point (reply from server)
!     integer, dimension(0:COMMAND_BUFFER_SIZE-1) :: buf     ! buffer for commaond to server PE
!     integer, dimension(0:COMMAND_BUFFER_SIZE-1) :: reply   ! buffer for replies from server PE
!     integer *1 , dimension(0:MAX_COMMANDS-1)    :: status ! command status (none/in/processing/done)
!   end type

!   type(command_channel), dimension(:,:), pointer, save :: cmem => NULL()  ! command channels between all Compute and Server processes
  type(cmd_channel), dimension(:,:), pointer, save :: cmem => NULL()  ! command channels between all Compute and Server processes
                                                                          ! dimension (nservers , ncompute)
!   type(command_channel), dimension(:), pointer, save :: cmem_me => NULL() ! command channels between this compute PE and Server PEs
  type(cmd_channel), dimension(:), pointer, save :: cmem_me => NULL() ! command channels between this compute PE and Server PEs

  contains
!=============== USE storage_size to find size in bits of a derived type ==========
  function init_server_channels(server, pe_me) result(status)
    implicit none
    integer, intent(IN), value :: server
    integer, intent(IN), value :: pe_me
    integer :: status
    integer :: j

    do j = 1, npegroup-nservers
      cmem(server,j)%id    = pe_me / nservers   ! server group ( 0 <= server group < number of groups )
      cmem(server,j)%tag   = 0
      cmem(server,j)%command%first = 1
      cmem(server,j)%command%in    = 1
      cmem(server,j)%command%out   = 1
      cmem(server,j)%command%limit = COMMAND_BUFFER_SIZE
!       cmem(server,j)%command%data  = C_NULL_PTR
      cmem(server,j)%reply%first   = 1
      cmem(server,j)%reply%in      = 1
      cmem(server,j)%reply%out     = 1
      cmem(server,j)%reply%limit   = COMMAND_BUFFER_SIZE
!       cmem(server,j)%reply%data    = C_NULL_PTR
      cmem(server,j)%status = 0
    enddo
    status = 0
  end function init_server_channels
end module rpn_comm_server_mod
