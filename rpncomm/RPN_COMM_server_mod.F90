module rpn_comm_server_mod                     ! setup data used by rpn_comm_init and server code
use ISO_C_BINDING 
implicit none
type(C_PTR), save :: internal_shared_mem = C_NULL_PTR       ! shared memery on SMP node for internal use
type(C_PTR), save :: server_shared_mem = C_NULL_PTR         ! shared memery on SMP node for server PEs
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
! in == out                  : buffer is empty
! in + 1 modulo limit == out : buffer is full
type :: command_channel
  integer :: id                                 ! server ID ( 0 <= ID <= nservers-1 )
  integer :: tag                                ! command tag (used for status)
  integer :: first, limit                       ! first, limit : index limits for buf 
  integer :: in, out                            ! in : insertion point : out : extraction point (commands to server)
  integer :: inr, outr                          ! inr : insertion point : outr : extraction point (reply from server)
  integer, dimension(0:COMMAND_BUFFER_SIZE-1) :: buf     ! buffer for commaond to server PE
  integer, dimension(0:COMMAND_BUFFER_SIZE-1) :: reply   ! buffer for replies from server PE
  integer *1 , dimension(0:MAX_COMMANDS-1)    :: status ! command status (none/in/processing/done)
end type
type(command_channel), dimension(:,:), pointer, save :: cmem => NULL()  ! command channels between all Compute and Server processes
                                                                        ! dimension (nservers , ncompute)
type(command_channel), dimension(:), pointer, save :: cmem_me => NULL() ! command channels between this compute PE and Server PEs
end module rpn_comm_server_mod
