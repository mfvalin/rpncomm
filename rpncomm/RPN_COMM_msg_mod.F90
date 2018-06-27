module rpn_comm_msg_mod                                     ! used for in node messaging
use ISO_C_BINDING 
implicit none
integer, parameter :: MSG_BUFFER_SIZE     = 2048

integer, save :: msg_shared_tag = -1                        ! shared memory tag for this region
type(C_PTR), save :: msg_shared_mem = C_NULL_PTR            ! shared memory pointer on SMP node for messaging (WORLD)

integer, save :: rank_on_node = -1             ! rank of this process on SMP node
integer, save :: pes_on_node = -1              ! number of processes on SMP node
integer, save :: node_rank_zero = -1           ! communicator used by processes of rank 0 on SMP nodes
integer, save :: node_comm = -1                ! communicator for ALL PEs on SMP node, regardless of grid
integer, save :: node_grid_comm = -1           ! communicator for PEs on SMP node, belonging to same grid

type :: msg_buffer
  integer :: in
  integer :: out
  integer, dimension(0:MSG_BUFFER_SIZE-1) :: data
end type

type :: msg_channel
  type(msg_buffer) :: input
  type(msg_buffer) :: output
end type

type(msg_channel), save, dimension(:), pointer :: msg_channels =>NULL()
end module rpn_comm_msg_mod
