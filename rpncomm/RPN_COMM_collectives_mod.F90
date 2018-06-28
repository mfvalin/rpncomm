module rpn_comm_collectives_mod   ! used for intra node fast short collectives
use ISO_C_BINDING 
implicit none
integer, parameter :: COLLECTIVE_BUFFER_SIZE = 32          ! COLLECTIVE_BUFFER_SIZE in KBytes

integer, save :: collective_shared_tag = -1                ! shared memory tag for this region
type(C_PTR), save :: collective_shared_mem = C_NULL_PTR    ! shared memory pointer on SMP node for small fast in node collectives

type, bind(C) :: collective_node_buffer
  integer(C_INT) :: count
  integer(C_INT), dimension(COLLECTIVE_BUFFER_SIZE*256) :: data   ! translate KBytes into integers
end type

type(collective_node_buffer), save, pointer :: collective_buffer => NULL()
end module rpn_comm_collectives_mod
