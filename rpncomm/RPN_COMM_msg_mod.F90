module rpn_comm_msg_mod                                     ! used for in node messaging
  use ISO_C_BINDING
  use rpn_comm_fiol_buf_mod
  implicit none
  integer, parameter :: MSG_BUFFER_SIZE     = 2048            ! size in kilobytes

  integer, save     :: msg_shared_tag = -1                    ! shared memory tag for this region
  type(C_PTR), save :: msg_shared_mem = C_NULL_PTR            ! shared memory pointer on SMP node for messaging (WORLD)

  integer, save :: rank_on_node = -1             ! rank of this process on SMP node in node_comm (origin 0)
  integer, save :: pes_on_node = -1              ! number of processes on SMP node
  integer, save :: node_rank_zero = -1           ! communicator used by processes of rank 0 on SMP nodes, regardless of grid
  integer, save :: node_comm = -1                ! communicator for ALL PEs on SMP node, regardless of grid
  integer, save :: node_grid_comm = -1           ! communicator for PEs on SMP node, belonging to same grid

  type, bind(C):: msg_channel
    type(fiol_buf) :: inbound
    type(fiol_buf) :: outbound
  end type

  type(msg_channel), save, dimension(:), pointer :: msg_channels =>NULL()

  contains
  ! initialize msg_channels table, create table with  pes_on_node entries
  ! initialize first, in, out, limit, and data pointer for each entry
  ! data pointer will point to integer array capable of containing nslots 4 byte values
  ! status = -1 : no memory allocated for table or invalid rank and size on node
  ! status = -2 : table already initialized
  ! status = -3 : not enough memory available for table
  ! status : number of 4Byte items used for table and data buffers
  ! the data buffers are allocated at addresses above the table itself
  function init_msg_channels() result(status) bind(C,name='InitMsgChannels')
    implicit none
    integer(C_INT) :: status
    integer(C_INT) :: nslots    ! number of 4 byte slots in data buffers
    integer :: i

    status = -1
    if(.not. c_associated(msg_shared_mem)) return  ! no shared memory was allocated , OUCH !!
    if(pes_on_node < 0)return                      ! bad value, OUCH !!
    if(rank_on_node < 0)return                     ! bad value, OUCH !!
    status = -2
    if(associated(msg_channels))           return  ! initialization already done , OOPS !!

    status = 2 * FIOL_BUF_SIZE * pes_on_node        ! total size of msg_channels table in 4Byte items
    nslots = ( (MSG_BUFFER_SIZE * pes_on_node * 1024) - (status * 4) ) / pes_on_node   ! available number of slots per data buffer 
    if(nslots < 128) then
      status = -3
      return
    endif

    call C_F_Pointer(msg_shared_mem,msg_channels,[pes_on_node])

    if(rank_on_node > 0)   return    ! job done if not rank 0, barrier at return point highly recommended

    do i = 1, pes_on_node
      msg_channels(i)%inbound%first = status    ! configure inbound buffer for entry i
      msg_channels(i)%inbound%in    = msg_channels(i)%inbound%first
      msg_channels(i)%inbound%out   = msg_channels(i)%inbound%first
      msg_channels(i)%inbound%limit = msg_channels(i)%inbound%first + nslots
      status = status + nslots              ! increment address by nslots * 4 bytes (nslots integers)

      msg_channels(i)%outbound%first = status   ! configure outbound buffer for entry i
      msg_channels(i)%outbound%in    = msg_channels(i)%outbound%first
      msg_channels(i)%outbound%out   = msg_channels(i)%outbound%first
      msg_channels(i)%outbound%limit = msg_channels(i)%outbound%first + nslots
      status = status + nslots              ! increment address by nslots * 4 bytes (nslots integers)
    enddo
  end function init_msg_channels

  function msg_from_pe(pe, data, ndata) result(status) ! get message data from a PE  (used by MASTER)
    implicit none
    type(C_PTR), intent(IN), value    :: data      ! pointer to data
    integer(C_INT), intent(IN), value :: pe        ! target pe (ordinal in node communicator) (origin 0)
    integer(C_INT), intent(IN), value :: ndata     ! number of 4 byte items to read from buffer
    integer(C_INT) :: status

    status = fiol_extract(msg_shared_mem, msg_channels(pe+1)%outbound, data, ndata, 1)   ! +1  because indexing msg_channels from 1

  end function msg_from_pe

  function msg_from_me(data, ndata)  result(status) ! get message data from MY inbound buffer (used by WORKER)
    implicit none
    type(C_PTR), intent(IN), value    :: data      ! pointer to data
    integer(C_INT), intent(IN), value :: ndata     ! number of 4 byte items to read from buffer
    integer(C_INT) :: status

    status = fiol_extract(msg_shared_mem, msg_channels(rank_on_node+1)%outbound, data, ndata, 1)   ! +1  because indexing msg_channels from 1

  end function msg_from_me

  function msg_to_pe(pe, data, ndata) result(status)   ! send message data to a PE  (used by MASTER)
    implicit none
    type(C_PTR), intent(IN), value    :: data      ! pointer to data
    integer(C_INT), intent(IN), value :: pe        ! target pe (ordinal in node communicator) (origin 0)
    integer(C_INT), intent(IN), value :: ndata     ! number of 4 byte items to write into buffer
    integer(C_INT) :: status

    status = fiol_insert(msg_shared_mem, msg_channels(pe+1)%inbound, data, ndata, 1)   ! +1  because indexing msg_channels from 1

  end function msg_to_pe

  function msg_to_me(data, ndata) result(status)    ! post message data into MY outbound buffer (used by WORKER)
    implicit none
    type(C_PTR), intent(IN), value    :: data      ! pointer to data
    integer(C_INT), intent(IN), value :: ndata     ! number of 4 byte items to write into buffer
    integer(C_INT) :: status

    status = fiol_insert(msg_shared_mem, msg_channels(rank_on_node+1)%inbound, data, ndata, 1)   ! +1  because indexing msg_channels from 1

  end function msg_to_me

end module rpn_comm_msg_mod
