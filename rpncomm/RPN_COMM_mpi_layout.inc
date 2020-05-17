! grids and supergrids
  type, BIND(C) :: mpigrid
    integer(C_INT) :: all              ! all PEs, compute and service
    integer(C_INT) :: compute          ! compute PEs
    integer(C_INT) :: service          ! service (Input/Output) PEs
    integer(C_INT) :: same_node        ! PEs in same SMP node
    integer(C_INT) :: same_numa        ! PEs insame numa domain
    integer(C_INT) :: node_peer        ! node to node peer (PEs with same rank in node)
    integer(C_INT) :: numa_peer        ! numa space to numa space peer (PEs with same rank in numa space)
    integer(C_INT) :: grid_peer        ! grid to grid / supergrid to supergrid peer (PEs with same rank in grid/spuergrid)
    integer(C_INT) :: row              ! undefined in a supergrid, information for a grid row
    integer(C_INT) :: column           ! undefined in a supergrid, information for a grid column
  end type
! applications (domains)
  type, BIND(C) :: application
    integer(C_INT) :: all              ! all PEs, compute and service
    integer(C_INT) :: same_node        ! PEs in same SMP node
    integer(C_INT) :: same_numa        ! PEs insame numa domain
  end type
  type, BIND(C) :: grid_hierarchy
    type(application) :: appl          ! application (domain), contains (M identical supergrids)
    type(mpigrid)     :: sgrd          ! supergrid (contains N identical grids)
    type(mpigrid)     :: grid          ! basic computation grid
  end type
! the three usual associates, communicators, rank, and size
  type, BIND(C) :: mpi_layout
    integer(C_INT) :: version = 123456
    type(grid_hierarchy) :: comm       ! communicators
    type(grid_hierarchy) :: rank       ! rank in communicator
    type(grid_hierarchy) :: size       ! size of communicator
  end type

