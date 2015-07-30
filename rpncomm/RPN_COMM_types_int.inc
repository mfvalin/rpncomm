
  type, bind(C)    :: rpncomm_bnds
    integer(c_int) :: lo, hi, lnp, gnp  ! low bound, high bound, local dimension, global dimension
  end type rpncomm_bnds
  type(rpncomm_bnds), parameter :: NULL_rpncomm_bnds = rpncomm_bnds(0,0,0,0)

  type, bind(C)    :: rpncomm_field     ! 1/2/3D array, represented as as 3D array
    type(c_ptr)    :: p                 ! pointer to data
    type(c_ptr)    :: m                 ! pointer to metadata
    integer(c_int) :: sz                ! base data element size in bytes
    type(rpncomm_bnds) :: x, y, z       ! dimension bounds (up to 3 dimensions)
    integer(c_int) :: hx, hy, hz        ! halos for exchange
    integer(c_int) :: ew3d, ns3d        ! user defined data types for halo exchange
  end type rpncomm_field
  type(rpncomm_field), parameter :: NULL_rpncomm_field = rpncomm_field(         &
                       C_NULL_PTR, C_NULL_PTR, 0,                               &
                       NULL_rpncomm_bnds, NULL_rpncomm_bnds, NULL_rpncomm_bnds, &
                       0, 0, 0, 0, 0 )

  type, bind(C)    :: rpncomm_vector    ! 1D vector
    type(c_ptr)    :: p                 ! pointer to data
    type(c_ptr)    :: m                 ! pointer to metadata
    integer(c_int) :: x                 ! dimension
  end type rpncomm_vector
  type(rpncomm_vector), parameter :: NULL_rpncomm_vector = rpncomm_vector(C_NULL_PTR, C_NULL_PTR, 0)

  type, bind(C)    :: rpncomm_matrix     ! 2D array
    type(c_ptr)    :: p                 ! pointer to data
    type(c_ptr)    :: m                 ! pointer to metadata
    integer(c_int) :: x, y              ! dimensions
  end type rpncomm_matrix
  type(rpncomm_matrix), parameter :: NULL_rpncomm_matrix = rpncomm_matrix(C_NULL_PTR, C_NULL_PTR, 0, 0)

  type, bind(C)    :: rpncomm_block     ! 3D array
    type(c_ptr)    :: p                 ! pointer to data
    type(c_ptr)    :: m                 ! pointer to metadata
    integer(c_int) :: x, y, z           ! dimensions
  end type rpncomm_block
  type(rpncomm_block), parameter :: NULL_rpncomm_block = rpncomm_block(C_NULL_PTR, C_NULL_PTR, 0, 0, 0)
!
! Domain type
!
  type domm
     sequence
     character(len=12) nom
     character(len=1024) path
     integer npex, npey
  end type domm
  integer domm_size, domm_num
  type(domm), allocatable, dimension(:) :: pe_domains
!
!  symbol tables
!
  type :: SYMTAB
  integer :: number
  character (len=32) :: string
  end type
