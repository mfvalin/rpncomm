!
! typed blind C pointers used for user application interface
!
  integer, parameter :: RPNCOMM_types_version = 1

  type, bind(C)    :: rpncomm_context
    type(c_ptr)    :: p
  end type rpncomm_context
  type(rpncomm_context), parameter :: NULL_rpncomm_context = rpncomm_context(C_NULL_PTR)

  type, bind(C)    :: rpncomm_communicator
    type(c_ptr)    :: p
  end type rpncomm_communicator
  type(rpncomm_communicator), parameter :: NULL_rpncomm_communicator = rpncomm_communicator(C_NULL_PTR)

  type, bind(C)    :: rpncomm_datatype
    type(c_ptr)    :: p
  end type rpncomm_datatype
  type(rpncomm_datatype), parameter :: NULL_rpncomm_datatype = rpncomm_datatype(C_NULL_PTR)

  type, bind(C)    :: rpncomm_operator
    type(c_ptr)    :: p
  end type rpncomm_operator
  type(rpncomm_operator), parameter :: NULL_rpncomm_operator = rpncomm_operator(C_NULL_PTR)

  type, bind(C)    :: rpncomm_group
    type(c_ptr)    :: p
  end type rpncomm_group
  type(rpncomm_group), parameter :: NULL_rpncomm_group = rpncomm_group(C_NULL_PTR)

  type, bind(C)    :: rpncomm_ptr
    type(c_ptr)    :: p
  end type rpncomm_ptr
  type(rpncomm_ptr), parameter :: NULL_rpncomm_ptr = rpncomm_ptr(C_NULL_PTR)

  type, bind(C)    :: rpncomm_bnds
    integer(c_int)    :: lo, hi, np     ! low index, high index, nb of useful points
  end type rpncomm_bnds
  type(rpncomm_bnds), parameter :: NULL_rpncomm_bnds = rpncomm_bnds(0,0,0)

  type, bind(C)    :: rpncomm_field     ! 1/2/3D array, represented as as 3D array
    type(c_ptr)    :: p                 ! pointer to data
    type(c_ptr)    :: m                 ! pointer to metadata
    type(rpncomm_bnds) :: x, y, z       ! dimension bounds (up to 3 dimensions)
    integer(c_int) :: hx, hy, hz        ! halos after exchange
  end type rpncomm_field
  type(rpncomm_field), parameter :: NULL_rpncomm_field = rpncomm_field(         &
                       C_NULL_PTR, C_NULL_PTR,                                  &
                       NULL_rpncomm_bnds, NULL_rpncomm_bnds, NULL_rpncomm_bnds, &
                       0, 0, 0 )

  type, bind(C)    :: rpncomm_vector    ! 1D vector
    type(c_ptr)    :: p                 ! pointer to data
    type(c_ptr)    :: m                 ! pointer to metadata
    integer(c_int) :: x                 ! dimension
  end type rpncomm_vector
  type(rpncomm_vector), parameter :: NULL_rpncomm_vector = rpncomm_vector(C_NULL_PTR, C_NULL_PTR, 0)

  type, bind(C)    :: rpncomm_array     ! 2D array
    type(c_ptr)    :: p                 ! pointer to data
    type(c_ptr)    :: m                 ! pointer to metadata
    integer(c_int) :: x, y              ! dimensions
  end type rpncomm_array
  type(rpncomm_array), parameter :: NULL_rpncomm_array = rpncomm_array(C_NULL_PTR, C_NULL_PTR, 0, 0)

  type, bind(C)    :: rpncomm_block     ! 3D array
    type(c_ptr)    :: p                 ! pointer to data
    type(c_ptr)    :: m                 ! pointer to metadata
    integer(c_int) :: x, y, z           ! dimensions
  end type rpncomm_block
  type(rpncomm_block), parameter :: NULL_rpncomm_block = rpncomm_block(C_NULL_PTR, C_NULL_PTR, 0, 0, 0)
