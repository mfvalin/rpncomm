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

  type, bind(C)    :: rpncomm_array
    type(c_ptr)    :: p     ! pointer to data
    type(c_int)    :: t1    ! tag1 (reserved for rpn_comm usage)
    type(c_int)    :: t2    ! tag2 (reserved for rpn_comm usage)
  end type rpncomm_array
  type(rpncomm_array), parameter :: NULL_rpncomm_array = rpncomm_array(C_NULL_PTR,0,0)
