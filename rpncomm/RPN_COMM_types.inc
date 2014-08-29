!
! typed blind C pointers used for user application interface
!
  type, bind(C)    :: rpncomm_context
    type(c_ptr)    :: p
  end type rpncomm_context
  type, bind(C)    :: rpncomm_communicator
    type(c_ptr)    :: p
  end type rpncomm_communicator
  type, bind(C)    :: rpncomm_datatype
    type(c_ptr)    :: p
  end type rpncomm_datatype
  type, bind(C)    :: rpncomm_group
    type(c_ptr)    :: p
  end type rpncomm_group
  type, bind(C)    :: rpncomm_ptr
    type(c_ptr)    :: p
  end type rpncomm_ptr
