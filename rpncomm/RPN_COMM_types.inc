!
! typed blind C pointers used for user application interface
!
  type, bind(C)    :: rc_context
    type(c_ptr)    :: p
  end type rc_context
  type, bind(C)    :: rc_communicator
    type(c_ptr)    :: p
  end type rc_communicator
  type, bind(C)    :: rc_datatype
    type(c_ptr)    :: p
  end type rc_datatype
  type, bind(C)    :: rc_group
    type(c_ptr)    :: p
  end type rc_group
