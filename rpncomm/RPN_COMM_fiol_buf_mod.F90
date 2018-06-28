module rpn_comm_fiol_buf_mod
  use ISO_C_BINDING
  implicit none

  type, bind(C) :: fiol_buf   ! circular buffer, usable without locks if there is only one writer and one reader
    integer(C_INT) :: first   ! and the roles are never swapped (only one thread/process accessing in and out)
    integer(C_INT) :: in
    integer(C_INT) :: out
    integer(C_INT) :: limit
  end type
  integer, parameter :: FIOL_BUF_SIZE = 4  ! 4 integers

  contains

  function fiol_insert(membase, buf, data, nw, blocking) result(status) bind(C,name='FiolBufInsert')  ! insert data into fiol buffer
    implicit none
    type(fiol_buf), intent(INOUT)     :: buf            ! fiol buffer
    type(C_PTR), intent(IN), value    :: membase        ! base for circular buffer indexing
    type(C_PTR), intent(IN), value    :: data           ! user data
    integer(C_INT), intent(IN), value :: nw             ! number of 4 byte tokens
    integer(C_INT), intent(IN), value :: blocking       ! > 0 : blocking/synchronous ; == 0 :  non blocking/asynchronous
    integer(C_INT) :: status                            ! > 0 : number of tokens inserted, < nw if not enough room and non blocking

    integer(C_INT) :: nextin
    integer(C_INT), dimension(:), pointer :: fdata      ! Fortran version of C pointer data
    integer(C_INT), dimension(:), pointer :: fbuf       ! used for circular buffer indexing, Fortran version of C pointer membase

    status = 0
    call c_f_pointer(data,fdata,[nw])
    call c_f_pointer(membase,fbuf,[buf%limit])

    nextin = buf%in + 1
    if(nextin >= buf%limit) nextin = buf%first      ! wrap around
    do while(status < nw)                 ! while user request not satisfied
      do while(nextin == buf%out )      ! wait while buffer full
        if(blocking == 0) return          ! return current count if non blocking call
      enddo
      status = status + 1
      fbuf(nextin) = fdata(status)
      nextin = nextin + 1
      if(nextin >= buf%limit) nextin = buf%first      ! wrap around
      buf%in = nextin
    enddo
  end function fiol_insert

  function fiol_extract(membase, buf, data, nw, blocking) result(status)  bind(C,name='FiolBufExtract')  ! extract data from fiol buffer
    implicit none
    type(fiol_buf), intent(INOUT)     :: buf            ! fiol buffer
    type(C_PTR), intent(IN), value    :: membase        ! base for circular buffer indexing
    type(C_PTR), intent(IN), value    :: data           ! user data
    integer(C_INT), intent(IN), value :: nw             ! number of 4 byte tokens
    integer(C_INT), intent(IN), value :: blocking       ! > 0 : blocking/synchronous ; == 0 :  non blocking/asynchronous
    integer(C_INT) :: status                            ! > 0 : number of tokens extracted, < nw if not enough data and non blocking

    integer(C_INT), dimension(:), pointer :: fdata      ! Fortran version of C pointer data
    integer(C_INT), dimension(:), pointer :: fbuf       ! used for circular buffer indexing, Fortran version of C pointer membase

    status = 0
    call c_f_pointer(data,fdata,[nw])
    call c_f_pointer(membase,fbuf,[buf%limit])

    do while(status < nw)                 ! while user request not satisfied
      do while(buf%in .ne. buf%out)     ! wait while buffer empty
        if(blocking == 0) return          ! return current count if non blocking call
      enddo
      status = status + 1
      fdata(status) = fbuf(buf%out)
      buf%out = buf%out + 1
      if(buf%out >= buf%limit) buf%out = buf%first     ! wrap around
    enddo
  end function fiol_extract

  function fiol_status(buf) result(status)  bind(C,name='FiolBufStatus')  ! is buffer (empty / neither full nor empty / full) ( >0 / <0 / 0 )
    implicit none
    type(fiol_buf), intent(INOUT)     :: buf            ! fiol buffer
    integer(C_INT) :: status
    integer(C_INT) :: nextin   ! next insertion position

    nextin = buf%in + 1
    if(buf%in == buf%out) then
      status = buf%limit - 1   !  buffer is empty, return buffer capacity
    else
      if(buf%in > buf%out) then           ! in - out tokens in buffer
        status = -(buf%limit - 1 - (buf%in - buf%out))  ! number of free slots
      else                                    ! out - in - 1 free slots
        status = -(buf%out - buf%in - 1)  ! number of free slots
      endif
    endif
  end function fiol_status

end module rpn_comm_fiol_buf_mod
