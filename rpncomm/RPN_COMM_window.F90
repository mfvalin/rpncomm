!/! RPN_COMM - Library of useful routines for C and FORTRAN programming
! ! Copyright (C) 1975-2015  Division de Recherche en Prevision Numerique
! !                          Environnement Canada
! !
! ! This library is free software; you can redistribute it and/or
! ! modify it under the terms of the GNU Lesser General Public
! ! License as published by the Free Software Foundation,
! ! version 2.1 of the License.
! !
! ! This library is distributed in the hope that it will be useful,
! ! but WITHOUT ANY WARRANTY; without even the implied warranty of
! ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! ! Lesser General Public License for more details.
! !
! ! You should have received a copy of the GNU Lesser General Public
! ! License along with this library; if not, write to the
! ! Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! ! Boston, MA 02111-1307, USA.
! !/
! one sided communication window management code
!
module RPN_COMM_windows
  use ISO_C_BINDING
  implicit none
  include 'mpif.h'
  include 'RPN_COMM_types.inc'
  include 'RPN_COMM_types_int.inc'
  include 'RPN_COMM_constants.inc'
  integer, parameter :: RPN_COMM_MAX_WINDOWS = 64
  type(rpncomm_windef), dimension(:), pointer, save :: win_tab => NULL()
contains
  subroutine init_win_tab
    implicit none
    integer :: i

    if(ASSOCIATED(win_tab)) return                  ! already initialized
    allocate(win_tab(RPN_COMM_MAX_WINDOWS))
    do i = 1 , RPN_COMM_MAX_WINDOWS
      win_tab(i) = NULL_rpncomm_windef
    enddo
  end subroutine init_win_tab

  function valid_win_entry(win_ptr) result(is_valid)
    implicit none
    type(C_PTR), intent(IN), value :: win_ptr               ! pointer to entry in win_tab
    logical :: is_valid
    type(rpncomm_windef), pointer :: win_entry

    is_valid = .false.
    if(.not. c_associated(win_ptr)) return           ! pointer to table entry is not valid_win
    if(.not. associated(win_tab)) return                    ! win_tab does not exist yet

    call c_f_ptr(win_ptr,win_entry)
    if(.not. c_associated(win_entry%base)) return    ! no array associated with entry

    is_valid = .true.

  end function valid_win_entry

  subroutine create_win_entry(base,type,size,comm,indx,ierr)
    implicit none
    type(C_PTR), intent(IN), value :: base        ! base address of array exposed through window
    integer, intent(IN) :: type            ! MPI data type
    integer, intent(IN) :: size            ! number of elements
    integer, intent(IN) :: comm            ! communicator for one sided window
    integer, intent(OUT) :: indx           ! index into wintab of created window
    integer, intent(OUT) :: ierr           ! return status RPN_COMM_ERROR or RPN_COMM_OK
    integer :: i, extent, ierror

    if(.not. associated(win_tab)) call init_win_tab
    ierr = RPN_COMM_ERROR                  ! preset for failure
    indx = -1

    do i = 1 , RPN_COMM_MAX_WINDOWS
      if(.not. c_associated(win_tab(i)%base) ) then
        if(c_associated(base)) then    ! user has already allocated space
          win_tab(i)%base = base
        else                           ! allocate needed space
          call MPI_TYPE_EXTENT(type, extent, ierror)
          if(ierror .ne. MPI_SUCCESS) return     ! invalid type ? 
          win_tab(i)%base = base
        endif
        win_tab(i)%is_open = .false.
        win_tab(i)%is_user = c_associated(base)
        win_tab(i)%type = type
        win_tab(i)%size = size
        ierr = RPN_COMM_OK
        return
      endif
    enddo
  end subroutine create_win_entry
end module RPN_COMM_windows

!InTf!
!===============================================================================
! create a one sided communication window
!
! window (OUT)     rpn_comm window type returned to caller (see RPN_COMM_types.inc)
! dtype  (IN)      rpn_comm datatype descriptor (see RPN_COMM_types.inc)
! siz    (IN)      number of elements of type dtype in window
! com    (IN)      rpn_comm communicator used for window (see RPN_COMM_types.inc)
! array  (IN)      C pointer to array associated with window
!                  if defined (not equal to C_NULL_PTR), this user array will be used
!                  if not defined (equal to C_NULL_PTR), an internal array will be allocated and used
! ierr   (OUT)     RPN_COMM_OK or RPN_COMM_ERROR will be returned
!===============================================================================
subroutine RPN_COMM_i_win_create(window,dtype,siz,com,array,ierr)  !InTf!
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window, rpncomm_datatype, rpncomm_communicator  !InTf!
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(OUT) :: window                         !InTf!
  type(rpncomm_datatype), intent(IN) :: dtype                         !InTf!
  integer, intent(IN) :: siz                                          !InTf!
  type(rpncomm_communicator), intent(IN) :: com                       !InTf!
  type(C_PTR), intent(IN), value :: array
! !  integer :: array                                                  !InTf!
! !  !DEC$ ATTRIBUTES NO_ARG_CHECK :: array                            !InTf!
! !  !GCC$ ATTRIBUTES NO_ARG_CHECK :: array                            !InTf!
! !  !IBM* ignore_tkr array                                            !InTf!
! !  !DIR$ ignore_tkr array                                            !InTf!
! !  !$PRAGMA ignore_tkr array                                         !InTf!
  integer :: indx

  ierr = RPN_COMM_ERROR
  window = NULL_rpncomm_window

  call create_win_entry(array,dtype,siz,com,indx,ierr)
  if(ierr .ne. RPN_COMM_OK) return

  window%p = c_loc(win_tab[indx])
  window%t1 = indx
  window%t2 = -indx

  ierr = RPN_COMM_OK
  return

end subroutine RPN_COMM_i_win_create                                  !InTf!

!InTf!
subroutine RPN_COMM_i_win_free(window,ierr)                           !InTf!
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(INOUT) :: window                       !InTf!

  ierr = RPN_COMM_ERROR
end subroutine RPN_COMM_i_win_free                                    !InTf!

!InTf!
subroutine RPN_COMM_i_win_open(window,ierr)                           !InTf!
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!

  ierr = RPN_COMM_ERROR
end subroutine RPN_COMM_i_win_open                                    !InTf!

!InTf!
subroutine RPN_COMM_i_win_close(window,ierr)                          !InTf!
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!

  ierr = RPN_COMM_ERROR
end subroutine RPN_COMM_i_win_close                                   !InTf!

!InTf!
function RPN_COMM_i_win_valid(window,ierr) result(is_valid)           !InTf!
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  logical :: is_valid                                                 !InTf!
  type(C_PTR) :: temp

  ierr = RPN_COMM_ERROR
  is_valid = .false.

  if(.not. c_associated(window%p)) return    ! no win_tab entry pointer
  if(window%t1 + window%t2 .ne. 0) return    ! bad tags
  if(window%t1 < 0 .or. window%t1>RPN_COMM_MAX_WINDOWS) return  ! invalid index

  temp = c_loc(win_tab(window%t1))
  if( .not.c_associated(window%p,temp)) return     ! window%p must point to proper entry in window table

  is_valid = valid_win_entry(window%p)             ! check that window description is good

end function RPN_COMM_i_win_valid                                     !InTf!

!InTf!
function RPN_COMM_i_win_check(window,ierr) result(is_open)            !InTf!
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  logical :: is_open                                                  !InTf!

  ierr = RPN_COMM_ERROR
  is_open = .false.
end function RPN_COMM_i_win_check                                     !InTf!

!InTf!
subroutine RPN_COMM_i_win_put_r(window,larray,target,offset,nelem,ierr) !InTf!
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  type(C_PTR), intent(IN), value :: larray                            !InTf!
  integer, intent(IN) :: target                                       !InTf!
  integer, intent(IN) :: offset                                       !InTf!
  integer, intent(IN) :: nelem                                        !InTf!

  ierr = RPN_COMM_ERROR
end subroutine RPN_COMM_i_win_put_r                                   !InTf!

!InTf!
subroutine RPN_COMM_i_win_put_l(window,larray,offset,nelem,ierr)      !InTf!
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  type(C_PTR), intent(IN), value :: larray                            !InTf!
  integer, intent(IN) :: offset                                       !InTf!
  integer, intent(IN) :: nelem                                        !InTf!

  ierr = RPN_COMM_ERROR
end subroutine RPN_COMM_i_win_put_l                                   !InTf!

!InTf!
subroutine RPN_COMM_i_win_get_r(window,larray,target,offset,nelem,ierr) !InTf!
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  type(C_PTR), intent(IN), value :: larray                            !InTf!
  integer, intent(IN) :: target                                       !InTf!
  integer, intent(IN) :: offset                                       !InTf!
  integer, intent(IN) :: nelem                                        !InTf!

  ierr = RPN_COMM_ERROR
end subroutine RPN_COMM_i_win_get_r                                   !InTf!

!InTf!
subroutine RPN_COMM_i_win_get_l(window,larray,offset,nelem,ierr)      !InTf!
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  type(C_PTR), intent(IN), value :: larray                            !InTf!
  integer, intent(IN) :: offset                                       !InTf!
  integer, intent(IN) :: nelem                                        !InTf!

  ierr = RPN_COMM_ERROR
end subroutine RPN_COMM_i_win_get_l                                   !InTf!
