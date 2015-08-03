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
    type(C_PTR), intent(IN) :: win_ptr               ! pointer to entry in win_tab
    logical :: is_valid
    type(rpncomm_windef), pointer :: win_entry

    is_valid = .false.
    if(.not. c_associated(win_ptr)) return           ! pointer to table entry is not valid_win
    if(.not. associated(win_tab)) return                    ! win_tab does not exist yet

    call c_f_ptr(win_ptr,win_entry)
    if(.not. c_associated(win_entry%base)) return    ! no array associated with entry

  end function valid_win_entry

  subroutine create_win_entry(base,type,size,comm,ierr)
    implicit none
    type(C_PTR), intent(IN) :: base        ! base address of array exposed through window
    integer, intent(IN) :: type            ! MPI data type
    integer, intent(IN) :: size            ! number of elements
    integer, intent(IN) :: comm            ! communicator for one sided window
    integer, intent(OUT) :: ierr           ! return status RPN_COMM_ERROR or RPN_COMM_OK
    integer :: i, extent, ierror

    if(.not. associated(win_tab)) call init_win_tab
    ierr = RPN_COMM_ERROR

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

subroutine RPN_COMM_is_win_create(window,type,size,communicator,array,ierr)
  use RPN_COMM_windows
  implicit none
!  include 'RPN_COMM_types.inc'
!  include 'RPN_COMM_constants.inc'
  integer, intent(OUT) :: ierr
  type(rpncomm_window), intent(OUT) :: window
  type(rpncomm_datatype), intent(IN) :: type
  integer, intent(IN) :: size
  type(rpncomm_communicator), intent(IN) :: communicator
  type(C_PTR), intent(IN) :: array

  ierr = RPN_COMM_ERROR
  window = NULL_rpncomm_window
end subroutine RPN_COMM_is_win_create

subroutine RPN_COMM_is_win_free(window,ierr)
  use RPN_COMM_windows
  implicit none
!  include 'RPN_COMM_types.inc'
!  include 'RPN_COMM_constants.inc'
  integer, intent(OUT) :: ierr
  type(rpncomm_window), intent(INOUT) :: window

  ierr = RPN_COMM_ERROR
end subroutine RPN_COMM_is_win_free

subroutine RPN_COMM_is_win_open(window,ierr)
  use RPN_COMM_windows
  implicit none
!  include 'RPN_COMM_types.inc'
!  include 'RPN_COMM_constants.inc'
  integer, intent(OUT) :: ierr
  type(rpncomm_window), intent(IN) :: window

  ierr = RPN_COMM_ERROR
end subroutine RPN_COMM_is_win_open

subroutine RPN_COMM_is_win_close(window,ierr)
  use RPN_COMM_windows
  implicit none
!  include 'RPN_COMM_types.inc'
!  include 'RPN_COMM_constants.inc'
  integer, intent(OUT) :: ierr
  type(rpncomm_window), intent(IN) :: window

  ierr = RPN_COMM_ERROR
end subroutine RPN_COMM_is_win_close

function RPN_COMM_is_win_check(window,ierr) result(is_open)
  use RPN_COMM_windows
  implicit none
!  include 'RPN_COMM_types.inc'
!  include 'RPN_COMM_constants.inc'
  integer, intent(OUT) :: ierr
  type(rpncomm_window), intent(IN) :: window
  logical :: is_open

  ierr = RPN_COMM_ERROR
  is_open = .false.
end function RPN_COMM_is_win_check

subroutine RPN_COMM_is_win_put_r(window,larray,target,offset,nelem,ierr)
  use RPN_COMM_windows
  implicit none
!  include 'RPN_COMM_types.inc'
!  include 'RPN_COMM_constants.inc'
  integer, intent(OUT) :: ierr
  type(rpncomm_window), intent(IN) :: window
  type(C_PTR), intent(IN) :: larray
  integer, intent(IN) :: target
  integer, intent(IN) :: offset
  integer, intent(IN) :: nelem

  ierr = RPN_COMM_ERROR
end subroutine RPN_COMM_is_win_put_r

subroutine RPN_COMM_is_win_put_l(window,larray,offset,nelem,ierr)
  use RPN_COMM_windows
  implicit none
!  include 'RPN_COMM_types.inc'
!  include 'RPN_COMM_constants.inc'
  integer, intent(OUT) :: ierr
  type(rpncomm_window), intent(IN) :: window
  type(C_PTR), intent(IN) :: larray
  integer, intent(IN) :: offset
  integer, intent(IN) :: nelem

  ierr = RPN_COMM_ERROR
end subroutine RPN_COMM_is_win_put_l

subroutine RPN_COMM_is_win_get_r(window,larray,target,offset,nelem,ierr)
  use RPN_COMM_windows
  implicit none
!  include 'RPN_COMM_types.inc'
!  include 'RPN_COMM_constants.inc'
  integer, intent(OUT) :: ierr
  type(rpncomm_window), intent(IN) :: window
  type(C_PTR), intent(IN) :: larray
  integer, intent(IN) :: target
  integer, intent(IN) :: offset
  integer, intent(IN) :: nelem

  ierr = RPN_COMM_ERROR
end subroutine RPN_COMM_is_win_get_r

subroutine RPN_COMM_is_win_get_l(window,larray,offset,nelem,ierr)
  use RPN_COMM_windows
  implicit none
!  include 'RPN_COMM_types.inc'
!  include 'RPN_COMM_constants.inc'
  integer, intent(OUT) :: ierr
  type(rpncomm_window), intent(IN) :: window
  type(C_PTR), intent(IN) :: larray
  integer, intent(IN) :: offset
  integer, intent(IN) :: nelem

  ierr = RPN_COMM_ERROR
end subroutine RPN_COMM_is_win_get_l
