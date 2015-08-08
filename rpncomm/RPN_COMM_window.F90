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
!===============================================================================
! one sided communication window management code
!===============================================================================
module RPN_COMM_windows
  use ISO_C_BINDING
  implicit none
  include 'mpif.h'
  include 'RPN_COMM_types.inc'
  include 'RPN_COMM_types_int.inc'
  include 'RPN_COMM_constants.inc'
  integer, parameter :: RPN_COMM_MAX_WINDOWS = 64
  integer, save :: integer_size = 0
  type(rpncomm_windef), dimension(:), pointer, save :: win_tab => NULL()  ! rpn comm window table
contains
!===============================================================================
! allocate and initialize internal single sided communication wwindows table
!===============================================================================
  subroutine init_win_tab
    implicit none
    integer :: i

    if(ASSOCIATED(win_tab)) return                  ! already initialized
    allocate(win_tab(RPN_COMM_MAX_WINDOWS))
    do i = 1 , RPN_COMM_MAX_WINDOWS
      win_tab(i) = NULL_rpncomm_windef              ! null entry
      win_tab(i)%com = MPI_COMM_NULL                ! with MPI null communicator
      win_tab(i)%typ = MPI_DATATYPE_NULL            ! and MPI null datatype
      win_tab(i)%win = MPI_WIN_NULL                 ! MPI null window
    enddo
    call MPI_TYPE_EXTENT(MPI_INTEGER, integer_size, i) ! get size of MPI_INTEGER
  end subroutine init_win_tab

!===============================================================================
! check if win_ptr (pointing into win_tab) points to a valid entry
!===============================================================================
  function valid_win_entry(win_ptr) result(is_valid)
    implicit none
    type(C_PTR), intent(IN), value :: win_ptr        ! pointer to entry in win_tab
    logical :: is_valid

    type(rpncomm_windef), pointer :: win_entry
    type(C_PTR) :: temp

    is_valid = .false.
    if(.not. c_associated(win_ptr)) return           ! pointer to table entry is not valid
    if(.not. associated(win_tab)) return             ! win_tab does not exist yet

    call c_f_ptr(win_ptr,win_entry)                  ! make pointer to win_entry from win_ptr

    if(win_entry%indx <=0 .or. win_entry%indx>RPN_COMM_MAX_WINDOWS) return  ! impossible index into window table
    temp = c_loc(win_tab(win_entry%indx))            ! address pointed to by entry index
    if(.not. c_associated(temp,win_ptr)) return      ! entry index should point to itself
    if(.not. c_associated(win_entry%base)) return    ! no array associated with entry
    if(win_entry%com == MPI_COMM_NULL) return        ! no communicator
    if(win_entry%typ == MPI_DATATYPE_NULL) return    ! no datatype
    if(win_entry%win == MPI_WIN_NULL) return         ! no window
    if(win_entry%siz <= 0) return                    ! invalid size

    is_valid = .true.                                ! no more reason to think entry is not valid

  end function valid_win_entry

!===============================================================================
! create a new one sided communication window
!
! if argument base is a NULL pointer (C_NULL_PTR), allocate an internal array
! using the recommended MPI_alloc_mem MPI routine
!===============================================================================
  subroutine create_win_entry(base,typ,siz,comm,indx,ierr)
    implicit none
    type(C_PTR), intent(IN), value :: base ! base address of array exposed through window (may be C_NULL_PTR)
    integer, intent(IN) :: typ             ! MPI data type
    integer, intent(IN) :: siz             ! number of elements in array associated to this window
    integer, intent(IN) :: comm            ! MPI communicator for this one sided window
    integer, intent(OUT) :: indx           ! index into wintab of created window (used for consistency test)
    integer, intent(OUT) :: ierr           ! return status (RPN_COMM_ERROR or RPN_COMM_OK)
    integer :: i, extent, ierror

    if(.not. associated(win_tab)) call init_win_tab  ! create and initialize window table if necessary
    ierr = RPN_COMM_ERROR                  ! preset for failure
    indx = -1

    do i = 1 , RPN_COMM_MAX_WINDOWS        ! loop over window table entries (we are looking for a free entry)
      if(.not. c_associated(win_tab(i)%base) ) then  ! this entry is free
        call MPI_TYPE_EXTENT(typ, extent, ierror)  ! determine size associated with MPI datatype
        if( mod(extent,integer_size) .ne. 0 ) return    ! extent of data type not a multiple of integer
        if(ierror .ne. MPI_SUCCESS) return         ! invalid type ? other error ?
        if(c_associated(base)) then    ! user has already allocated space
          win_tab(i)%base = base
        else                           ! allocate needed space with MPI_alloc_mem
          call MPI_alloc_mem(extent*siz, MPI_INFO_NULL, win_tab(i)%base,ierr)
          if(ierr .ne. MPI_SUCCESS) return       ! could not allocate memory
        endif
        indx = i
        call c_f_ptr(win_tab(i)%base,win_tab(i)%remote,extent*siz)
        win_tab(i)%is_open = .false.             ! window is not "exposed"
        win_tab(i)%is_user = c_associated(base)  ! user supplied space ?
        win_tab(i)%typ = typ                     ! MPI data type associated to window
        win_tab(i)%ext = extent/integer_size     ! size (extent) of MPI data type in MPI_INTEGER units
        win_tab(i)%indx = indx                   ! entry index points to itself
        win_tab(i)%siz = siz                     ! number of elements in window
        win_tab(i)%com = comm                    ! communicator associated with window
        ierr = RPN_COMM_OK                       ! all is well
        return
      endif
    enddo
    return  ! if we fall through the loop, the table was full, ew will return RPN_COMM_ERROR
  end subroutine create_win_entry
!===============================================================================
! check if window is indeed a valid one
!===============================================================================
  function win_valid(window) result(is_valid)
  implicit none
  type(rpncomm_window), intent(IN) :: window
  logical :: is_valid
  type(C_PTR) :: temp

  is_valid = .false.
  if(.not. associated(win_tab)) call init_win_tab  ! create and initialize window table if necessary

  if( ieor(window%t1,RPN_COMM_MAGIC) .ne. window%t2 )   return  ! inconsistent tags
  if(window%t2 < 0 .or. window%t2>RPN_COMM_MAX_WINDOWS) return  ! invalid index

  temp = c_loc(win_tab(window%t2))                 ! address of relevant win_tab entry
  if( .not.c_associated(window%p,temp)) return     ! window%p must point to win_tab(window%t2)

  is_valid = valid_win_entry(window%p)             ! check that window description in table is good
  return

end function win_valid

end module RPN_COMM_windows

!InTf!
!===============================================================================
! create a one sided communication window (user exposed interface)
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

  call create_win_entry(array,dtype%t2,siz,com%t2,indx,ierr)
  if(ierr .ne. RPN_COMM_OK) return

  window%p = c_loc(win_tab(indx))        ! point to entry in window table
  window%t1 = ieor(indx,RPN_COMM_MAGIC)  ! xor with magic token
  window%t2 = indx                       ! index into table

  ierr = RPN_COMM_OK
  return

end subroutine RPN_COMM_i_win_create                                  !InTf!

!===============================================================================
! delete a previously created one sided communication window (see RPN_COMM_i_win_create)
!===============================================================================
!InTf!
subroutine RPN_COMM_i_win_free(window,ierr)                           !InTf!
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(INOUT) :: window                       !InTf!

  integer :: indx

  ierr = RPN_COMM_ERROR
  if(.not. win_valid(window) ) return

  indx = window%t2
  window = NULL_rpncomm_window

  if(.not. win_tab(indx)%is_user) then   ! internal storage allocation, release it
    call MPI_free_mem(win_tab(indx)%base,ierr)
  endif
  win_tab(indx) = NULL_rpncomm_windef              ! blank entry
  win_tab(indx)%com = MPI_COMM_NULL                ! MPI null communicator
  win_tab(indx)%typ = MPI_DATATYPE_NULL            ! MPI null datatype
  win_tab(indx)%win = MPI_WIN_NULL                 ! MPI null window

  ierr = RPN_COMM_OK
  return

end subroutine RPN_COMM_i_win_free                                    !InTf!

!===============================================================================
! "expose" a one sided communication window (see RPN_COMM_i_win_create)
!===============================================================================
!InTf!
subroutine RPN_COMM_i_win_open(window,ierr)                           !InTf!
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!

  integer :: ierr2, indx
  logical ::is_open
  logical, external :: RPN_COMM_i_win_check

  ierr = RPN_COMM_ERROR
  indx = window%t2                   ! window table entry for this window
  is_open = RPN_COMM_i_win_check(window,ierr2)
  if(is_open .or. ierr2 .ne. RPN_COMM_OK) return    ! ERROR: window is already open (exposed) or not valid

  call MPI_win_fence(0,win_tab(indx)%win,ierr2)

  if(ierr2 .eq. MPI_SUCCESS) ierr = RPN_COMM_OK
  return

end subroutine RPN_COMM_i_win_open                                    !InTf!

!===============================================================================
! stop "exposing" a one sided communication window (see RPN_COMM_i_win_create)
! the result of all remotely performed get/put operations may now be used
!===============================================================================
!InTf!
subroutine RPN_COMM_i_win_close(window,ierr)                          !InTf!
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!

  integer :: ierr2, indx
  logical :: is_not_open
  logical, external :: RPN_COMM_i_win_check

  ierr = RPN_COMM_ERROR
  indx = window%t2                   ! window table entry for this window
  is_not_open = .not. RPN_COMM_i_win_check(window,ierr2)
  if(is_not_open .or. ierr2 .ne. RPN_COMM_OK) return    ! ERROR: window is not open (exposed) or not valid

  call MPI_win_fence(0,win_tab(indx)%win,ierr2)

  if(ierr2 .eq. MPI_SUCCESS) ierr = RPN_COMM_OK
  return

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
  is_valid = win_valid(window)
  if(.not. is_valid ) return             ! check that window description is valid
  ierr = RPN_COMM_OK
  return

end function RPN_COMM_i_win_valid                                     !InTf!

!===============================================================================
! check if a one sided communication window (see RPN_COMM_i_win_create) is "exposed"
!===============================================================================
!InTf!
function RPN_COMM_i_win_check(window,ierr) result(is_open)            !InTf!
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  logical :: is_open                                                  !InTf!

  integer :: indx

  ierr = RPN_COMM_ERROR
  is_open = .false.
  if(.not. win_valid(window) ) return

  indx = window%t2                   ! window table entry for this window
  is_open = win_tab(indx)%is_open    ! get open (exposed) flag from window table entry

  ierr = RPN_COMM_OK
  return

end function RPN_COMM_i_win_check                                     !InTf!

!===============================================================================
! one sided communication remote put (write) into one sided communication window
! from a local array
!===============================================================================
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

  logical :: is_open
  integer :: ierr2, indx
  logical, external :: RPN_COMM_i_win_check
  integer, dimension(:), pointer :: local
  type(rpncomm_windef), pointer :: win_entry

  ierr = RPN_COMM_ERROR
  is_open = RPN_COMM_i_win_check(window,ierr2)
  if( (.not. is_open) .or. (ierr2 .ne. MPI_SUCCESS) )  return  ! bad window reference or window not open (exposed)

  indx = window%t2
  call c_f_ptr( window%p , win_entry )                   ! pointer to win_tab entry (rpncomm_windef)
  if(offset+nelem > win_entry%siz) return                ! out of bounds condition for "remote" array
  call c_f_ptr( larray , local, nelem* win_entry%ext )   ! pointer to local array

  call MPI_put(local,nelem,win_entry%typ,target,offset,nelem,win_entry%typ,win_entry%win,ierr2)

  if(ierr2 == MPI_SUCCESS) ierr = RPN_COMM_OK
  return

end subroutine RPN_COMM_i_win_put_r                                   !InTf!

!===============================================================================
! one sided communication local put (write) into one sided communication window
! from a local array
!===============================================================================
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

  logical :: is_open
  integer :: ierr2, i, indx, extent
  type(rpncomm_windef), pointer :: win_entry
  integer, dimension(:), pointer :: local
  logical, external :: RPN_COMM_i_win_check

  ierr = RPN_COMM_ERROR
  is_open = RPN_COMM_i_win_check(window,ierr2)
  if( (is_open) .or. (ierr2 .ne. MPI_SUCCESS) )  return  ! bad window reference or window open (exposed)

  indx = window%t2
  call c_f_ptr( window%p , win_entry )                   ! pointer to win_tab entry (rpncomm_windef)
  if(offset+nelem > win_entry%siz) return                ! out of bounds condition for "remote" array
  extent = win_entry%ext
  call c_f_ptr( larray , local, nelem*extent )                   ! pointer to local array
  do i = 1, nelem*extent
    win_entry%remote(i+offset*extent) = local(i)
  enddo

  ierr = RPN_COMM_OK
  return

end subroutine RPN_COMM_i_win_put_l                                   !InTf!

!===============================================================================
! one sided communication remote get (read) from one sided communication window
! into a local array
!===============================================================================
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

  logical :: is_open
  integer :: ierr2, indx
  logical, external :: RPN_COMM_i_win_check
  integer, dimension(:), pointer :: local
  type(rpncomm_windef), pointer :: win_entry

  ierr = RPN_COMM_ERROR
  is_open = RPN_COMM_i_win_check(window,ierr2)
  if( (.not. is_open) .or. (ierr2 .ne. MPI_SUCCESS) )  return  ! bad window reference or window not open (exposed)

  indx = window%t2
  call c_f_ptr( window%p , win_entry )                   ! pointer to win_tab entry (rpncomm_windef)
  if(offset+nelem > win_entry%siz) return                ! out of bounds condition for "remote" array
  call c_f_ptr( larray , local, nelem* win_entry%ext )   ! pointer to local array

  call MPI_get(local,nelem,win_entry%typ,target,offset,nelem,win_entry%typ,win_entry%win,ierr2)

  if(ierr2 == MPI_SUCCESS) ierr = RPN_COMM_OK
  return

end subroutine RPN_COMM_i_win_get_r                                   !InTf!

!===============================================================================
! one sided communication local get (read) from one sided communication window
! into a local array
!===============================================================================
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

  logical :: is_open
  integer :: ierr2, i, indx, extent
  type(rpncomm_windef), pointer :: win_entry
  integer, dimension(:), pointer :: local
  logical, external :: RPN_COMM_i_win_check

  ierr = RPN_COMM_ERROR
  is_open = RPN_COMM_i_win_check(window,ierr2)
  if( (.not. is_open) .or. (ierr2 .ne. MPI_SUCCESS) )  return  ! bad window reference or window open (exposed)

  indx = window%t2
  call c_f_ptr( window%p , win_entry )                   ! pointer to win_tab entry (rpncomm_windef)
  if(offset+nelem > win_entry%siz) return                ! out of bounds condition for "remote" array
  extent = win_entry%ext
  call c_f_ptr( larray , local, nelem*extent )                   ! pointer to local array
  do i = 1, nelem*extent
     local(i) = win_entry%remote(i+offset*extent)
  enddo

  ierr = RPN_COMM_OK
  return

end subroutine RPN_COMM_i_win_get_l                                   !InTf!
