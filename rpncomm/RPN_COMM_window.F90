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
#if defined(ROBODOC)
!****P* rpn_comm/windows  simplified/restricted version of MPI one sided communications
! DESCRIPTION
!   This is a simplified and restricted version of MPI-2 one sided communications.
!
!   When creating the "communication window", some attributes are 
!   determined for said window and will not change during its life.
!   1 - the MPI communicator
!   2 - the MPI data type of the elements contained in the window
!   3 - the size of the window (in number of elements)
!   4 - an array large enough to contain these elements
!       (this array can either be supplied by the caller or allocated internall.
!
!   all further operations refer to the window by its "identifier" 
!   Fortran type : type(rpncomm_window)  (include 'RPN_COMM_types.inc')
!
!   the window creation is a COLLECTIVE operation, all members of the communicator group
!   must call RPN_COMM_i_win_create
!
!   "exposing" a window and terminating the "exposure" of a window are also a COLLECTIVE operation
!
!   remote get/put operations, i.e. sending read/write requests targeting the window
!   belonging to a remote PE may only happen when the window is "exposed"
!   remote operations are NOT ALLOWED when the window is "not exposed"
!
!   local get/put operations, i.e. reading/writing from/into the window
!   belonging the local PE may only happen when the window is "not exposed"
!   local operations are NOT ALLOWED whe the window is "exposed"
!
!   a typical sequence of operations would be
!   1 - create a one sided communication window (COLLECTIVE operation)
!    repeat as needed
!    {
!      2a - modify the local copy of the window (if and as needed)
!      2b - "expose" the window (COLLECTIVE operation)
!      2c - perform get/put operations targeting remote PEs (as needed on each PE)
!      2d - "end exposing" the window (COLLECTIVE operation)
!      2e - get modifications from the local copy of the window (if and as needed)
!    }
!   3 - free the one sided communication window (COLLECTIVE operation)
!
!   window creation/destruction : RPN_COMM_i_win_create, RPN_COMM_i_win_free
!   window status queries       : RPN_COMM_i_win_valid, RPN_COMM_i_win_check
!   window info queries         : RPN_COMM_i_win_get_ptr
!   window "exposition"         : RPN_COMM_i_win_open, RPN_COMM_i_win_close
!   window get operations       : RPN_COMM_i_win_get_r, RPN_COMM_i_win_get_l
!   window put operations       : RPN_COMM_i_win_put_r, RPN_COMM_i_win_put_l
!
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
!******
#endif
!****iP* rpn_comm/windows  one sided communication window management module
! SYNOPSIS
module RPN_COMM_windows
!===============================================================================
! one sided communication window management code
!===============================================================================
  use ISO_C_BINDING
  implicit none
  include 'mpif.h'
  include 'RPN_COMM_types.inc'
  include 'RPN_COMM_types_int.inc'
  include 'RPN_COMM_constants.inc'
!******
  integer, parameter :: RPN_COMM_MAX_WINDOWS = 64
  integer, save :: integer_size = 0
  type(rpncomm_windef), dimension(:), pointer, save :: win_tab => NULL()  ! rpn comm window table
contains
!****if* RPN_COMM_windows/init_win_tab
! SYNOPSIS
  subroutine init_win_tab
!===============================================================================
! allocate and initialize internal single sided communication windows table
!===============================================================================
!******
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

!****if* RPN_COMM_windows/win_ptr
! SYNOPSIS
  function valid_win_entry(win_ptr) result(is_valid)
!===============================================================================
! check if win_ptr (pointing into win_tab) points to a valid entry
!===============================================================================
!******
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

!****if* RPN_COMM_windows/create_win_entry
! SYNOPSIS
  subroutine create_win_entry(base,typ,siz,comm,indx,ierr)
!===============================================================================
! create a new one sided communication window
!
! if argument base is a NULL pointer (C_NULL_PTR), allocate an internal array
! using the recommended MPI_alloc_mem MPI routine
!===============================================================================
!******
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

!****if* RPN_COMM_windows/win_valid
! SYNOPSIS
  function win_valid(window) result(is_valid)
!===============================================================================
! check if window is indeed a valid one
!===============================================================================
!******
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

!===============================================================================
! beginning of USER CALLABLE routines/functions
!===============================================================================

!InTf!
!****f* rpn_comm/RPN_COMM_i_win_create create a one sided communication window
! SYNOPSIS
subroutine RPN_COMM_i_win_create(window,dtype,siz,com,array,ierr)  !InTf!
!===============================================================================
! create a one sided communication window (user exposed interface)
!===============================================================================
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
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window, rpncomm_datatype, rpncomm_communicator  !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(OUT) :: window                         !InTf!
  type(rpncomm_datatype), intent(IN) :: dtype                         !InTf!
  integer, intent(IN) :: siz                                          !InTf!
  type(rpncomm_communicator), intent(IN) :: com                       !InTf!
  type(C_PTR), intent(IN), value :: array                             !InTf!
! IGNORE
!!! VOID$ array                                                         !InTf!
!******
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
!InTf!
!****f* rpn_comm/RPN_COMM_i_win_free free a one sided communication window
! SYNOPSIS
subroutine RPN_COMM_i_win_free(window,ierr)                           !InTf!
!===============================================================================
! delete a previously created one sided communication window (see RPN_COMM_i_win_create)
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! ierr   (OUT)    error status, RPN_COMM_OK or RPN_COMM_ERROR
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(INOUT) :: window                       !InTf!
!******

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

!InTf!
!****f* rpn_comm/RPN_COMM_i_win_open open a one sided communication window
! SYNOPSIS
subroutine RPN_COMM_i_win_open(window,ierr)                           !InTf!
!===============================================================================
! "expose" a one sided communication window (see RPN_COMM_i_win_create)
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! ierr   (OUT)    error status, RPN_COMM_OK or RPN_COMM_ERROR
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
!******

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

!InTf!
!****f* rpn_comm/RPN_COMM_i_win_close close a one sided communication window
! SYNOPSIS
subroutine RPN_COMM_i_win_close(window,ierr)                          !InTf!
!===============================================================================
! stop "exposing" a one sided communication window (see RPN_COMM_i_win_create)
! the result of all remotely performed get/put operations may now be used
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! ierr   (OUT)    error status, RPN_COMM_OK or RPN_COMM_ERROR
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
!******

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
!****f* rpn_comm/RPN_COMM_i_win_valid check if a one sided communication window is valid
! SYNOPSIS
function RPN_COMM_i_win_valid(window,ierr) result(is_valid)           !InTf!
!===============================================================================
! find if a one sided communication window description is valid
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! ierr   (OUT)    error status, RPN_COMM_OK or RPN_COMM_ERROR
!
! function value : .true. (window description is valid) or .false. (not valid)
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  logical :: is_valid                                                 !InTf!
!******
  type(C_PTR) :: temp

  ierr = RPN_COMM_ERROR
  is_valid = win_valid(window)
  if(.not. is_valid ) return             ! check that window description is valid
  ierr = RPN_COMM_OK
  return

end function RPN_COMM_i_win_valid                                     !InTf!

!InTf!
!****f* rpn_comm/RPN_COMM_i_win_check check if a one sided communication window is "exposed"
! SYNOPSIS
function RPN_COMM_i_win_check(window,ierr) result(is_open)            !InTf!
!===============================================================================
! check if a one sided communication window (see RPN_COMM_i_win_create) is "exposed"
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! ierr   (OUT)    error status, RPN_COMM_OK or RPN_COMM_ERROR
!
! function value : .true. (window "exposed") or .false. (window not "exposed")
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  logical :: is_open                                                  !InTf!
!******

  integer :: indx

  ierr = RPN_COMM_ERROR
  is_open = .false.
  if(.not. win_valid(window) ) return

  indx = window%t2                   ! window table entry for this window
  is_open = win_tab(indx)%is_open    ! get open (exposed) flag from window table entry

  ierr = RPN_COMM_OK
  return

end function RPN_COMM_i_win_check                                     !InTf!

!InTf!
!****f* rpn_comm/RPN_COMM_i_win_get_ptr get data pointer associated to a one sided communication window
! SYNOPSIS
function RPN_COMM_i_win_get_ptr(window,ierr) result(ptr)                 !InTf!
!===============================================================================
! get a one sided communication window (see RPN_COMM_i_win_create) data pointer
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! ierr   (OUT)    error status, RPN_COMM_OK or RPN_COMM_ERROR
!
! function value : C compatible (type(C_PTR)) pointer to the data array associated with window
!                  in case of error, C_NULL_PTR is returned (null pointer)
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  type(C_PTR) :: ptr                                                  !InTf!
!******

  integer :: indx

  ierr = RPN_COMM_ERROR
  ptr = C_NULL_PTR
  if(.not. win_valid(window) ) return

  indx = window%t2            ! window table entry for this window
  ptr = win_tab(indx)%base    ! get data pointer from window table entry

  ierr = RPN_COMM_OK
  return

end function RPN_COMM_i_win_get_ptr                                      !InTf!

!InTf!
!****f* rpn_comm/RPN_COMM_i_win_put_r write into a remote one sided communication window
! SYNOPSIS
subroutine RPN_COMM_i_win_put_r(window,larray,target,offset,nelem,ierr) !InTf!
!===============================================================================
! one sided communication remote put (write) into one sided communication window
! from a local array
! it is an error to attempt a "remote" put when the window is not "exposed"
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! larray (IN)     C compatible pointer (type(C_PTR)) to local array (source of put)
! target (IN)     ordinal in window communicator of remote PE
! offset (IN)     displacement (origin 0) into remote PE window data array
! nelem  (IN)     number of elements to transfer (type of element was defined at window creation)
! ierr   (OUT)    error status, RPN_COMM_OK or RPN_COMM_ERROR
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  type(C_PTR), intent(IN), value :: larray                            !InTf!
  integer, intent(IN) :: target                                       !InTf!
  integer, intent(IN) :: offset                                       !InTf!
  integer, intent(IN) :: nelem                                        !InTf!
!******

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

!InTf!
!****f* rpn_comm/RPN_COMM_i_win_put_l write into a local one sided communication window
! SYNOPSIS
subroutine RPN_COMM_i_win_put_l(window,larray,offset,nelem,ierr)      !InTf!
!===============================================================================
! one sided communication local put (write) into one sided communication window
! from a local array
! it is an error to attempt a "local" put when the window is "exposed"
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! larray (IN)     C compatible pointer (type(C_PTR)) to local array (source of put)
! offset (IN)     displacement (origin 0) into this PE window data array
! nelem  (IN)     number of elements to transfer (type of element was defined at window creation)
! ierr   (OUT)    error status, RPN_COMM_OK or RPN_COMM_ERROR
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  type(C_PTR), intent(IN), value :: larray                            !InTf!
  integer, intent(IN) :: offset                                       !InTf!
  integer, intent(IN) :: nelem                                        !InTf!
!******

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

!InTf!
!****f* rpn_comm/RPN_COMM_i_get_r read a remote one sided communication window
! SYNOPSIS
subroutine RPN_COMM_i_win_get_r(window,larray,target,offset,nelem,ierr) !InTf!
!===============================================================================
! one sided communication remote get (read) from one sided communication window
! into a local array
! it is an error to attempt a "remote" get when the window is not "exposed"
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! larray (IN)     C compatible pointer (type(C_PTR)) to local array (destination of get)
! target (IN)     ordinal in window communicator of remote PE
! offset (IN)     displacement (origin 0) into remote PE window data array
! nelem  (IN)     number of elements to transfer (type of element was defined at window creation)
! ierr   (OUT)    error status, RPN_COMM_OK or RPN_COMM_ERROR
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  type(C_PTR), intent(IN), value :: larray                            !InTf!
  integer, intent(IN) :: target                                       !InTf!
  integer, intent(IN) :: offset                                       !InTf!
  integer, intent(IN) :: nelem                                        !InTf!
!******

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

!InTf!
!****f* rpn_comm/RPN_COMM_i_win_get_l read a local one sided communication window
! SYNOPSIS
subroutine RPN_COMM_i_win_get_l(window,larray,offset,nelem,ierr)      !InTf!
!===============================================================================
! one sided communication local get (read) from one sided communication window
! into a local array
! it is an error to attempt a "local" get when the window is "exposed"
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! larray (IN)     C compatible pointer (type(C_PTR)) to local array (destination of get)
! offset (IN)     displacement (origin 0) into this PE window data array
! nelem  (IN)     number of elements to transfer (type of element was defined at window creation)
! ierr   (OUT)    error status, RPN_COMM_OK or RPN_COMM_ERROR
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  type(C_PTR), intent(IN), value :: larray                            !InTf!
  integer, intent(IN) :: offset                                       !InTf!
  integer, intent(IN) :: nelem                                        !InTf!
!******

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
