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
#define IN_RPN_COMM_window
#if defined(WITH_DOC)
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
!       (this array will either be supplied by the caller or allocated internally.
!
!   all further operations refer to the window by its "identifier" 
!   Fortran type : type(rpncomm_window)  (include 'RPN_COMM_types.inc')
!
!   the window creation is a COLLECTIVE operation, all members of the communicator group
!   must call RPN_COMM_i_win_create
!
!   "exposing" a window and terminating the "exposure" of a window are also a COLLECTIVE operation
!
!   remote get/put operations, i.e. sending read/write requests targeting the communication
!   window belonging to a remote PE may only happen when the window is "exposed"
!   remote operations are NOT ALLOWED when the window is "not exposed"
!
!   two modes of one sided communication are available.
!   - active mode
!      a) for the whole communicator group  (fence)
!      b) for a subset of the communicator group (scales better when said group is large)
!         (start/complete/post/wait)
!   - passive mode
!   this mode is selected when calling the routine that "exposes" the communication window
!
!   local get/put operations, i.e. reading/writing from/into the window
!   belonging the local PE may only happen when the window is "not exposed"
!   local operations are NOT ALLOWED whe the window is "exposed" as the result of
!   such an operation would be "undefined"
!
!   a typical sequence of operations would be
!   1 - create a one sided communication window (COLLECTIVE operation)
!    repeat as needed
!    {
!      2a - modify the local copy of the window (if and as needed)
!      2b - "expose" the window and select active or passive mode (COLLECTIVE operation)
!      2c - perform get/put operations targeting remote PEs (as needed on each PE)
!      2d - "end exposing" the window (COLLECTIVE operation)
!      2e - get modifications from the local copy of the window (if and as needed)
!    }
!   3 - free the one sided communication window (COLLECTIVE operation)
!
!   window creation/destruction : RPN_COMM_i_win_create, RPN_COMM_i_win_create_secondary, RPN_COMM_i_win_free
!   window status queries       : RPN_COMM_i_valid_win, RPN_COMM_i_win_check
!   window info queries         : RPN_COMM_i_win_get_ptr, RPN_COMM_i_win_get_size
!   window "exposition"         : RPN_COMM_i_win_open, RPN_COMM_i_win_close
!   window get operations       : RPN_COMM_i_win_get_r, RPN_COMM_i_win_get_l
!   window put operations       : RPN_COMM_i_win_put_r, RPN_COMM_i_win_put_l
!
!  -------------- one sided request/reply system (RPN_COMM_i_win_post) --------------
!
!  in this special case, two (2) MPI one sided windows are used in internal type rpncomm_windef
!
!  the primary window MUST be created by a call to RPN_COMM_i_win_create
!    the storage associated with the primary window may be allocated by the user or internally.
!
!  the secondary window MUST be created by a call to RPN_COMM_i_win_create_secondary after creation
!    of rpncomm_window by a previous call to RPN_COMM_i_win_create.
!  (the storage associated with the secondary window  will be allocated internally)
!
!  primary window layout for request/reply  on PE with ordinal PEK  (client)
!
!  +-------------+-------------------+-------------------------------------+-------------------+--------------------
!  |             |  request to PEJ   |                                     |  reply from PEJ   |
!  +-------------+-------------------+-------------------------------------+-------------------+--------------------
!                ^                                                         ^
!        offseti |                                                 offseto |
!                <--    n_in       -->                                     <--    n_out     -->
!
!
!  secondary window layout for request/reply on PE with ordinal PEJ  (supplier)
!
!  +-----------------+-----------------+.............+-----------------+.............+-----------------+
!  | request 1       |  request 2      |             |  request n      |             +  request MAXREQ |
!  +-----------------+-----------------+.............+-----------------+.............+-----------------+
!
!  request format : 5 integers (PEK, offseti, n_in, offseto, n_out)  (request from PEK)
!  request 1 is special : (last valid request, 0, 0, 0, 0) (initial value of last valid request is 1)
!  PE ordinals are ranks in the communicator associated with both primary and secondary windows
!
!  PEJ is expected to "get" n_in words at displacement offseti in primary window of PEK
!  PEJ is expected to "put" n_out words at displacement offseto in primary window of PEK
!  a "word" may be anything having the size of MPI_INTEGER (integer, real, ...)
!
!  typical sequence of operations for request/reply:
!
!  1 - call RPN_COMM_i_win_open(some_window,.false.,ierr)
!  2 - call RPN_COMM_i_win_post(some_window,pe(:),offseti(:),nelemi(:),offseto(:),nelemo(:),nreq,ierr)
!      (step 2 may be repeated)
!  3 - call RPN_COMM_i_win_close(some_window,ierr)
!  4 - call RPN_COMM_i_win_open(some_window,.false.,ierr)
!  5 - remote "get" operations to get request(s)
!  6 - remote "put" operations to put replies
!  7 - call RPN_COMM_i_win_close(some_window,ierr)
!
!  RPN_COMM_i_win_open and RPN_COMM_i_win_close are COLLECTIVE operations
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
! (used internally by user callable routines)
!===============================================================================
  use ISO_C_BINDING
  implicit none
  include 'mpif.h'
  include 'RPN_COMM_types.inc'
  include 'RPN_COMM_types_int.inc'
  include 'RPN_COMM_constants.inc'
!******
  integer, parameter :: RPN_COMM_MAX_WINDOWS = 64
  integer, save :: integer_size = 0                                       ! will contain extent of MPI_INTEGER
  type(rpncomm_windef), dimension(:), pointer, save :: win_tab => NULL()  ! rpn comm window table
  logical, save :: debug_mode = .false.
  type(rpncomm_operator), save :: op = NULL_rpncomm_operator
  type :: slot
    integer :: pe  ! PE on which this PE has acquired a message slot
    integer :: n   ! slot number
  end type slot
  integer, parameter :: WINDEF_SLOT_SIZE = 2     ! number of integers in type slot
   interface
    function c_malloc(siz) result (p) bind(C,name='malloc')
      import C_SIZE_T, C_PTR
      integer(C_SIZE_T), intent(IN), value :: siz
      type(C_PTR) :: p
    end function c_malloc
    subroutine c_free(p) bind(C,name='free')
      import C_PTR
      type(C_PTR), intent(IN), value :: p
    end subroutine c_free
  end interface
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

    if(ASSOCIATED(win_tab)) return                  ! already initialized, nothing to do
    allocate(win_tab(RPN_COMM_MAX_WINDOWS))         ! allocate for up to RPN_COMM_MAX_WINDOWS windows
    do i = 1 , RPN_COMM_MAX_WINDOWS
      win_tab(i) = NULL_rpncomm_windef              ! null entry (see RPN_COMM_types_int.inc)
      win_tab(i)%com = MPI_COMM_NULL                ! with MPI null communicator
      win_tab(i)%typ = MPI_DATATYPE_NULL            ! and MPI null datatype
      win_tab(i)%win = MPI_WIN_NULL                 ! MPI null window
    enddo
    call MPI_TYPE_EXTENT(MPI_INTEGER, integer_size, i) ! get size(extent) of MPI_INTEGER
  end subroutine init_win_tab

!****if* RPN_COMM_windows/win_ptr
! SYNOPSIS
  function valid_win_entry(win_ptr) result(is_valid)
!===============================================================================
! check if win_ptr (C pointer pointing into win_tab) points to a valid entry
!===============================================================================
!******
    implicit none
    type(C_PTR), intent(IN), value :: win_ptr        ! pointer to entry in win_tab
    logical :: is_valid

    type(rpncomm_windef), pointer :: win_entry
    type(C_PTR) :: temp

    is_valid = .false.
    if(.not. c_associated(win_ptr)) return           ! win_ptr C pointer is not valid
    if(.not. associated(win_tab)) return             ! win_tab does not exist yet

    call c_f_pointer(win_ptr,win_entry)              ! make Fortran pointer to win_entry from win_ptr

    if(win_entry%indx <=0 .or. win_entry%indx>RPN_COMM_MAX_WINDOWS) return  ! invalid index into window table
    temp = c_loc(win_tab(win_entry%indx))            ! address pointed to by win_entry index component
    if(.not. c_associated(temp,win_ptr)) return      ! component index should point to win_entry itself
    if(.not. c_associated(win_entry%base)) return    ! no array associated with entry
    if(win_entry%com == MPI_COMM_NULL) return        ! no valid communicator
    if(win_entry%typ == MPI_DATATYPE_NULL) return    ! no valid datatype
    if(win_entry%win == MPI_WIN_NULL) return         ! no valid window
    if(win_entry%siz <= 0) return                    ! invalid size (MUST be > 0)

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
    integer, intent(IN) :: typ             ! MPI data type (as defined in mpif.h)
    integer, intent(IN) :: siz             ! number of elements in array associated to this window
    integer, intent(IN) :: comm            ! MPI communicator for this one sided window
    integer, intent(OUT) :: indx           ! index into wintab of created window (used for consistency test)
    integer, intent(OUT) :: ierr           ! return status (RPN_COMM_ERROR or RPN_COMM_OK)

    integer :: i, extent, ierror, group
    integer(kind=MPI_ADDRESS_KIND) win_size  ! may be wider than a default integer
    integer, dimension(:), pointer :: remote   ! base translated as Fortran pointer to primary array

    if(.not. associated(win_tab)) call init_win_tab  ! create and initialize window table if necessary
    ierr = RPN_COMM_ERROR                  ! preset for failure
    indx = -1

    call MPI_comm_group(comm,group,ierror)
    if(ierror .ne. MPI_SUCCESS) return     ! error getting group associated with communicator (bad communicator ?)

    call MPI_TYPE_EXTENT(typ, extent, ierror)  ! determine size associated with MPI datatype
    if(ierror .ne. MPI_SUCCESS) return         !  invalid type ? other error ?
    if( mod(extent,integer_size) .ne. 0 ) return    ! extent of data type must be a multiple of integer size

    win_size = extent * siz    ! possibly wider integer needed by MPI_win_create

    do i = 1 , RPN_COMM_MAX_WINDOWS        ! loop over window table entries (we are looking for a free entry)
      if(.not. c_associated(win_tab(i)%base) ) then  ! this entry is free
        if(c_associated(base)) then    ! if user has already allocated space
          if(debug_mode) print *,'DEBUG: user supplied space, size=',win_size
          win_tab(i)%base = base
          win_tab(i)%is_user = .true.  ! user supplied space
        else                           ! otherwise allocate needed space with MPI_alloc_mem
          if(debug_mode) print *,'DEBUG: allocating with MPI_alloc_mem, size=',win_size
          call MPI_alloc_mem(win_size, MPI_INFO_NULL, win_tab(i)%base,ierr)
          if(ierr .ne. MPI_SUCCESS) return       ! could not allocate memory, OUCH !!
          win_tab(i)%is_user = .false.  ! internally supplied space
        endif
        win_tab(i)%base2 = C_NULL_PTR
        win_tab(i)%nbase2 = 0
        win_tab(i)%slots = C_NULL_PTR
        win_tab(i)%nslots = 0
        win_tab(i)%win2 = MPI_WIN_NULL
!         call c_f_pointer(win_tab(i)%base,win_tab(i)%remote,[extent*siz])  ! make Fortran pointer to array associated with window
        call c_f_pointer(win_tab(i)%base,remote,[extent*siz])  ! make Fortran pointer to array associated with window
!         win_tab(i)%remote2 => NULL()

!         call MPI_win_create(win_tab(i)%remote, win_size, extent, MPI_INFO_NULL, comm, win_tab(i)%win, ierror) ! create window
        call MPI_win_create(remote, win_size, extent, MPI_INFO_NULL, comm, win_tab(i)%win, ierror) ! create window
        if(ierror .ne. MPI_SUCCESS) then
          print *,'ERROR: window creation unsuccessful'
          return
        endif
        indx = i
        win_tab(i)%nbase = extent*siz
        win_tab(i)%is_open = .false.             ! window is not "exposed"
        win_tab(i)%opr = MPI_OP_NULL             ! MPI accumulate operator associated to window
        win_tab(i)%typ = typ                     ! MPI data type associated to window
        win_tab(i)%ext = extent/integer_size     ! size (extent) of MPI data type in MPI_INTEGER units
        win_tab(i)%indx = indx                   ! entry index points to itself
        win_tab(i)%siz = siz                     ! number of elements in window
        win_tab(i)%com = comm                    ! communicator associated with window
        call mpi_comm_size(comm,win_tab(i)%npe,ierror)   ! number of PEs in communicator
        call mpi_comm_rank(comm,win_tab(i)%rank,ierror)  ! my rank in communicator
        win_tab(i)%grp = group                   ! group associated with said communicator
        win_tab(i)%active_mode = .true.          ! use active mode by default
        win_tab(i)%s_group = MPI_GROUP_NULL      ! no group by default for remote targets
        win_tab(i)%r_group = MPI_GROUP_NULL      ! no group by default for exposure to remote accesses
        ierr = RPN_COMM_OK                       ! all is well
        return
      endif
    enddo
    return  ! if we fall through the loop, the table was full, we will return RPN_COMM_ERROR
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

  if( ieor(window%t1,RPN_COMM_MAGIC) .ne. window%t2 )    return  ! inconsistent tags
  if(window%t2 <= 0 .or. window%t2>RPN_COMM_MAX_WINDOWS) return  ! invalid index

  temp = c_loc(win_tab(window%t2))                 ! address of relevant win_tab entry
  if( .not. c_associated(window%p,temp)) return    ! window%p must point to win_tab(window%t2)

  is_valid = valid_win_entry(window%p)             ! check that window description itself in table is valid
  return

end function win_valid

end module RPN_COMM_windows

!===============================================================================
! test code for one sided request/reply window package
!===============================================================================
!****f* rpn_comm/DEMO demo code for one sided request/reply
! SYNOPSIS
subroutine RPN_COMM_i_win_req_test(nparams,params)
  use ISO_C_BINDING
  implicit none
  include 'RPN_COMM.inc'
  include 'mpif.h'
  integer, intent(IN) :: nparams
  integer, intent(IN), dimension(nparams) :: params

  integer, parameter :: REQUEST_SIZE = 2
  real, dimension(REQUEST_SIZE), target :: request
  integer, parameter :: REPLY_SIZE = 3
  real, dimension(REPLY_SIZE), target :: reply
  real :: a, b
  integer :: status, Me, Me_x, Me_y, npes, ierr, i, indx, nreq, nreq2, errors, toterrors
  type(rpncomm_window) :: window
  type(rpncomm_datatype) :: dtype
  type(rpncomm_communicator) :: com
  integer, dimension(:), pointer :: pe, offseti, n_in, offseto, n_out
  type(C_PTR) :: cptr
  real, dimension(:), pointer :: request_reply1, request_reply2
  integer, dimension(:,:), pointer :: requests
  integer :: pe_from, pe_offset, pe_nwds 
  real *8 t1, t2

  status = RPN_COMM_mype(Me,Me_x,Me_y)
  call RPN_COMM_size( RPN_COMM_GRID, npes ,ierr )
  call RPN_COMM_i_comm(RPN_COMM_GRID,com)
  call RPN_COMM_i_datyp(RPN_COMM_INTEGER,dtype)

  allocate(pe(npes),offseti(npes), n_in(npes), offseto(npes), n_out(npes))
  nreq = 0
  if(npes < 4) print *,'requests to PE(n) from PE',Me
  do i = 0 , npes-1
    if(i .ne. Me) then
      nreq = nreq + 1
      pe(nreq) = i
      offseti(nreq) = (nreq - 1) * REQUEST_SIZE
      n_in(nreq) = REQUEST_SIZE
      offseto(nreq) = npes*REQUEST_SIZE + (nreq - 1) * REPLY_SIZE
      n_out(nreq) = REPLY_SIZE
      if(npes < 4) print 100,nreq,pe(nreq),offseti(nreq),n_in(nreq),offseto(nreq),n_out(nreq)
    endif
  enddo

  allocate(request_reply1(npes * (REPLY_SIZE + REQUEST_SIZE)))
  allocate(request_reply2(npes * (REPLY_SIZE + REQUEST_SIZE)))
  request_reply1 = 999.999

  do i = 1 , nreq                    ! requests
    request_reply1(2*i-1) = i+5*Me
    request_reply1(2*i) = i+1+5*Me
  enddo

  request_reply2 = request_reply1
  do i = 1 , nreq
    a = request_reply1(2*i-1)         ! requests
    b = request_reply1(2*i)
    request_reply2(offseto(i) + 1) = a + b       ! expected replies
    request_reply2(offseto(i) + 2) = a * b
    request_reply2(offseto(i) + 3)   = a*a + b*b
  enddo
  if(npes < 4) print 101,request_reply1
101 format(15F8.1)

  call RPN_COMM_i_win_create(window, dtype, size(request_reply1), com, C_LOC(request_reply1(1)), ierr)
  call RPN_COMM_i_win_create_secondary(window,npes-1,ierr)   ! deliberately wrong allocation size
  call RPN_COMM_i_win_create_secondary(window,-npes,ierr)    ! reallocate with correct size

  cptr = RPN_COMM_i_win_get_ptr(window,2,ierr)    ! get pointer to secondary window
  call c_f_pointer(cptr,requests,[5,npes])        ! fortran pointer to table of client requests
  t1 = MPI_wtime()
  call RPN_COMM_i_win_open(window,.false.,ierr)   ! open in passive mode
  call RPN_COMM_i_win_post(window,pe,offseti,n_in,offseto,n_out,nreq,ierr)  ! post requests
  call RPN_COMM_i_win_close(window,ierr)          ! close window
  t1 = MPI_wtime() - t1

  nreq2 = requests(1,1)
  print *,'number of requests received = ',nreq2-1
  do i = 2, npes
    if(npes < 4) print 100,i,requests(:,i)
100 format(I3,15I6)
  enddo

  t2 = MPI_wtime()
  do i = 2 , nreq2
    pe_from   = requests(1,i)
    pe_offset = requests(2,i)
    pe_nwds   = requests(3,i)
    call RPN_COMM_i_win_open(window,.false.,ierr)   ! open in passive mode
    call RPN_COMM_i_win_get_r(window,C_LOC(request),pe_from,pe_offset,pe_nwds,ierr)  ! get request (a,b)
    call RPN_COMM_i_win_close(window,ierr)          ! close window
    reply(1) = request(1) + request(2)     ! a + b
    reply(2) = request(1) * request(2)     ! a * b
    reply(3) = request(1)*request(1) + request(2)*request(2)  ! a*a + b*b
    pe_offset = requests(4,i)
    pe_nwds   = requests(5,i)
    if(npes < 4) print *,'request+reply :',pe_from,pe_offset,pe_nwds
    if(npes < 4) print 101,request,reply
    call RPN_COMM_i_win_open(window,.false.,ierr)   ! open in passive mode
    call RPN_COMM_i_win_put_r(window,C_LOC(reply),pe_from,pe_offset,pe_nwds,ierr)  ! put reply (a+b , a*b, a*a+b*b)
    call RPN_COMM_i_win_close(window,ierr)          ! close window
  enddo
  t2 = MPI_wtime() - t2

  indx = 2*npes
  errors = 0
  do i = indx+1 , indx+3*nreq   ! check replies against what was expected
    if( request_reply2(i) .ne. request_reply1(i) ) errors = errors + 1
  enddo
  call MPI_allreduce(errors,toterrors,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  print *,"local errors, total errors =",errors,toterrors,indx+1 , indx+3*nreq
  print 102,"time to post requests=",t1," time to reply=",t2
102 format(A,F9.6,A,F9.6)
  if(npes < 4) print 101,request_reply1
  if(npes < 4) print 101,request_reply2

999 continue
  deallocate(pe, offseti, n_in, offseto, n_out, request_reply1, request_reply2)
  print *,"freeing window"
  call RPN_COMM_i_win_free(window,ierr)
  return
end subroutine RPN_COMM_i_win_req_test
!******
!===============================================================================
! test code for one sided communication window package
!===============================================================================
subroutine RPN_COMM_i_win_test(nparams,params)
  use ISO_C_BINDING
  implicit none
  include 'RPN_COMM.inc'
  include 'mpif.h'
  integer, intent(IN) :: nparams
  integer, intent(IN), dimension(nparams) :: params

  integer, parameter :: DATA_SIZE = 1024
  type(rpncomm_window) :: window, window2
  type(rpncomm_datatype) :: dtype
  type(rpncomm_communicator) :: com
  type(rpncomm_operator) :: oper
  integer :: siz
  type(C_PTR) :: array, array2, cptr, cptr2
  integer, dimension(:), pointer :: fptr, fptr2
!   integer, dimension(DATA_SIZE), target :: my_data
  integer, dimension(:), pointer :: my_data
  integer(kind=MPI_ADDRESS_KIND) my_win_size  ! may be wider than a default integer
  type(C_PTR) :: my_win_base
  integer :: ierr, i, nerrors, nval, offset, offset_r, expected, got
  integer :: me, me_x, me_y, status, wsiz, npes, target_pe, from_pe
  integer, dimension(100), target :: local, local2
  integer, dimension(:), pointer :: errors

!   debug_mode = .true.
!
  my_win_size = DATA_SIZE * 4   ! size in bytes of window
!   call MPI_Alloc_mem(my_win_size, MPI_INFO_NULL, my_win_base,ierr)   ! use MPI library allocator
!   if(ierr .ne. MPI_SUCCESS) goto 888
!   call c_f_pointer(my_win_base,my_data,[DATA_SIZE])
  allocate(my_data(DATA_SIZE))                                       ! use regular memory
!
  siz = DATA_SIZE
  status = RPN_COMM_mype(Me,Me_x,Me_y)
  call RPN_COMM_size( RPN_COMM_GRID, npes ,ierr )
  print *,'TEST INFO: this is PE',me+1,' of',npes
  allocate(errors(npes))
  errors = 0
  call RPN_COMM_i_comm(RPN_COMM_GRID,com)
  call RPN_COMM_i_datyp(RPN_COMM_INTEGER,dtype)

  do i = 1 , siz
    my_data(i) = -i
  enddo
! ================================== create windows ======================================
  array = c_loc(my_data(1))
  call RPN_COMM_i_win_create(window,dtype,siz,com,array,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: created window, PE=',me
  else
    print *,'TEST ERROR: cannot create window, PE=',me
    goto 888
  endif

  array2 = C_NULL_PTR
  call RPN_COMM_i_win_create(window2,dtype,siz*2,com,array2,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: created window2, PE=',me
  else
    print *,'TEST ERROR: cannot create window2, PE=',me
    goto 888
  endif
! ================================== get windows properties ======================================
  wsiz = RPN_COMM_i_win_get_size(window,ierr)  ! get size of data array associated to window
  if(me == 0) then
    print *,'TEST INFO: size associated to window =',wsiz
  endif
  cptr = RPN_COMM_i_win_get_ptr(window,1,ierr)   ! get C pointer to primary local window data
  call c_f_pointer(cptr,fptr,[wsiz])           ! convert to Fortran pointer
  nerrors = 0
  do i = 1 , wsiz                              ! check values in array (consistency check)
    if(fptr(i) .ne. my_data(i)) nerrors = nerrors + 1
  enddo
  print *,'TEST INFO: nerrors for window =',nerrors
  if(nerrors > 0) goto 888

  wsiz = RPN_COMM_i_win_get_size(window2,ierr)  ! get size of data array associated to window
  if(me == 0) then
    print *,'TEST INFO: size associated to window2 =',wsiz
  endif
  cptr2 = RPN_COMM_i_win_get_ptr(window2,1,ierr)   ! get C pointer to primary local window data
  call c_f_pointer(cptr2,fptr2,[wsiz])           ! convert to Fortran pointer
  fptr2 = -1                                     ! initialize contents of array associated with windows2
! ===================================== remote PUT test ===================================
  target_pe = mod(me+1,npes)       ! remote target
  from_pe = mod(me+npes-1,npes)    ! remote accessor
  print *,'TEST INFO: target PE is',target_pe,' accessing PE is',from_pe
! ===================================== active mode ===================================
  do i = 1 , siz
    my_data(i) = -i
  enddo
  nval = 5                         ! number of values for remote put
  offset = 10                      ! displacement into remote window
! ==================================
goto 600
!  if(me == 0) then
    print *,'=============================================================================='
    print *,'TEST INFO: PUT active fence mode test start'
!  endif

  call RPN_COMM_i_win_open(window,.true.,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: accessing/exposing window, PE=',me
  else
    print *,'TEST ERROR: cannot access/expose window, PE=',me
    goto 888
  endif

  local(1:5) = [2,4,6,8,10] + 100*me
  call RPN_COMM_i_win_put_r(window,c_loc(local(1)),target_pe,offset,nint(nval/2.0),ierr)  ! (window,larray,target,offset,nelem,ierr)
  do i = nint(nval/2.0),nval
    call RPN_COMM_i_win_put_r(window,c_loc(local(i)),target_pe,offset+i-1,1,ierr)
  enddo

  call RPN_COMM_i_win_close(window,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: ending access/exposure to window, PE=',me
  else
    print *,'TEST ERROR: cannot end access/exposure to window, PE=',me
    goto 888
  endif

  nerrors = 0
  local(1:5) = local(1:5) + 100*from_pe - 100*me
  do i = 1 , size(my_data)
    if(i>offset .and. i<=offset+nval) then
      if( my_data(i) .ne. local(i-offset)) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',my_data(i),' expected',local(i-offset),' error'
      endif
    else
      if(my_data(i) .ne. -i) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',my_data(i),' expected',-i,' error'
      endif
    endif
  enddo
  errors = 0
  call RPN_COMM_allgather(nerrors,1,RPN_COMM_INTEGER,errors,1,RPN_COMM_INTEGER,RPN_COMM_GRID,ierr)
  print *,'TEST INFO:',sum(errors),' unexpected values found from PUT',errors
  if(sum(errors) > 0) goto 777
!  if(nerrors .ne. 0) then
!    print *,'TEST INFO:',nerrors),' unexpected values found from PUT'
!  endif

!  if(me == 0) then
    print *,'TEST INFO: PUT active fence mode test end'
!  endif
! ===================================== active group mode ===================================
  do i = 1 , siz
    my_data(i) = -i
  enddo
  nval = 5                         ! number of values for remote put
  offset = 35                      ! displacement into remote window
! ==================================
!  if(me == 0) then
    print *,'=============================================================================='
    print *,'TEST INFO: PUT active group mode test start'
!  endif

  call RPN_COMM_i_win_group(window,[target_pe],[from_pe],ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: group creation OK'
  else
    print *,'TEST ERROR: group creation failure'
    goto 888
  endif
  call RPN_COMM_i_win_open(window,.true.,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: accessing/exposing window, PE=',me
  else
    print *,'TEST ERROR: cannot access/expose window, PE=',me
    goto 888
  endif

  local(1:5) = [12,14,16,18,20] + 100*me
  call RPN_COMM_i_win_put_r(window,c_loc(local(1)),target_pe,offset,nint(nval/2.0),ierr)  ! (window,larray,target,offset,nelem,ierr)
  do i = nint(nval/2.0),nval
    call RPN_COMM_i_win_put_r(window,c_loc(local(i)),target_pe,offset+i-1,1,ierr)
  enddo

  call RPN_COMM_i_win_close(window,ierr)

  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: ending access/exposure to window, PE=',me
  else
    print *,'TEST ERROR: cannot end access/exposure to window, PE=',me
    goto 888
  endif

  nerrors = 0
  local(1:5) = local(1:5) + 100*from_pe - 100*me
  do i = 1 , size(my_data)
    if(i>offset .and. i<=offset+nval) then
      if( my_data(i) .ne. local(i-offset)) then
        nerrors = nerrors + 1
        print *,'at i=',i,' got',my_data(i),' expected',local(i-offset),' error'
      endif
    else
      if(my_data(i) .ne. -i) then
        nerrors = nerrors + 1
        print *,'at i=',i,' got',my_data(i),' expected',-i,' error'
      endif
    endif
  enddo
  errors = 0
  call RPN_COMM_allgather(nerrors,1,RPN_COMM_INTEGER,errors,1,RPN_COMM_INTEGER,RPN_COMM_GRID,ierr)
  print *,'TEST INFO:',sum(errors),' unexpected values found from PUT',errors
  if(sum(errors) > 0) goto 777
!   if(nerrors .ne. 0) then
!     print *,'TEST ERROR:',nerrors,' unexpected values found from PUT'
!   endif

  call RPN_COMM_i_win_group(window,[-1],[-1],ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: window group deletion OK'
  else
    print *,'TEST ERROR: window group deletion failure'
    goto 888
  endif

!  if(me == 0) then
    print *,'TEST INFO: PUT active group mode test end'
!  endif
! ================================== passive mode ======================================
600 continue
  do i = 1 , siz
    my_data(i) = -i
  enddo
  nval = 5                         ! number of values for remote put
  offset = 17                      ! displacement into remote window
! ==================================
!  if(me == 0) then
    print *,'=============================================================================='
    print *,'TEST INFO: PUT passive mode test start'
!  endif

  target_pe = mod(me+1,npes)

  call RPN_COMM_i_win_open(window,.false.,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: accessing/exposing window, PE=',me
  else
    print *,'TEST ERROR: cannot access/expose window, PE=',me
    goto 888
  endif

  local(1:5) = [1,3,5,7,9] + 100*me
  nerrors = 0
  call RPN_COMM_i_win_put_r(window,c_loc(local(1)),target_pe,offset,nint(nval/2.0),ierr)  ! (window,larray,target,offset,nelem,ierr)
  if(ierr .ne. RPN_COMM_OK) nerrors = nerrors + 1
! !   call RPN_COMM_i_win_put_r(window,c_loc(local(1)),target_pe,offset,nint(nval/2.0),ierr)  ! (window,larray,target,offset,nelem,ierr)
!   if(ierr .ne. RPN_COMM_OK) nerrors = nerrors + 1
  do i = nint(nval/2.0),nval
    call RPN_COMM_i_win_put_r(window,c_loc(local(i)),target_pe,offset+i-1,1,ierr)
!     call RPN_COMM_i_win_put_r(window,c_loc(local(i)),target_pe,offset+i-1,nval-i+1,ierr)
  if(ierr .ne. RPN_COMM_OK) nerrors = nerrors + 1
  enddo

  call RPN_COMM_i_win_close(window,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: ending access/exposure to window, PE=',me
  else
    print *,'TEST ERROR: cannot end access/exposure to window, PE=',me
    goto 888
  endif

  local(1:5) = local(1:5) + 100*from_pe - 100*me
  do i = 1 , size(my_data)
    if(i>offset .and. i<=offset+nval) then
      got = my_data(i)
      expected = local(i-offset)
      if( got .ne. expected)  then
        nerrors = nerrors + 1
        print *,'i=',i,' got',got,' expected',expected,' ERROR',nerrors
      else
        print *,'i=',i,' got',got,' expected',expected,' OK'
      endif
    else
      if(my_data(i) .ne. -i) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',my_data(i),' expected',-i,' ERROR',nerrors
      endif
    endif
  enddo
  errors = 0
  call RPN_COMM_allgather(nerrors,1,RPN_COMM_INTEGER,errors,1,RPN_COMM_INTEGER,RPN_COMM_GRID,ierr)
  print *,'TEST INFO:',sum(errors),' unexpected values found from PUT',errors
  if(sum(errors) > 0) goto 777
!   if(nerrors .ne. 0) then
!     print *,'TEST ERROR:',nerrors,' unexpected values found from PUT'
!   endif

!  if(me == 0) then
    print *,'TEST INFO: PUT passive mode test end'
!  endif
! ===================================== remote GET test ===================================
  target_pe = mod(me+1,npes)       ! remote target
  from_pe = mod(me+npes-1,npes)    ! remote accessor
goto 601
! ===================================== active mode ===================================
  do i = 1 , siz
    my_data(i) = -i
  enddo
  nval = 5                         ! number of values for remote put
  offset = 10                      ! displacement into remote window
! ==================================
!  if(me == 0) then
    print *,'=============================================================================='
    print *,'TEST INFO: GET active fence mode test start'
!  endif

!  local = -(me+1)
  local(1:5) = [2,4,6,8,10]
  do i = 1,nval
    my_data(i+offset) = local(i) + 100*me
  enddo
  local2 = -1
  print *,'TEST INFO: my_data',my_data(offset+1:offset+nval)

  call RPN_COMM_i_win_open(window,.true.,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: accessing/exposing window, PE=',me
  else
    print *,'TEST ERROR: cannot access/expose window, PE=',me
    goto 888
  endif

  call RPN_COMM_i_win_get_r(window,c_loc(local2(1)),target_pe,offset,nint(nval/2.0),ierr)  ! (window,larray,target,offset,nelem,ierr)
  do i = nint(nval/2.0)+1,nval
    call RPN_COMM_i_win_get_r(window,c_loc(local2(i)),target_pe,offset+i-1,1,ierr)
  enddo

  call RPN_COMM_i_win_close(window,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: ending access/exposure to window, PE=',me
  else
    print *,'TEST ERROR: cannot end access/exposure to window, PE=',me
    goto 888
  endif

  nerrors = 0
  do i = 1 , nval
      expected = my_data(i+offset)-100*me+100*target_pe
      got = local2(i)
      if( got .ne. expected ) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',got,' expected',expected,' ERROR',nerrors
      else
        print *,'i=',i,' got',got,' expected',expected,' OK'
      endif
  enddo
  errors = 0
  call RPN_COMM_allgather(nerrors,1,RPN_COMM_INTEGER,errors,1,RPN_COMM_INTEGER,RPN_COMM_GRID,ierr)
  print *,'TEST INFO:',sum(errors),' unexpected values found from GET',errors
  if(sum(errors) > 0) goto 777
!   if(nerrors .ne. 0) then
!     print *,'TEST INFO:',nerrors,' local unexpected values found from GET'
!   endif

!   if(me == 0) then
    print *,'TEST INFO: GET active fence mode test end'
!   endif
! ===================================== active group mode ===================================
  do i = 1 , siz
    my_data(i) = -i
  enddo
  nval = 5                         ! number of values for remote put
  offset = 10                      ! displacement into remote window
! ==================================
!  if(me == 0) then
    print *,'=============================================================================='
    print *,'TEST INFO: GET active group mode test start'
!  endif

!  local = -(me+1)
  local(1:5) = [12,14,16,18,20]
  do i = 1,nval
    my_data(i+offset) = local(i) + 100*me
  enddo
  local2 = -1
  print *,'TEST INFO: my_data',my_data(offset+1:offset+nval)

  call RPN_COMM_i_win_group(window,[target_pe],[from_pe],ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: group creation OK'
  else
    print *,'TEST ERROR: group creation failure'
    goto 888
  endif
  call RPN_COMM_i_win_open(window,.true.,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: accessing/exposing window, PE=',me
  else
    print *,'TEST ERROR: cannot access/expose window, PE=',me
    goto 888
  endif

  call RPN_COMM_i_win_get_r(window,c_loc(local2(1)),target_pe,offset,nint(nval/2.0),ierr)  ! (window,larray,target,offset,nelem,ierr)
  do i = nint(nval/2.0)+1,nval
    call RPN_COMM_i_win_get_r(window,c_loc(local2(i)),target_pe,offset+i-1,1,ierr)
  enddo

  call RPN_COMM_i_win_close(window,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: ending access/exposure to window, PE=',me
  else
    print *,'TEST ERROR: cannot end access/exposure to window, PE=',me
    goto 888
  endif

  nerrors = 0
  do i = 1 , nval
      expected = my_data(i+offset)-100*me+100*target_pe
      got = local2(i)
      if( got .ne. expected ) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',got,' expected',expected,' ERROR',nerrors
      else
        print *,'i=',i,' got',got,' expected',expected,' OK'
      endif
  enddo
  errors = 0
  call RPN_COMM_allgather(nerrors,1,RPN_COMM_INTEGER,errors,1,RPN_COMM_INTEGER,RPN_COMM_GRID,ierr)
  print *,'TEST INFO:',sum(errors),' unexpected values found from GET',errors
  if(sum(errors) > 0) goto 777
!   if(nerrors .ne. 0) then
!     print *,'TEST INFO:',nerrors,' local unexpected values found from GET'
!   endif

  call RPN_COMM_i_win_group(window,[-1],[-1],ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: window group deletion OK'
  else
    print *,'TEST ERROR: window group deletion failure'
    goto 888
  endif

!   if(me == 0) then
    print *,'TEST INFO: GET active group mode test end'
!   endif
! ===================================== passive mode ===================================
601 continue
  do i = 1 , siz
    my_data(i) = -i
  enddo
  nval = 5                         ! number of values for remote put
  offset = 10                      ! displacement into remote window
! ==================================
!  if(me == 0) then
    print *,'=============================================================================='
    print *,'TEST INFO: GET passive mode test start'
!  endif

!  local = -(me+1)
  local(1:5) = [1,3,5,7,9]
  do i = 1,nval
    my_data(i+offset) = local(i) + 100*me
  enddo
  local2 = -1
  print *,'TEST INFO: my_data',my_data(offset+1:offset+nval)

  call RPN_COMM_i_win_open(window,.false.,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: accessing/exposing window, PE=',me
  else
    print *,'TEST ERROR: cannot access/expose window, PE=',me
    goto 888
  endif

  nerrors = 0
  call RPN_COMM_i_win_get_r(window,c_loc(local2(1)),target_pe,offset,nint(nval/2.0),ierr)  ! (window,larray,target,offset,nelem,ierr)
  if(ierr .ne. RPN_COMM_OK) nerrors = nerrors + 1
  do i = nint(nval/2.0)+1,nval
    call RPN_COMM_i_win_get_r(window,c_loc(local2(i)),target_pe,offset+i-1,1,ierr)
    if(ierr .ne. RPN_COMM_OK) nerrors = nerrors + 1
  enddo

  call RPN_COMM_i_win_close(window,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: ending access/exposure to window, PE=',me
  else
    print *,'TEST ERROR: cannot end access/exposure to window, PE=',me
    goto 888
  endif

  do i = 1 , nval
      expected = my_data(i+offset)-100*me+100*target_pe
      got = local2(i)
      if( got .ne. expected ) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',got,' expected',expected,' ERROR',nerrors
      else
        print *,'i=',i,' got',got,' expected',expected,' OK'
      endif
  enddo
  errors = 0
  call RPN_COMM_allgather(nerrors,1,RPN_COMM_INTEGER,errors,1,RPN_COMM_INTEGER,RPN_COMM_GRID,ierr)
  print *,'TEST INFO:',sum(errors),' unexpected values found from GET',errors
  if(sum(errors) > 0) goto 777
!   if(nerrors .ne. 0) then
!     print *,'TEST INFO:',nerrors,' local unexpected values found from GET'
!   endif

!   if(me == 0) then
    print *,'TEST INFO: GET passive mode test end'
!   endif
! ===================================== remote GET/PUT test ===================================
  target_pe = mod(me+1,npes)       ! remote target
  from_pe = mod(me+npes-1,npes)    ! remote accessor
  print *,'TEST INFO: target PE is',target_pe,' accessing PE is',from_pe
! ===================================== active mode ===================================
  do i = 1 , siz
    my_data(i) = -i - 100*me
  enddo
  nval = 5                         ! number of values for remote put
  offset = 10                      ! displacement into remote window
  offset_r = 20
goto 602
! ==================================
!  if(me == 0) then
    print *,'=============================================================================='
    print *,'TEST INFO: GET/PUT active fence mode test start'
!  endif

  call RPN_COMM_i_win_open(window,.true.,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: accessing/exposing window, PE=',me
  else
    print *,'TEST ERROR: cannot access/expose window, PE=',me
    goto 888
  endif

  local(1:5) = [2,4,6,8,10] + 100*me
  call RPN_COMM_i_win_put_r(window,c_loc(local(1)),target_pe,offset,nint(nval/2.0),ierr)  ! (window,larray,target,offset,nelem,ierr)
  call RPN_COMM_i_win_get_r(window,c_loc(local2(1)),from_pe,offset_r,nint(nval/2.0),ierr)  ! (window,larray,target,offset,nelem,ierr)
  do i = nint(nval/2.0)+1,nval
    call RPN_COMM_i_win_put_r(window,c_loc(local(i)),target_pe,offset+i-1,1,ierr)
    call RPN_COMM_i_win_get_r(window,c_loc(local2(i)),from_pe,offset_r+i-1,1,ierr)
  enddo

  call RPN_COMM_i_win_close(window,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: ending access/exposure to window, PE=',me
  else
    print *,'TEST ERROR: cannot end access/exposure to window, PE=',me
    goto 888
  endif

  nerrors = 0
  local(1:5) = local(1:5) + 100*from_pe - 100*me
  do i = 1 , size(my_data)
    if(i>offset .and. i<=offset+nval) then
      if( my_data(i) .ne. local(i-offset)) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',my_data(i),' expected',local(i-offset),' error'
      endif
    else
      if(my_data(i) .ne. -i - 100*me) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',my_data(i),' expected',-i - 100*me,' error'
      endif
    endif
  enddo
  do i = 1 , nval
      expected = my_data(i+offset_r) + 100*me - 100*from_pe
      got = local2(i)
      if( got .ne. expected ) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',got,' expected',expected,' ERROR',nerrors
      else
        print *,'i=',i,' got',got,' expected',expected,' OK'
      endif
  enddo
  errors = 0
  call RPN_COMM_allgather(nerrors,1,RPN_COMM_INTEGER,errors,1,RPN_COMM_INTEGER,RPN_COMM_GRID,ierr)
  print *,'TEST INFO:',sum(errors),' unexpected values found from GET/PUT',errors
  if(sum(errors) > 0) goto 777
!  if(nerrors .ne. 0) then
!    print *,'TEST INFO:',nerrors),' unexpected values found from PUT'
!  endif

!  if(me == 0) then
    print *,'TEST INFO: GET/PUT active fence mode test end'
!  endif
! ===================================== group mode ===================================
  do i = 1 , siz
    my_data(i) = -i - 100*me
  enddo
  nval = 5                         ! number of values for remote put
  offset = 10                      ! displacement into remote window
  offset_r = 20
! ==================================
!  if(me == 0) then
    print *,'=============================================================================='
    print *,'TEST INFO: GET/PUT active group mode test start'
!  endif
  if(target_pe == from_pe) then
    call RPN_COMM_i_win_group(window,[target_pe],[target_pe],ierr)
  else
    call RPN_COMM_i_win_group(window,[target_pe,from_pe],[target_pe,from_pe],ierr)
  endif
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: group creation OK'
  else
    print *,'TEST ERROR: group creation failure'
    goto 888
  endif
  call RPN_COMM_i_win_open(window,.true.,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: accessing/exposing window, PE=',me
  else
    print *,'TEST ERROR: cannot access/expose window, PE=',me
    goto 888
  endif

  local(1:5) = [2,4,6,8,10] + 100*me
  call RPN_COMM_i_win_put_r(window,c_loc(local(1)),target_pe,offset,nint(nval/2.0),ierr)  ! (window,larray,target,offset,nelem,ierr)
  call RPN_COMM_i_win_get_r(window,c_loc(local2(1)),from_pe,offset_r,nint(nval/2.0),ierr)  ! (window,larray,target,offset,nelem,ierr)
  do i = nint(nval/2.0)+1,nval
    call RPN_COMM_i_win_put_r(window,c_loc(local(i)),target_pe,offset+i-1,1,ierr)
    call RPN_COMM_i_win_get_r(window,c_loc(local2(i)),from_pe,offset_r+i-1,1,ierr)
  enddo

  call RPN_COMM_i_win_close(window,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: ending access/exposure to window, PE=',me
  else
    print *,'TEST ERROR: cannot end access/exposure to window, PE=',me
    goto 888
  endif

  nerrors = 0
  local(1:5) = local(1:5) + 100*from_pe - 100*me
  do i = 1 , size(my_data)
    if(i>offset .and. i<=offset+nval) then
      if( my_data(i) .ne. local(i-offset)) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',my_data(i),' expected',local(i-offset),' error'
      endif
    else
      if(my_data(i) .ne. -i - 100*me) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',my_data(i),' expected',-i - 100*me,' error'
      endif
    endif
  enddo
  do i = 1 , nval
      expected = my_data(i+offset_r) + 100*me - 100*from_pe
      got = local2(i)
      if( got .ne. expected ) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',got,' expected',expected,' ERROR',nerrors
      else
        print *,'i=',i,' got',got,' expected',expected,' OK'
      endif
  enddo
  errors = 0
  call RPN_COMM_allgather(nerrors,1,RPN_COMM_INTEGER,errors,1,RPN_COMM_INTEGER,RPN_COMM_GRID,ierr)
  print *,'TEST INFO:',sum(errors),' unexpected values found from GET/PUT',errors
  if(sum(errors) > 0) goto 777
!  if(nerrors .ne. 0) then
!    print *,'TEST INFO:',nerrors),' unexpected values found from PUT'
!  endif

  call RPN_COMM_i_win_group(window,[-1],[-1],ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: window group deletion OK'
  else
    print *,'TEST ERROR: window group deletion failure'
    goto 888
  endif
!  if(me == 0) then
    print *,'TEST INFO: GET/PUT active group mode test end'
!  endif
! ===================================== passive mode ===================================
602 continue
  do i = 1 , siz
    my_data(i) = -i - 100*me
  enddo
  nval = 5                         ! number of values for remote put
  offset = 10                      ! displacement into remote window
  offset_r = 20
! ==================================
!  if(me == 0) then
    print *,'=============================================================================='
    print *,'TEST INFO: GET/PUT passive mode test start'
!  endif

  call RPN_COMM_i_win_open(window,.false.,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: accessing/exposing window, PE=',me
  else
    print *,'TEST ERROR: cannot access/expose window, PE=',me
    goto 888
  endif

  local(1:5) = [2,4,6,8,10] + 100*me
  nerrors = 0
  call RPN_COMM_i_win_put_r(window,c_loc(local(1)),target_pe,offset,nint(nval/2.0),ierr)  ! (window,larray,target,offset,nelem,ierr)
  if(ierr .ne. RPN_COMM_OK) nerrors = nerrors + 1
  call RPN_COMM_i_win_get_r(window,c_loc(local2(1)),from_pe,offset_r,nint(nval/2.0),ierr)  ! (window,larray,target,offset,nelem,ierr)
  if(ierr .ne. RPN_COMM_OK) nerrors = nerrors + 1
  do i = nint(nval/2.0)+1,nval
    call RPN_COMM_i_win_put_r(window,c_loc(local(i)),target_pe,offset+i-1,1,ierr)
    if(ierr .ne. RPN_COMM_OK) nerrors = nerrors + 1
    call RPN_COMM_i_win_get_r(window,c_loc(local2(i)),from_pe,offset_r+i-1,1,ierr)
    if(ierr .ne. RPN_COMM_OK) nerrors = nerrors + 1
  enddo

  call RPN_COMM_i_win_close(window,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: ending access/exposure to window, PE=',me
  else
    print *,'TEST ERROR: cannot end access/exposure to window, PE=',me
    goto 888
  endif

  local(1:5) = local(1:5) + 100*from_pe - 100*me
  do i = 1 , size(my_data)
    if(i>offset .and. i<=offset+nval) then
      if( my_data(i) .ne. local(i-offset)) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',my_data(i),' expected',local(i-offset),' error'
      endif
    else
      if(my_data(i) .ne. -i - 100*me) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',my_data(i),' expected',-i - 100*me,' error'
      endif
    endif
  enddo
  do i = 1 , nval
      expected = my_data(i+offset_r) + 100*me - 100*from_pe
      got = local2(i)
      if( got .ne. expected ) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',got,' expected',expected,' ERROR',nerrors
      else
        print *,'i=',i,' got',got,' expected',expected,' OK'
      endif
  enddo
  errors = 0
  call RPN_COMM_allgather(nerrors,1,RPN_COMM_INTEGER,errors,1,RPN_COMM_INTEGER,RPN_COMM_GRID,ierr)
  print *,'TEST INFO:',sum(errors),' unexpected values found from GET/PUT',errors
  if(sum(errors) > 0) goto 777
!  if(nerrors .ne. 0) then
!    print *,'TEST INFO:',nerrors),' unexpected values found from PUT'
!  endif

!  if(me == 0) then
    print *,'TEST INFO: GET/PUT passive mode test end'
!  endif
! ===================================== remote ACC test ===================================
  target_pe = mod(me+1,npes)       ! remote target
  from_pe = mod(me+npes-1,npes)    ! remote accessor
  print *,'TEST INFO: target PE is',target_pe,' accessing PE is',from_pe
! ================================== set windows properties ======================================
  call RPN_COMM_i_oper(RPN_COMM_MAX,oper)
  call RPN_COMM_i_win_oper(window,oper,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: associated operator to window, PE=',me
  else
    print *,'TEST ERROR: cannot associate operator to window, PE=',me
    goto 777
  endif
goto 603
! ===================================== active mode ===================================
! initial values for local window are negative, the accumulate operator used is MPI_MAX
! this is then equivalent to a straight PUT
  do i = 1 , siz
    my_data(i) = -i
  enddo
  nval = 5                         ! number of values for remote put
  offset = 10                      ! displacement into remote window
! ==================================
!  if(me == 0) then
    print *,'=============================================================================='
    print *,'TEST INFO: ACC active fence mode test start'
!  endif

  call RPN_COMM_i_win_open(window,.true.,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: accessing/exposing window, PE=',me
  else
    print *,'TEST ERROR: cannot access/expose window, PE=',me
    goto 888
  endif

  local(1:5) = [2,4,6,8,10] + 100*me
  call RPN_COMM_i_win_getacc_r(window,c_loc(local(1)),C_NULL_PTR,.false.,target_pe,offset,nint(nval/2.0),oper,ierr)  ! (window,larray,target,offset,nelem,ierr)
  do i = nint(nval/2.0),nval
    call RPN_COMM_i_win_getacc_r(window,c_loc(local(i)),C_NULL_PTR,.false.,target_pe,offset+i-1,1,oper,ierr)
  enddo

  call RPN_COMM_i_win_close(window,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: ending access/exposure to window, PE=',me
  else
    print *,'TEST ERROR: cannot end access/exposure to window, PE=',me
    goto 888
  endif

  nerrors = 0
  local(1:5) = local(1:5) + 100*from_pe - 100*me
  do i = 1 , size(my_data)
    if(i>offset .and. i<=offset+nval) then
      if( my_data(i) .ne. local(i-offset)) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',my_data(i),' expected',local(i-offset),' error'
      endif
    else
      if(my_data(i) .ne. -i) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',my_data(i),' expected',-i,' error'
      endif
    endif
  enddo
  errors = 0
  call RPN_COMM_allgather(nerrors,1,RPN_COMM_INTEGER,errors,1,RPN_COMM_INTEGER,RPN_COMM_GRID,ierr)
  print *,'TEST INFO:',sum(errors),' unexpected values found from ACC',errors
  if(sum(errors) > 0) goto 777
!  if(nerrors .ne. 0) then
!    print *,'TEST INFO:',nerrors),' unexpected values found from ACC'
!  endif

!  if(me == 0) then
    print *,'TEST INFO: ACC active fence mode test end'
!  endif
! ================================== set windows properties ======================================
603 continue
  print *,'=============================================================================='
!  call RPN_COMM_i_oper(RPN_COMM_BXOR,oper)
  call RPN_COMM_i_oper(RPN_COMM_MAX,oper)
  print *,'TEST INFO: operator=',oper%t2
  call RPN_COMM_i_win_oper(window,oper,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: associated operator to window, PE=',me
  else
    print *,'TEST ERROR: cannot associate operator to window, PE=',me
    goto 777
  endif
! ===================================== passive mode ===================================
! initial values for local window are negative, the accumulate operator used is MPI_MAX
! this is then equivalent to a straight PUT
  do i = 1 , siz
    my_data(i) = -i
  enddo
  nval = 5                         ! number of values for remote put
  offset = 10                      ! displacement into remote window
! ==================================
!  if(me == 0) then
    print *,'TEST INFO: ACC passive mode test start'
!  endif

  call RPN_COMM_i_win_open(window,.false.,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: accessing/exposing window, PE=',me
  else
    print *,'TEST ERROR: cannot access/expose window, PE=',me
    goto 888
  endif

  local(1:5) = [2,4,6,8,10] + 100*me
  nerrors = 0
  call RPN_COMM_i_win_getacc_r(window,c_loc(local(1)),C_NULL_PTR,.false.,target_pe,offset,nint(nval/2.0),oper,ierr)  ! (window,larray,target,offset,nelem,ierr)
  if(ierr .ne. RPN_COMM_OK) nerrors = nerrors + 1
  do i = nint(nval/2.0)+1,nval
    call RPN_COMM_i_win_getacc_r(window,c_loc(local(i)),C_NULL_PTR,.false.,target_pe,offset+i-1,1,oper,ierr)
    if(ierr .ne. RPN_COMM_OK) nerrors = nerrors + 1
  enddo

  call RPN_COMM_i_win_close(window,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: ending access/exposure to window, PE=',me
  else
    print *,'TEST ERROR: cannot end access/exposure to window, PE=',me
    goto 888
  endif

  local(1:5) = local(1:5) + 100*from_pe - 100*me
  do i = 1 , size(my_data)
    got = my_data(i)
    if(i>offset .and. i<=offset+nval) then
      expected = local(i-offset)
!       expected = -i
!      expected = ieor(local(i-offset),expected)
      if( got .ne. expected)  then
        nerrors = nerrors + 1
        print *,'i=',i,' got',got,' expected',expected,'  ERROR',nerrors
      else
        print *,'i=',i,' got',got,' expected',expected,'  OK'
      endif
    else
      expected = -i
      if(got .ne. -i) then
        nerrors = nerrors + 1
        print *,'i=',i,' got',got,' expected',expected,' ERROR',nerrors
      endif
    endif
  enddo
  if(nerrors .ne. 0) print *,'TEST INFO:',nerrors,' unexpected local values found from ACC'
  errors = 0
  call RPN_COMM_allgather(nerrors,1,RPN_COMM_INTEGER,errors,1,RPN_COMM_INTEGER,RPN_COMM_GRID,ierr)
  print *,'TEST INFO:',sum(errors),' unexpected values found from ACC',errors
  if(sum(errors) > 0) goto 777
!  if(nerrors .ne. 0) then
!    print *,'TEST INFO:',nerrors),' unexpected values found from ACC'
!  endif

!  if(me == 0) then
    print *,'TEST INFO: ACC passive mode test end'
!  endif
! ================================== free windows ======================================
777 continue
  call RPN_COMM_i_win_free(window,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: freed window, PE=',me
  else
    print *,'TEST ERROR: cannot free window, PE=',me
    goto 888
  endif

  call RPN_COMM_i_win_free(window2,ierr)
  if(ierr == RPN_COMM_OK) then
    print *,'TEST INFO: freed window2, PE=',me
  else
    print *,'TEST ERROR: cannot free window2, PE=',me
    goto 888
  endif
888 continue
  if(ierr .ne. RPN_COMM_OK .or. sum(errors) > 0) then
    print *,'TEST INFO: ==================== ERRORS IN TEST  ===================='
    if(sum(errors) > 0)print *,'TEST INFO: nerrors=',sum(errors)
    if(ierr .ne. RPN_COMM_OK) print *,'TEST INFO: ierr .ne. RPN_COMM_OK'
  else
    print *,'TEST INFO: ==================== TEST SUCCESSFUL ===================='
  endif
  return
end subroutine RPN_COMM_i_win_test
!===============================================================================
! beginning of USER CALLABLE routines/functions
!===============================================================================
!InTf!
!****f* rpn_comm/RPN_COMM_i_win_group control membership of get and put groups for this PE and this window
! SYNOPSIS
subroutine RPN_COMM_i_win_group(window,pes_to,pes_from,ierr)  !InTf!
!===============================================================================
! for this one sided communication window,
! provide a list of PEs into which remote operations (get/put/acc) will be performed
! provide a list of PEs from which remote operations (get/put/acc) will be accepted
!===============================================================================
!
! window   (IN)      rpn_comm window (see RPN_COMM_types.inc)
! pes_to   (IN)      array containing the list of pe_s who's window will be 
!                    the target of remote operations from this PE (get/put/acc)
! pes_from (IN)      array containing the list of pe_s that will be accessing
!                    via remote operations this PE's window (get/put/acc)
! ierr     (OUT)     will be set to RPN_COMM_OK or RPN_COMM_ERROR
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: rpncomm_window                                       !InTf!
! ARGUMENTS
  type(rpncomm_window), intent(IN) :: window                       !InTf!
  integer, dimension(:), intent(IN) :: pes_to                      !InTf!
  integer, dimension(:), intent(IN) :: pes_from                    !InTf!
  integer, intent(OUT) :: ierr                                     !InTf!
!******

  integer :: ierr2, index
  type(rpncomm_windef), pointer :: win_entry
  integer :: n_to, n_from, i, j, ref
  logical :: duplicates

  ierr = RPN_COMM_ERROR
  if( .not. win_valid(window) )  return  ! bad window reference
  n_to = size(pes_to)
  n_from = size(pes_from)
  if(n_to <= 0 .or. n_from <= 0) return  ! bad list of PEs
  duplicates = .false.
  do i = 1,n_to-1
    ref = pes_to(i)
    do j = i+1 , n_to
      if(pes_to(j) == ref) then
        print *,"RPN_COMM_i_win_group: ERROR, duplicates detected in pes_to"
        duplicates = .true.
        exit
      endif
    enddo
  enddo
  do i = 1,n_from-1
    ref = pes_from(i)
    do j = i+1 , n_from
      if(pes_from(j) == ref) then
        print *,"RPN_COMM_i_win_group: ERROR, duplicates detected in pes_from"
        duplicates = .true.
        exit
      endif
    enddo
  enddo
  if(duplicates) return

!  indx = window%t2                   ! window table entry for this window is win_tab(indx)
  call c_f_pointer( window%p , win_entry )                   ! pointer to win_tab entry (rpncomm_windef)
  if(win_entry%is_open) return           ! no change allowed if window is "exposed"
  if(win_entry%s_group .ne. MPI_GROUP_NULL) then   ! free group, it will be re-creataed
    call MPI_group_free(win_entry%s_group,ierr2)
    if(ierr2 .ne. MPI_SUCCESS) return
  endif
  if(win_entry%r_group .ne. MPI_GROUP_NULL) then   ! free group, it will be re-creataed
    call MPI_group_free(win_entry%r_group,ierr2)
    if(ierr2 .ne. MPI_SUCCESS) return
  endif

  if(pes_to(1)==-1 .and. pes_from(1)==-1) then     ! cancel groups are desired 
    win_entry%s_group = MPI_GROUP_NULL
    win_entry%r_group = MPI_GROUP_NULL
    ierr = RPN_COMM_OK
    return
  endif

  call MPI_group_incl(win_entry%grp, n_to, pes_to, win_entry%s_group, ierr2)
  if(ierr2 .ne. MPI_SUCCESS) return                ! group include failed

  call MPI_group_incl(win_entry%grp, n_from, pes_from, win_entry%r_group, ierr2)
  if(ierr2 .ne. MPI_SUCCESS) return                ! group include failed

  ierr = RPN_COMM_OK
  return
end subroutine RPN_COMM_i_win_group                                   !InTf!
!InTf!
!****f* rpn_comm/RPN_COMM_i_win_oper
! SYNOPSIS
subroutine RPN_COMM_i_win_oper(window,oper,ierr)  !InTf!
!===============================================================================
! for this one sided communication window,
! associate an operator to accumulate operations
!===============================================================================
!
! window   (IN)      rpn_comm window (see RPN_COMM_types.inc)
! oper     (IN)      rpncomm_operator (see RPN_COMM_types.inc)
! ierr     (OUT)     will be set to RPN_COMM_OK or RPN_COMM_ERROR
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
#include <RPN_COMM_interfaces_int.inc>
!!  import :: rpncomm_window                                       !InTf!
!!  import :: rpncomm_operator                                     !InTf!
! ARGUMENTS
  type(rpncomm_window), intent(IN) :: window                       !InTf!
  type(rpncomm_operator), intent(IN) :: oper                       !InTf!
  integer, intent(OUT) :: ierr                                     !InTf!
!******

  integer :: ierr2, index
  type(rpncomm_windef), pointer :: win_entry
  integer :: n_to, n_from

  ierr = RPN_COMM_ERROR
  if(.not. RPN_COMM_i_valid_oper(oper) ) then
    print *,'ERROR: invalid operator in RPN_COMM_i_win_oper'
    return  ! bad operator
  endif
  if( .not. win_valid(window) )  then
    print *,'ERROR: invalid window in RPN_COMM_i_win_oper'
    return  ! bad window reference
  endif

  call c_f_pointer( window%p , win_entry )   ! pointer to win_tab entry (rpncomm_windef)
  win_entry%opr = oper%t2

  ierr = RPN_COMM_OK
  return
end subroutine RPN_COMM_i_win_oper                                   !InTf!
!InTf!
!****f* rpn_comm/RPN_COMM_i_win_post post a "shopping list"
! SYNOPSIS
subroutine RPN_COMM_i_win_post(window,pe,offseti,nelemi,offseto,nelemo,nreq,ierr)          !InTf!
!===============================================================================
! post a set of messages in secondary window pointing to data in primary window
! the secondary window (see RPN_COMM_i_win_create_secondary) of the target PE
! is used to post a pointer to a message in this PE's primary window
!
! first an accumulate is used to bump the target PE's "message list index"
! then the message information (sender/message offset/message length) is
! remotely stored in the target PE's primary window at "message list index"
!
! the secondary window contains groups of WINDEF_MESSAGE_SIZE integers
! the first WINDEF_MESSAGE_SIZE numbers are (last_message_index,....)
!
! offseti(1) = -1 indicates a special reinitialization mode (offseto, nelemo are ignored)
!    nelemi(1) = -1 : reinitialize secondary window, keep slot assignments
!    nelemi(1) = -2 : fully reinitialize secondary window, delete slot assignments
!===============================================================================
!
! window  (IN)      rpn_comm window type from call to RPN_COMM_i_win_create (see RPN_COMM_types.inc)
! pe      (IN)      target pe (array of dimension nreq)
! offseti (IN)      offset in primary window of request message (array of dimension nreq)
! nelemi  (IN)      length of request message (array of dimension nreq)
! offseto (IN)      offset in primary window of reply message (array of dimension nreq)
! nelemo  (IN)      length of reply message (array of dimension nreq)
! nreq    (IN)      number of requests (dimension of arrays pe, offseti, nelemi, offseto, nelemo)
! ierr    (OUT)     RPN_COMM_OK or RPN_COMM_ERROR will be returned
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2016
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: rpncomm_window                                          !InTf!
! ARGUMENTS
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  integer, dimension(nreq), intent(IN) :: pe                          !InTf!
  integer, dimension(nreq), intent(IN) :: offseti                     !InTf!
  integer, dimension(nreq), intent(IN) :: nelemi                      !InTf!
  integer, dimension(nreq), intent(IN) :: offseto                     !InTf!
  integer, dimension(nreq), intent(IN) :: nelemo                      !InTf!
  integer, intent(IN) :: nreq                                         !InTf!
  integer, intent(OUT) :: ierr                                        !InTf!
!******

  integer :: indx, ierr2, incr, i, insert, m
  integer *8 :: offset_8
  type(rpncomm_request_message) :: message
  type(rpncomm_request_message), dimension(:), pointer :: remote2  ! base translated as Fortran pointer to secondary array (messages)
  type(slot), dimension(:), pointer :: slottab
  integer :: pe1, offseti1, nelemi1, offseto1, nelemo1

  ierr = RPN_COMM_ERROR
  if(.not. win_valid(window) ) return

  indx = window%t2

  if(win_tab(indx)%active_mode) THEN  ! OOPS, not passive mode
    print *,'ERROR: rpncomm_window must be open in passive mode to use RPN_COMM_i_win_post'
    return
  ENDIF

  if(C_ASSOCIATED(win_tab(indx)%base2)) then  ! secondary window
    call c_f_pointer(win_tab(indx)%base2,remote2,[win_tab(indx)%nbase2])
  else
    print *,'ERROR: secondary window does not exist for this rpncomm_window'
    return
  endif

  if(C_ASSOCIATED(win_tab(indx)%slots)) then  ! slot table
    call c_f_pointer(win_tab(indx)%slots,slottab,[win_tab(indx)%nslots])
  else
    print *,'ERROR: slot table does not exist for this rpncomm_window'
    return
  endif

  if(offseti(1) == -1 .and. nelemi(1) == -1) then         ! reinitialize secondary window, keep slot assignments
    ierr = RPN_COMM_OK
    remote2(2:) = NULL_rpncomm_request_message        ! clear messages
    return
  endif
  if(offseti(1) == -1 .and. nelemi(1) == -2) then         ! fully reinitialize secondary window
    ierr = RPN_COMM_OK
    remote2(1) = rpncomm_request_message(1,0,0,0,0)   ! counters in slot 1
    remote2(2:) = NULL_rpncomm_request_message        ! clear messages
    return
  endif
  if(offseti(1) == -1) then   ! nelemi(1) not -1 or -2
    print *,'ERROR: invalid secondary window reinitialization mode'
    return
  endif

  do m = 1, nreq   ! loop over requests

    pe1      = pe(m)             ! target PE
    offseti1 = offseti(m)        ! offset in secondary window of input request list
    nelemi1  = nelemi(m)         ! number of items in input request list
    offseto1 = offseto(m)        ! offset in secondary window of output reply list
    nelemo1  = nelemo(m)         ! number of items in output reply list
    if(pe1 >= win_tab(indx)%npe) then
      print *,'ERROR: target pe',pe,' out of range, maxpe=',win_tab(indx)%npe-1
      return
    endif
    incr = 0
    do i = 1, win_tab(indx)%nslots      ! look into slot cache for PE pe1
      if(slottab(i)%pe == pe1) then     ! found, set incr
        incr = slottab(i)%n
        exit
      endif
    enddo

    if(incr == 0) then     ! pe1 not found in slot cache table, get a slot from pe1 (through an atomic index increment)
      call MPI_win_lock(MPI_LOCK_EXCLUSIVE, pe1, 0, win_tab(indx)%win2, ierr2)   ! lock target for get/accumulate
      incr = 1
      offset_8 = 0
      call MPI_accumulate(incr,1,MPI_INTEGER,pe1,offset_8,1,MPI_INTEGER,MPI_SUM,win_tab(indx)%win2,ierr2) ! add 1 to index
      call MPI_get(incr,1,MPI_INTEGER,pe1,offset_8,1,MPI_INTEGER,win_tab(indx)%win2, ierr2)               ! get new index
      call MPI_win_unlock(pe1, win_tab(indx)%win2, ierr2)   ! make sure incr is usable locally
      if(incr > win_tab(indx)%nbase2) then
        print *,"ERROR: no space left in target pe's (",pe,") message slots"
        return
      endif
      insert = 0
      do i = 1, win_tab(indx)%nslots  ! insert into slot cache table
        if(slottab(i)%pe == -1) then
          slottab(i)%n = incr
          slottab(i)%pe = pe1
          insert = i
          exit
        endif
      enddo
      if(insert == 0) then
        print *,"ERROR: slot table full"
        return
      endif
    endif

!     call mpi_win_lock(MPI_LOCK_EXCLUSIVE, pe1, 0, win_tab(indx)%win2, ierr2)   ! exclusive lock for put, no overlap, safer
    call mpi_win_lock(MPI_LOCK_SHARED, pe1, 0, win_tab(indx)%win2, ierr2)        ! shared lock for put, no overlap, faster
    message = rpncomm_request_message( win_tab(indx)%rank, offseti1, nelemi1, offseto1, nelemo1)
    offset_8 = WINDEF_MESSAGE_SIZE*(incr - 1)    ! WINDEF_MESSAGE_SIZE integers per entry , first WINDEF_MESSAGE_SIZE integers=[index,.....]
    call MPI_put(message,WINDEF_MESSAGE_SIZE,MPI_INTEGER,pe1,offset_8,WINDEF_MESSAGE_SIZE,MPI_INTEGER,win_tab(indx)%win2, ierr2)  ! send "shopping list" to target pe
    call mpi_win_unlock(pe1, win_tab(indx)%win2, ierr2)

  enddo

  ierr = RPN_COMM_OK
  return
end subroutine RPN_COMM_i_win_post                                    !InTf!
!InTf!
!****f* rpn_comm/RPN_COMM_i_win_create_secondary create a one sided secondary communication window
! SYNOPSIS
subroutine RPN_COMM_i_win_create_secondary(window,slots,ierr)  !InTf!
!===============================================================================
! add one sided communication secondary window to existing window (user exposed interface)
!===============================================================================
!
! window (IN)      rpn_comm window type from call to RPN_COMM_i_win_create (see RPN_COMM_types.inc)
! slots  (IN)      number of message slots (max message targets)
!                  slots < 0 means reallocate window with size abs(nslots)
! ierr   (OUT)     RPN_COMM_OK or RPN_COMM_ERROR will be returned
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2016
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window, rpncomm_datatype, rpncomm_communicator  !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  integer, intent(IN) :: slots                                        !InTf!
!******
  integer :: indx, extent, siz, nslots
  integer(kind=MPI_ADDRESS_KIND) win_size  ! may be wider than a default integer
  type(rpncomm_request_message), dimension(:), pointer :: remote2  ! base translated as Fortran pointer to secondary array (messages)
  type(slot), dimension(:), pointer :: slottab
  integer :: ierr2
  integer(C_SIZE_T) :: msize

  nslots = abs(slots)
  ierr = RPN_COMM_ERROR
  if(.not. win_valid(window) ) return

  indx = window%t2
  if(C_ASSOCIATED(win_tab(indx)%base2)) then  ! secondary window exists
    if(slots > 0) then
      print *,'ERROR: secondary window already exists'
      return
    else
      call MPI_win_free(win_tab(indx)%win2,ierr2)      ! free secondary window
      if(ierr2 .ne. MPI_SUCCESS) then
        print *,'ERROR: error freeing secondary window'
      endif
      call c_f_pointer(win_tab(indx)%base2,remote2,[win_tab(indx)%nbase2])
      call MPI_free_mem(remote2,ierr2)                 ! free secondary window memory
      win_tab(indx)%win2 = MPI_WIN_NULL                ! MPI null window
!       print *,'DEBUG: deallocating slot table, size=',win_tab(indx)%nslots
      call c_f_pointer(win_tab(indx)%slots,slottab,[win_tab(indx)%nslots])
!       print *,slottab
      call c_free(win_tab(indx)%slots)
      win_tab(indx)%base2 = C_NULL_PTR
      win_tab(indx)%nbase2 = 0
!       deallocate(slottab)                              ! free slot table
      win_tab(indx)%slots = C_NULL_PTR
      win_tab(indx)%nslots = 0
    endif
  endif

  call MPI_TYPE_EXTENT(MPI_INTEGER, extent, ierr2)  ! determine size associated with MPI integer
  siz = nslots + 1
  win_size = WINDEF_MESSAGE_SIZE * extent * siz
  call MPI_alloc_mem(win_size, MPI_INFO_NULL, win_tab(indx)%base2,ierr2)  ! allocate window memory
  if(ierr2 .ne. MPI_SUCCESS) then
    print *,'ERROR: cannot allocate memory for secondary window'
    return
  endif

  win_tab(indx)%nbase2 = siz
  call c_f_pointer(win_tab(indx)%base2,remote2,[siz])    ! make fortran pointer
  remote2(1)  = rpncomm_request_message(1,0,0,0,0)                   ! slot 1 is used as counter
  remote2(2:) = NULL_rpncomm_request_message                         ! initialize
  call MPI_win_create(remote2, win_size, extent, MPI_INFO_NULL, win_tab(indx)%com, win_tab(indx)%win2,ierr2)  ! create secondary window
  if(ierr2 .ne. MPI_SUCCESS) then
    print *,'ERROR: creation of secondary window unsuccessful'
    return
  endif
!   print *,'INFO: allocated slot table, size=',nslots
!   allocate(slottab(nslots))   ! allocate slot table
  
!   win_tab(indx)%slots = C_LOC(slottab(1))
  msize = nslots*4*WINDEF_SLOT_SIZE
  win_tab(indx)%slots = c_malloc(msize)
  win_tab(indx)%nslots = nslots
  call c_f_pointer(win_tab(indx)%slots,slottab,[win_tab(indx)%nslots])
  slottab(1:nslots) = slot(-1,-1)

!   if(slots < 0) print *,'INFO: reallocated secondary window'
  ierr = RPN_COMM_OK
  return
end subroutine RPN_COMM_i_win_create_secondary                                  !InTf!
!InTf!
!****f* rpn_comm/RPN_COMM_i_win_create create a one sided primary communication window
! SYNOPSIS
subroutine RPN_COMM_i_win_create(window,dtype,siz,com,array,ierr)  !InTf!
!===============================================================================
! create a one sided communication primary window (user exposed interface)
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
!******
  integer :: indx

  ierr = RPN_COMM_ERROR
  window = NULL_rpncomm_window

  call create_win_entry(array,dtype%t2,siz,com%t2,indx,ierr)
  if(ierr .ne. RPN_COMM_OK) return
!print *,'DEBUG: window created indx =',indx
!print *,'DEBUG: window created win =',win_tab(indx)%win
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

  integer :: indx, ierr1, ierr2
  integer, dimension(:), pointer :: slottab   ! slot table
  integer, dimension(:), pointer :: remote   ! base translated as Fortran pointer to primary array
  type(rpncomm_request_message), dimension(:), pointer :: remote2  ! base translated as Fortran pointer to secondary array (messages)

  ierr = RPN_COMM_ERROR
  if(.not. win_valid(window) ) return
  indx = window%t2

!print *,'DEBUG: window is valid, index =',indx
!print *,'DEBUG: window is valid, win   =',win_tab(indx)%win
  call MPI_win_free(win_tab(indx)%win,ierr2)       ! free primary window
  if(ierr2 .ne. MPI_SUCCESS) then
    print *,'ERROR: error freeing primary window'
  endif
  if(debug_mode) print *,"DEBUG: freed primary window"

  if(win_tab(indx)%win2 .ne. MPI_WIN_NULL) then
    call MPI_win_free(win_tab(indx)%win2,ierr2)      ! free secondary window if it exists
    if(ierr2 .ne. MPI_SUCCESS) then
      print *,'ERROR: error freeing secondary window'
    endif
    if(debug_mode) print *,"DEBUG: freed secondary window"
  endif

  if(.not. win_tab(indx)%is_user) then   ! internal storage allocation for primary window, release it
    call c_f_pointer(win_tab(indx)%base,remote,[win_tab(indx)%nbase])
    call MPI_free_mem(remote,ierr1)
    if(ierr1 .ne. MPI_SUCCESS) then
      print *,'ERROR: error freeing primary window storage'
    endif
    win_tab(indx)%base = C_NULL_PTR
    if(debug_mode) print *,"DEBUG: freed primary window storage"
  else
    ierr1 = MPI_SUCCESS
  endif

  if(C_ASSOCIATED(win_tab(indx)%base2)) then   ! internal storage allocation for secondary window, release it if it exists
    call c_f_pointer(win_tab(indx)%base2,remote2,[win_tab(indx)%nbase2])
    call MPI_free_mem(remote2,ierr1)
    if(ierr1 .ne. MPI_SUCCESS) then
      print *,'ERROR: error freeing secondary window storage'
    endif
    win_tab(indx)%base2 = C_NULL_PTR
    if(debug_mode) print *,"DEBUG: freed secondary window storage"
  endif

  if(C_ASSOCIATED(win_tab(indx)%slots)) then   ! internal storage allocation for message slots, release it if it exists
    call c_free(win_tab(indx)%slots)
!     call c_f_pointer(win_tab(indx)%slots,slottab,[win_tab(indx)%nslots])
!     deallocate(slottab)
    win_tab(indx)%slots = C_NULL_PTR
    if(debug_mode) print *,"DEBUG: freed slot table"
  endif
  win_tab(indx) = NULL_rpncomm_windef              ! blank entry
  win_tab(indx)%win  = MPI_WIN_NULL                ! MPI null window
  win_tab(indx)%win2 = MPI_WIN_NULL                ! MPI null window
  win_tab(indx)%com = MPI_COMM_NULL                ! MPI null communicator
  win_tab(indx)%typ = MPI_DATATYPE_NULL            ! MPI null datatype

  window = NULL_rpncomm_window
  if(ierr1 == MPI_SUCCESS .and. ierr2 == MPI_SUCCESS) then
    ierr = RPN_COMM_OK
  endif

  return
end subroutine RPN_COMM_i_win_free                                    !InTf!

!InTf!
!****f* rpn_comm/RPN_COMM_i_win_open open a one sided communication window
! SYNOPSIS
subroutine RPN_COMM_i_win_open(window,active,ierr)                           !InTf!
!===============================================================================
! "expose" a one sided communication window (see RPN_COMM_i_win_create)
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! active(IN)      .true. : use active mode, .false.: use passive mode
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
  logical, intent(IN) :: active                                       !InTf!
!******

  integer :: ierr1, ierr2, indx
  logical ::is_open
  logical, external :: RPN_COMM_i_win_check

  ierr = RPN_COMM_ERROR
  indx = window%t2                   ! window table entry for this window
  is_open = RPN_COMM_i_win_check(window,ierr2)
  if(is_open .or. ierr2 .ne. RPN_COMM_OK) return    ! ERROR: window is already open (exposed) or not valid

  win_tab(indx)%active_mode = active
  if(win_tab(indx)%active_mode) then           ! active mode
    if(win_tab(indx)%s_group == MPI_GROUP_NULL .and. win_tab(indx)%r_group == MPI_GROUP_NULL) then  ! fence mode
      ierr1 = MPI_SUCCESS
!      call MPI_win_fence(MPI_MODE_NOPRECEDE,win_tab(indx)%win,ierr2)  ! fence does not complete any sequence of locally issued RMA calls
      call MPI_win_fence(0,win_tab(indx)%win,ierr2)
      if(debug_mode) print *,'DEBUG: opening window, fence mode',win_tab(indx)%win,win_tab(indx)%grp
!      call MPI_win_start(win_tab(indx)%grp,0,win_tab(indx)%win,ierr1)    ! start RMA access epoch on win to members of s_group
!      call MPI_win_post (win_tab(indx)%grp,0,win_tab(indx)%win,ierr2)    ! start RMA exposure epoch on win from members of r_group
    else                                                                                         ! start/complete/post/wait mode
      call MPI_win_start(win_tab(indx)%s_group,0,win_tab(indx)%win,ierr1)    ! start RMA access epoch on win to members of s_group
      call MPI_win_post (win_tab(indx)%r_group,0,win_tab(indx)%win,ierr2)    ! start RMA exposure epoch on win from members of r_group
      if(debug_mode) print *,'DEBUG: opening window, groups mode',win_tab(indx)%win,win_tab(indx)%s_group,win_tab(indx)%r_group
    endif
  else                                         ! barrier call in passive mode
    if(debug_mode) print *,'DEBUG: opening window, passive mode',win_tab(indx)%win
    call MPI_barrier(win_tab(indx)%com,ierr1)
    ierr2 = MPI_SUCCESS
  endif

  if(ierr1 .eq. MPI_SUCCESS .and. ierr2 .eq. MPI_SUCCESS) then
    ierr = RPN_COMM_OK
    win_tab(indx)%is_open = .true.   ! set window open flag
  endif
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

  integer :: ierr1, ierr2, indx
  logical :: is_not_open
  logical, external :: RPN_COMM_i_win_check

  ierr = RPN_COMM_ERROR
  indx = window%t2                   ! window table entry for this window
  is_not_open = .not. RPN_COMM_i_win_check(window,ierr2)
  if(is_not_open .or. ierr2 .ne. RPN_COMM_OK) return    ! ERROR: window is not open (exposed) or not valid

  if(win_tab(indx)%active_mode) then           ! active mode
    if( win_tab(indx)%s_group == MPI_GROUP_NULL .and. win_tab(indx)%r_group == MPI_GROUP_NULL ) then  ! fence mode
      ierr1 = MPI_SUCCESS
!      call MPI_win_fence(MPI_MODE_NOSUCCEED,win_tab(indx)%win,ierr2)   ! fence does not start any sequence of locally issued RMA calls
      call MPI_win_fence(0,win_tab(indx)%win,ierr2)
      if(debug_mode) print *,'DEBUG: fence mode closing window',win_tab(indx)%win
!      call MPI_win_complete(win_tab(indx)%win, ierr1)    ! Complete RMA access epoch on win started MPI_WIN_START
!      call MPI_win_wait    (win_tab(indx)%win, ierr2)    ! Complete RMA exposure epoch started by MPI_WIN_POST on win
    else                                                                                            ! start/complete/post/wait mode
      call MPI_win_complete(win_tab(indx)%win, ierr1)    ! Complete RMA access epoch on win started MPI_WIN_START
      call MPI_win_wait    (win_tab(indx)%win, ierr2)    ! Complete RMA exposure epoch started by MPI_WIN_POST on win
      if(debug_mode) print *,'DEBUG: group mode closing window',win_tab(indx)%win
    endif
  else                                         ! nothing to do in passive mode
    if(debug_mode) print *,'DEBUG: passive mode closing window',win_tab(indx)%win
    ierr1 = MPI_SUCCESS
    call MPI_barrier(win_tab(indx)%com,ierr2)   ! barrier call in passive mode
  endif

  if(ierr1 .eq. MPI_SUCCESS .and. ierr2 .eq. MPI_SUCCESS) then
    ierr = RPN_COMM_OK
    win_tab(indx)%is_open = .false.   ! unset window open flag
  endif
  return

end subroutine RPN_COMM_i_win_close                                   !InTf!

!InTf!
!****f* rpn_comm/RPN_COMM_i_win_valid check if a one sided communication window is valid
! SYNOPSIS
function RPN_COMM_i_valid_win(window,ierr) result(is_valid)           !InTf!
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

end function RPN_COMM_i_valid_win                                     !InTf!

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
function RPN_COMM_i_win_get_ptr(window,np,ierr) result(ptr)              !InTf!
!===============================================================================
! get a one sided communication window (see RPN_COMM_i_win_create) data pointer
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! np     (IN)     1 : get pointer to primary window, 2: get pointer to secondary window
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
  integer, intent(IN)  :: np                                          !InTf!
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  type(C_PTR) :: ptr                                                  !InTf!
!******

  integer :: indx

  ierr = RPN_COMM_ERROR
  ptr = C_NULL_PTR
  if(.not. win_valid(window) ) return   ! invalid window
  if(np < 1 .or. np > 2) return         ! bad np code (must be 1 or 2)

  indx = window%t2            ! window table entry for this window
  if(np == 1) ptr = win_tab(indx)%base    ! get primary data pointer from window table entry
  if(np == 2) ptr = win_tab(indx)%base2   ! get secondary data pointer from window table entry

  ierr = RPN_COMM_OK
  return

end function RPN_COMM_i_win_get_ptr                                      !InTf!

!InTf!
!****f* rpn_comm/RPN_COMM_i_win_get_size get number of elements in primary window
! SYNOPSIS
function RPN_COMM_i_win_get_size(window,ierr) result(siz)                 !InTf!
!===============================================================================
! get size of a one sided communication window (see RPN_COMM_i_win_create) 
! (nb of elements in primary window)
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! ierr   (OUT)    error status, RPN_COMM_OK or RPN_COMM_ERROR
!
! function value : number of elements in the primary window
!                  in case of error, -1 is returned
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: rpncomm_window                                          !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  integer :: siz                                                  !InTf!
!******

  integer :: indx

  ierr = RPN_COMM_ERROR
  siz = -1
  if(.not. win_valid(window) ) return

  indx = window%t2            ! window table entry for this window
  siz = win_tab(indx)%siz     ! get size from window table entry

  ierr = RPN_COMM_OK
  return

end function RPN_COMM_i_win_get_size                                      !InTf!

!InTf!
!****f* rpn_comm/RPN_COMM_i_win_get_size get size of secondary window
! SYNOPSIS
function RPN_COMM_i_win_get_size2(window,ierr) result(siz)                 !InTf!
!===============================================================================
! get size of secondary window (see RPN_COMM_i_win_create_secondary)
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! ierr   (OUT)    error status, RPN_COMM_OK or RPN_COMM_ERROR
!
! function value : size of secondary window (request/reply) in "words"
!                  in case of error, -1 is returned
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: rpncomm_window                                          !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  integer :: siz                                                  !InTf!
!******

  integer :: indx

  ierr = RPN_COMM_ERROR
  siz = -1
  if(.not. win_valid(window) ) return

  indx = window%t2            ! window table entry for this window
  siz = win_tab(indx)%nbase2  ! get size from window table entry

  ierr = RPN_COMM_OK
  return

end function RPN_COMM_i_win_get_size2                                      !InTf!

!InTf!
!****f* rpn_comm/RPN_COMM_i_win_put_r write into a remote one sided communication window
! SYNOPSIS
subroutine RPN_COMM_i_win_put_r(window,larray,targetpe,offset,nelem,ierr) !InTf!
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
  integer, intent(IN) :: targetpe                                     !InTf!
  integer, intent(IN) :: offset                                       !InTf!
  integer, intent(IN) :: nelem                                        !InTf!
!******

  logical :: is_open
  integer :: ierr2, indx
  logical, external :: RPN_COMM_i_win_check
  integer, dimension(:), pointer :: local
  type(rpncomm_windef), pointer :: win_entry
  integer(kind=MPI_ADDRESS_KIND) :: offset_8

  ierr = RPN_COMM_ERROR
  is_open = RPN_COMM_i_win_check(window,ierr2)
  if( (.not. is_open) .or. (ierr2 .ne. MPI_SUCCESS) )  return  ! bad window reference or window not open (exposed)

  indx = window%t2
  call c_f_pointer( window%p , win_entry )                   ! pointer to win_tab entry (rpncomm_windef)
  if(offset+nelem > win_entry%siz) return                ! out of bounds condition for "remote" array
  call c_f_pointer( larray , local, [nelem* win_entry%ext] )   ! pointer to local array

  if(.not. win_entry%active_mode) then
    call mpi_win_lock(MPI_LOCK_EXCLUSIVE, targetpe, 0, win_entry%win, ierr2)
    if(ierr2 .ne. MPI_SUCCESS) then
      print *,'ERROR: failed to lock window',win_entry%win,' on PE',targetpe
      return
    else
      if(debug_mode) print *,'INFO: locked window',win_entry%win,' on PE',targetpe
    endif
  endif
  offset_8 = offset
  if(debug_mode) print *,'DEBUG: PUT to',targetpe,' at',offset_8+1,' :',local(1:nelem)
  call MPI_put(local,       nelem,        win_entry%typ,   targetpe,   offset_8,    nelem,        win_entry%typ,   win_entry%win,ierr2)
!              ORIGIN_ADDR, ORIGIN_COUNT, ORIGIN_DATATYPE, TARGET_RANK,TARGET_DISP, TARGET_COUNT, TARGET_DATATYPE, WIN,          IERROR)
  if(ierr2 .ne. MPI_SUCCESS) then
    print *,'ERROR: remote PUT failed in window',win_entry%win
    return
  endif
  if(.not. win_entry%active_mode) then
    call mpi_win_unlock(targetpe, win_entry%win, ierr2)
    if(ierr2 .ne. MPI_SUCCESS) then
      print *,'ERROR: failed to unlock window',win_entry%win,' on PE',targetpe
      return
    else
      if(debug_mode) print *,'INFO: unlocked window',win_entry%win,' on PE',targetpe
    endif
  endif

  ierr = RPN_COMM_OK
  return

end subroutine RPN_COMM_i_win_put_r                                   !InTf!

!InTf!
!****f* rpn_comm/RPN_COMM_i_win_getacc_r get/accumulate into a remote one sided communication window
! SYNOPSIS
subroutine RPN_COMM_i_win_getacc_r(window,larray,garray,before,targetpe,offset,nelem,oper,ierr) !InTf!
!===============================================================================
! one sided communication remote get/accumulate from/into a one sided communication window
! from a local array, into another local array
! it is an error to attempt a "remote" get/accumulate when the window is not "exposed"
! larray must not be a NULL pointer
! if garray is not a NULL pointer, a get and accumulate operation is performed
! if garray is a NULL pointer, a simple accumulate will be performed
!
! window (IN)     rpn_comm one sided window type(rpncomm_window) (see RPN_COMM_types.inc)
! larray (IN)     C compatible pointer (type(C_PTR)) to local array (source of accumulate)
! garray (IN)     C compatible pointer (type(C_PTR)) to local array (destination of get)
! before (IN)     if .true. get before accumulate, otherwise get after accumulate
! target (IN)     ordinal in window communicator of remote PE
! offset (IN)     displacement (origin 0) into remote PE window data array
! nelem  (IN)     number of elements to transfer (type of element was defined at window creation)
! oper   (IN)     rpn comm operator (see RPN_COMM_types.inc and RPN_COMM_i_oper)
! ierr   (OUT)    error status, RPN_COMM_OK or RPN_COMM_ERROR
!===============================================================================
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use RPN_COMM_windows
  implicit none
!!  import :: C_PTR                                                   !InTf!
!!  import :: rpncomm_window                                          !InTf!
!!  import :: rpncomm_operator                                        !InTf!
! ARGUMENTS
  integer, intent(OUT) :: ierr                                        !InTf!
  type(rpncomm_window), intent(IN) :: window                          !InTf!
  type(rpncomm_operator), intent(IN) :: oper                          !InTf!
  type(C_PTR), intent(IN), value :: larray                            !InTf!
  type(C_PTR), intent(IN), value :: garray                            !InTf!
  logical, intent(IN) :: before                                       !InTf!
  integer, intent(IN) :: targetpe                                     !InTf!
  integer, intent(IN) :: offset                                       !InTf!
  integer, intent(IN) :: nelem                                        !InTf!
!******

  logical :: is_open, after
  integer :: ierr2, indx
  logical, external :: RPN_COMM_i_win_check
  integer, dimension(:), pointer :: local, local2
  type(rpncomm_windef), pointer :: win_entry
  integer(kind=MPI_ADDRESS_KIND) :: offset_8

  after = .not. before
  ierr = RPN_COMM_ERROR
  is_open = RPN_COMM_i_win_check(window,ierr2)
  if( (.not. is_open) .or. (ierr2 .ne. MPI_SUCCESS) )  return  ! bad window reference or window not open (exposed)

  indx = window%t2
  call c_f_pointer( window%p , win_entry )                   ! pointer to win_tab entry (rpncomm_windef)
  if(offset+nelem > win_entry%siz) return                ! out of bounds condition for "remote" array
  call c_f_pointer( larray , local, [nelem* win_entry%ext] )   ! pointer to local array
  if(C_ASSOCIATED(garray)) call c_f_pointer( garray , local2, [nelem* win_entry%ext] )

  if(.not. win_entry%active_mode) then
    call mpi_win_lock(MPI_LOCK_EXCLUSIVE, targetpe, 0, win_entry%win, ierr2)
    if(ierr2 .ne. MPI_SUCCESS) then
      print *,'ERROR: failed to lock window',win_entry%win,' on PE',targetpe
      return
    else
      if(debug_mode) print *,'INFO: locked window',win_entry%win,' on PE',targetpe
    endif
  endif
  offset_8 = offset
  if(debug_mode) print *,'DEBUG: ACC to',targetpe,' at',offset_8+1,oper%t2,' :',local(1:nelem)
  if(C_ASSOCIATED(garray) .and. before) then
    call MPI_get(local2,nelem,win_entry%typ,targetpe,offset_8,nelem,win_entry%typ,win_entry%win,ierr2)
    if(ierr2 .ne. MPI_SUCCESS) then
      print *,'ERROR: GET before ACC failed in window',win_entry%win
      return
    endif
  endif
  if(win_entry%opr == MPI_REPLACE) then                  ! replace operation, use put instead of accumulate
    call MPI_put       (local,       nelem,        win_entry%typ,   targetpe,   offset_8,    nelem,        win_entry%typ,           win_entry%win, ierr2)
  else
    call MPI_accumulate(local,       nelem,        win_entry%typ,   targetpe,   offset_8,    nelem,        win_entry%typ, oper%t2,  win_entry%win, ierr2)
!              ORIGIN_ADDR, ORIGIN_COUNT, ORIGIN_DATATYPE, TARGET_RANK,TARGET_DISP, TARGET_COUNT, TARGET_DATATYPE, WIN,          IERROR)
  endif
  if(ierr2 .ne. MPI_SUCCESS) then
    print *,'ERROR: remote ACC failed in window',win_entry%win
    return
  endif
  if(C_ASSOCIATED(garray) .and. after) then
    call MPI_get(local2,nelem,win_entry%typ,targetpe,offset_8,nelem,win_entry%typ,win_entry%win,ierr2)
    if(ierr2 .ne. MPI_SUCCESS) then
      print *,'ERROR: GET after ACC failed in window',win_entry%win
      return
    endif
  endif
  if(.not. win_entry%active_mode) then
    call mpi_win_unlock(targetpe, win_entry%win, ierr2)
    if(ierr2 .ne. MPI_SUCCESS) then
      print *,'ERROR: failed to unlock window',win_entry%win,' on PE',targetpe
      return
    else
      if(debug_mode) print *,'INFO: unlocked window',win_entry%win,' on PE',targetpe
    endif
  endif

  ierr = RPN_COMM_OK
  return

end subroutine RPN_COMM_i_win_getacc_r                                   !InTf!

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
  integer, dimension(:), pointer :: remote   ! base translated as Fortran pointer to primary array

  ierr = RPN_COMM_ERROR
  is_open = RPN_COMM_i_win_check(window,ierr2)
  if( (is_open) .or. (ierr2 .ne. MPI_SUCCESS) )  return  ! bad window reference or window open (exposed)

  indx = window%t2
  call c_f_pointer( window%p , win_entry )                   ! pointer to win_tab entry (rpncomm_windef)
  if(offset+nelem > win_entry%siz) return                ! out of bounds condition for "remote" array
  extent = win_entry%ext
  call c_f_pointer( larray , local, [nelem*extent] )                   ! pointer to local array
  call c_f_pointer(win_entry%base,remote,[win_entry%nbase])
  do i = 1, nelem*extent
    remote(i+offset*extent) = local(i)
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
  integer(kind=MPI_ADDRESS_KIND) offset_8

  ierr = RPN_COMM_ERROR
  is_open = RPN_COMM_i_win_check(window,ierr2)
  if( (.not. is_open) .or. (ierr2 .ne. MPI_SUCCESS) )  return  ! bad window reference or window not open (exposed)

  indx = window%t2
  call c_f_pointer( window%p , win_entry )                   ! pointer to win_tab entry (rpncomm_windef)
  if(offset+nelem > win_entry%siz) return                ! out of bounds condition for "remote" array
  call c_f_pointer( larray , local, [nelem*win_entry%ext] )   ! pointer to local array

  if(.not. win_entry%active_mode) then
    call mpi_win_lock(MPI_LOCK_SHARED, target, 0, win_entry%win, ierr2) ! read mode, lock can be shared
    if(ierr2 .ne. MPI_SUCCESS) return
  endif

  offset_8 = offset
  call MPI_get(local,nelem,win_entry%typ,target,offset_8,nelem,win_entry%typ,win_entry%win,ierr2) 
  if(ierr2 .ne. MPI_SUCCESS) return

  if(.not. win_entry%active_mode) then
    call mpi_win_unlock(target, win_entry%win, ierr2)
    if(ierr2 .ne. MPI_SUCCESS) return
  endif

  ierr = RPN_COMM_OK
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
  integer, dimension(:), pointer :: remote   ! base translated as Fortran pointer to primary array

  ierr = RPN_COMM_ERROR
  is_open = RPN_COMM_i_win_check(window,ierr2)
  if( (.not. is_open) .or. (ierr2 .ne. MPI_SUCCESS) )  return  ! bad window reference or window open (exposed)

  indx = window%t2
  call c_f_pointer( window%p , win_entry )                   ! pointer to win_tab entry (rpncomm_windef)
  if(offset+nelem > win_entry%siz) return                ! out of bounds condition for "remote" array
  extent = win_entry%ext
  call c_f_pointer( larray , local, [nelem*extent] )                   ! pointer to local array
  call c_f_pointer(win_entry%base,remote,[win_entry%nbase])
  do i = 1, nelem*extent
     local(i) = remote(i+offset*extent)
  enddo

  ierr = RPN_COMM_OK
  return

end subroutine RPN_COMM_i_win_get_l                                   !InTf!
