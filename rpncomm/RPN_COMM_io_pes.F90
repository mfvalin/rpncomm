!/! RMNLIB - Library of useful routines for C and FORTRAN programming
! ! Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
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
#if defined(SELF_TEST)
#define STAND_ALONE
#endif
!
! STAND_ALONE mode is mostly used for debugging
! in normal mode, these routines rely on internal module rpn_comm
!
module RPN_COMM_io_pe_tables   ! this module may also be used by data distribution routines
  use ISO_C_BINDING
#if ! defined(STAND_ALONE)
  use rpn_comm
#endif
  type :: RPN_COMM_io_set
    integer, dimension(:), pointer :: x    ! x coordinate in grid of IO PEs in this set
    integer, dimension(:), pointer :: y    ! y coordinate in grid of IO PEs in this set
    integer :: ioset                       ! ID for this IO PE set
    integer :: comm                        ! MPI communicator  for this IO PE set     if it belongs to the set (MPI_COMM_NULL otherwise)
    integer :: me                          ! ordinal of this PE in above communicator if it belongs to the set (-1 otherwise)
    integer :: megrid                      ! ordinal of this PE in grid               if it belongs to the set (-1 otherwise)
    integer :: npe                         ! number of OP PEs in this set
    integer :: ngroups                     ! number of PE groups in this set (size of a group may not exceed MIN(pe_nx,pe_ny) )
    integer :: groupsize                   ! size of groups (last group may be shorter)
    integer :: reserved                    ! unused for now
  end type
  integer, save :: iosets=0
  integer, save :: scramblex = 0
  integer, save :: scrambley = 0
  integer, PARAMETER :: MAXIOSETS=16
  type(RPN_COMM_io_set), dimension(MAXIOSETS), save :: io_set
  integer, dimension(16) :: primes = [ &
         5,      7,      11,     13,   &
        17,     19,      23,     29,   &
        31,     37,      41,     43,   &
        47,     53,      59,     61  ]
contains
!
! helper routine used by create_ioset
! for now only method 0 is supported for IO PE dispersion
!
  subroutine make_io_pe_list(x,y,npes,pe_nx,pe_ny,setno,method)
    integer, dimension(npes), intent(OUT) :: x   ! x coordinates of PEs in set
    integer, dimension(npes), intent(OUT) :: y   ! y coordinates of PEs in set
    integer, intent(IN) :: npes                  ! number of PEs in set
    integer, intent(IN) :: pe_nx, pe_ny          ! number of PEs along x and y in PE grid
    integer, intent(IN) :: setno, method         ! set number, fill method
    integer :: i
    integer :: deltax, deltay, pe_nxy
!
    if(scramblex == 0) then    ! initialize
      scramblex = 1
      scrambley = 1
      if(pe_nx > pe_ny) then   ! scramblex = lowest number that is prime with respect to pe_nx
        do i = 1 , size(primes)
          scramblex = primes(i)
          if(mod(pe_nx,primes(i)) .ne. 0) exit
        enddo
      else                     ! scrambley = lowest number that is prime with respect to pe_ny
        do i = 1 , size(primes)
          scrambley = primes(i)
          if(mod(pe_ny,primes(i)) .ne. 0) exit
        enddo
      endif
!      print *,"DEBUG: scramblex, scrambley =",scramblex, scrambley
    endif
    x = -1
    y = -1
    if(method .ne. 0) return      ! method 0 is the only one supported for the time being
    deltax = 1
    deltay = 1
    pe_nxy = min(pe_nx,pe_ny)
    if(method == 0) then
      deltax = scramblex
      deltay = scrambley
    endif
    if( npes > pe_nxy * pe_nxy ) return
    safe = 0
    do i = 0 , npes-1
      if(npes > pe_nxy) then
        if(pe_nx > pe_ny) then
          x(i+1) = mod( i * deltax , pe_nx )
          y(i+1) = mod( mod( i , pe_ny) + i / pe_nx , pe_ny)
        else
          x(i+1) = mod( mod( i , pe_nx ) + i / pe_ny , pe_nx)
          y(i+1) = mod( i * deltay , pe_ny)
        endif
      else
        x(i+1) = mod( i * deltax , pe_nx )
        y(i+1) = mod( i * deltay , pe_ny )
      endif
    enddo
  end subroutine make_io_pe_list
!
  function is_io_pe(setno) result(ordinal)  ! is this pe part of IO PE set setno. if yes return rank in set, else return -1
    implicit none
    integer, intent(IN) :: setno             ! set number as returned by create_ioset
    integer :: ordinal
    integer :: i
#if defined(STAND_ALONE)
    integer :: pe_mex, pe_mey, pe_me, ierr
    integer, external :: RPN_COMM_mype
    ierr = RPN_COMM_mype(pe_me,pe_mex,pe_mey)  ! who am I, where am I ?
#endif
!
    ordinal = -1
    if(setno < 1 .or. setno > iosets) return    ! setno out of bounds
    if(io_set(setno)%ioset .ne. setno) return   ! set no longer valid
!
    do i = 1 , io_set(setno)%npe
      if(pe_mex == io_set(setno)%x(i) .and. pe_mey == io_set(setno)%y(i)) then  ! I am an IO pe
        ordinal =  io_set(setno)%me
        return
      endif
    enddo
    
  end function is_io_pe
!
  function io_pe_comm(setno) result(communicator)  ! is this pe part of IO PE set setno. if yes return set communicator, else return MPI_COMM_NULL
    implicit none
    include 'mpif.h'
    integer, intent(IN) :: setno             ! set number as returned by create_ioset
    integer :: communicator
    integer :: i
#if defined(STAND_ALONE)
    integer :: pe_mex, pe_mey, pe_me, ierr
    integer, external :: RPN_COMM_mype
    ierr = RPN_COMM_mype(pe_me,pe_mex,pe_mey)  ! who am I, where am I ?
#endif

    communicator = MPI_COMM_NULL
    if(setno < 1 .or. setno > iosets) return    ! setno out of bounds
    if(io_set(setno)%ioset .ne. setno) return   ! set no longer valid
!
    do i = 1 , io_set(setno)%npe
      if(pe_mex == io_set(setno)%x(i) .and. pe_mey == io_set(setno)%y(i)) then  ! I am an IO pe
        communicator =  io_set(setno)%comm
        return
      endif
    enddo
    
  end function io_pe_comm
!
  function io_pe_list(setno) result(list)  ! return pointer to list of PEs if set is valid, else return a null pointer
    implicit none
    integer, intent(IN) :: setno             ! set number as returned by create_ioset
    integer, dimension(:,:), pointer :: list
    integer :: npes

    nullify(list)
    if(setno < 1 .or. setno > iosets) return    ! setno out of bounds
    if(io_set(setno)%ioset .ne. setno) return   ! set no longer valid
!
    npes = io_set(setno)%npe
    allocate(list(npes,2))
    list(1:npes,1) = io_set(setno)%x            ! list(:,1) x coordinates of IO PEs
    list(1:npes,2) = io_set(setno)%y            ! list(:,2) y coordinates of IO PEs
  end function io_pe_list
!
  function io_pe_call(setno,callback,argv) result(status)
    implicit none
    integer, intent(IN) :: setno      ! set number as returned by create_ioset
    type(C_PTR) :: argv               ! pointer describing arguments to callback
    integer, external :: callback     ! user routine to be called
    integer :: status
    integer, dimension(:), pointer :: argvf

    status = -1
    if(setno < 1 .or. setno > iosets) return    ! setno out of bounds
    if(io_set(setno)%ioset .ne. setno) return   ! set no longer valid
!
    status = 0
    if(io_set(setno)%me < 0) return
!
    status = callback(argv,setno,io_set(setno)%comm,io_set(setno)%me,   &
                      io_set(setno)%x,io_set(setno)%y,io_set(setno)%npe)
  end function io_pe_call
!
  function io_pe_size(setno) result(population)  ! return set population if set is valid, else return -1
    implicit none
    integer, intent(IN) :: setno             ! set number as returned by create_ioset
    integer :: population

    population = -1
    if(setno < 1 .or. setno > iosets) return    ! setno out of bounds
    if(io_set(setno)%ioset .ne. setno) return   ! set no longer valid
!
    population = io_set(setno)%npe
  end function io_pe_size
!
  function io_pe_groups(setno) result(ngroups)  ! return number of groups if set is valid, else return -1
    implicit none
    integer, intent(IN) :: setno             ! set number as returned by create_ioset
    integer :: ngroups

    ngroups = -1
    if(setno < 1 .or. setno > iosets) return    ! setno out of bounds
    if(io_set(setno)%ioset .ne. setno) return   ! set no longer valid
!
    ngroups = io_set(setno)%ngroups
  end function io_pe_groups
!
  function free_ioset(setno) result(freed)   ! free table space associated with IO PE set setno
    implicit none
    include 'mpif.h'
    integer, intent(IN) :: setno             ! set number as returned by create_ioset
    integer :: freed
    integer :: ierr
!
    freed = -1
    if(setno < 1 .or. setno > iosets) return    ! setno out of bounds
    if(io_set(setno)%ioset .ne. setno) return   ! set no longer valid
!
    deallocate(io_set(setno)%x)
    nullify(io_set(setno)%x)
    deallocate(io_set(setno)%y)
    nullify(io_set(setno)%y)
    io_set(setno)%ioset = -1
    io_set(setno)%npe = -1
    if(io_set(setno)%comm .ne. MPI_COMM_NULL) call mpi_comm_free(io_set(setno)%comm,ierr)
    io_set(setno)%comm = MPI_COMM_NULL
    io_set(setno)%me = -1
    io_set(setno)%ngroups = 0
    io_set(setno)%groupsize = 0
    freed = setno
  end function free_ioset
!
! check that this PE set is valid (no duplicates) and print IO PE map
! also check that no column has 2 IO PEs in a group and neither has any row
! used to validate what create_ioset has produced
!
  function check_ioset(newset, x ,y, npes, pe_nx, pe_ny, pe_me, diag) result(status)
    implicit none
    integer, intent(IN) :: newset, npes         ! set number, number of PEs in set
    integer, intent(IN) :: pe_nx, pe_ny         ! number of PEs along x and y
    integer :: pe_me                            ! ordinal in grid of this PE
    logical, intent(IN) :: diag                 ! if .true., PE 0 prints the IO PE map for this set
    integer, intent(IN), dimension(npes) :: x   ! x coordinates of IO PEs
    integer, intent(IN), dimension(npes) :: y   ! y coordinates of IO PEs
    integer :: status                           ! 0 if set is OK, -1 otherwise
!
    integer*1, dimension(0:pe_nx-1,0:pe_ny-1) :: safe
    integer*1, dimension(0:pe_ny-1) :: row      ! there are pe_ny rows
    integer*1, dimension(0:pe_nx-1) :: col      ! there are pe_nx columns
    integer :: i, j, groupsize, low, high
!
    status = -1
    safe = 0
    groupsize = min(pe_nx, pe_ny)
    do i = 1, npes, groupsize  ! loop over groups
      row = 0
      col = 0
      do j = i , min(i+groupsize-1,npes)      ! for PE at coordinates ( x(j), y(j) )
        row(y(j)) = row(y(j)) + 1             ! add one to count for this row
        col(x(j)) = col(x(j)) + 1             ! add one to count for this column
      enddo
      if(any(row > 1) .or. any(col > 1) ) then
        print *,"ERROR: problem creating IO set, there are 2 or more PEs on row or column"
        print *,"ERROR: x = ",x
        print *,"ERROR: row = ",row
        print *,"ERROR: y = ",y
        print *,"ERROR: col = ",col
        status = -1
        return
      endif
    enddo
    do i = 1 , npes
      if(safe(x(i),y(i)) .ne. 0) then  ! OOPS, 2 PEs in this set are the same
        print *,"ERROR: problem creating IO set, there are duplicates"
        status = -1
        return
      else
        safe(x(i),y(i)) = 1 + mod(  ( (i-1) / min(pe_nx, pe_ny) ) , 9)  ! group number, 9 if group number > 9
      endif
    enddo
    if(pe_me == 0 .and. diag)then
      print 101,"===== IO PE map for set",newset," (",(npes+min(pe_nx, pe_ny)-1)/min(pe_nx, pe_ny)," groups) ====="
      do j = pe_ny-1 , 0 , -1
        print 100,j,safe(:,j)
      enddo
    endif
100   format(1X,I4,1x,128I1)
101   format(A,I3,A,I3,A)
    status = 0
    return
  end function check_ioset
!
  function create_ioset(npes,method) result(setno)  ! pe_indomm AKA "GRID" is assumed as a communicator
    implicit none
    integer, intent(IN) :: npes
    integer, intent(IN) :: method
    integer :: setno
    integer :: ierr, my_color, i, j, newset

#if defined(STAND_ALONE)
    include 'mpif.h'
    integer :: pe_mex, pe_mey, pe_indomm, pe_me, pe_nx, pe_ny
    integer, external :: RPN_COMM_comm, RPN_COMM_mype
    pe_indomm = RPN_COMM_comm('GRID')          ! grid communicator
    ierr = RPN_COMM_mype(pe_me,pe_mex,pe_mey)  ! who am I, where am I ?
    call mpi_allreduce(pe_mex,pe_nx,1,MPI_INTEGER,MPI_MAX,pe_indomm,ierr)  ! roundabout way to get pe_nx
    pe_nx = pe_nx + 1
    call mpi_allreduce(pe_mey,pe_ny,1,MPI_INTEGER,MPI_MAX,pe_indomm,ierr)  ! roundabout way to get pe_ny
    pe_ny = pe_ny + 1
#endif
!
    setno = -1                       ! precondition for failure
    if(iosets >= MAXIOSETS) return   ! OOPS, table is full
    newset = iosets + 1
    do i = 1 , iosets
      if(io_set(i)%ioset == -1) newset = i   ! recycle freed set
    enddo
    iosets = max(iosets,newset)
!
    allocate(io_set(newset)%x(npes))
    allocate(io_set(newset)%y(npes))
    io_set(newset)%ioset = -1
    io_set(newset)%npe = -1
    io_set(newset)%comm = MPI_COMM_NULL
    io_set(newset)%me = -1
    io_set(newset)%megrid = -1
    io_set(newset)%ngroups = 0
    io_set(newset)%groupsize = 0
    call make_io_pe_list(io_set(newset)%x, io_set(newset)%y, npes, pe_nx, pe_ny, newset, method)
    if(io_set(newset)%x(1) .ne. -1) then
      if( check_ioset(newset,io_set(newset)%x ,io_set(newset)%y, npes, pe_nx, pe_ny, pe_me, .true.) == -1 ) io_set(newset)%x(1) = 1
    endif
    if(io_set(newset)%x(1) == -1) then        ! miserable failure, list of IO pes could not be created
      deallocate(io_set(newset)%x,io_set(newset)%y)
      nullify(io_set(newset)%x)
      nullify(io_set(newset)%y)
!      setno = -1    ! not needed, it is already -1
      return
    endif
    io_set(newset)%ioset = newset
    io_set(newset)%npe = npes
    io_set(newset)%groupsize = min(pe_nx, pe_ny, npes)   ! size of groups for this IO PE set (last group may be shorter)
    io_set(newset)%ngroups = (npes+io_set(newset)%groupsize-1)/(io_set(newset)%groupsize)   ! number of groups in this IO PE set
    my_color = MPI_UNDEFINED
    do i = 1 , npes
      if(pe_mex == io_set(newset)%x(i) .and. pe_mey == io_set(newset)%y(i)) then  ! I am an IO pe
        my_color = 1
        exit
      endif
    enddo
    call MPI_COMM_SPLIT(pe_indomm,my_color,pe_me,io_set(newset)%comm,ierr)   ! communicator for this set
    if(io_set(newset)%comm .ne. MPI_COMM_NULL) then                          ! get rank in communicator if part of set
      call mpi_comm_rank(io_set(newset)%comm,io_set(newset)%me,ierr)
      io_set(newset)%megrid = pe_me                                          ! rank in grid
    endif
    setno = newset
#if defined(SELF_TEST)
    if(pe_me == 0) then
      print *,'DEBUG: creating IO set ',newset,' npes=',npes
      print *,'DEBUG: npe, groups, me',io_set(newset)%npe,io_set(newset)%ngroups,io_set(newset)%me
    endif
#endif
  end function create_ioset
!
end module RPN_COMM_io_pe_tables
!
!=========================   embedded self test   ===============================
!
#if defined(STAND_ALONE)
  subroutine RPN_COMM_io_pe_test(pe_nx,pe_ny,pe_me,pe_mex,pe_mey)
  use ISO_C_BINDING
  implicit none
  integer, external :: RPN_COMM_create_io_set, RPN_COMM_free_io_set, RPN_COMM_io_pe_call
  integer, external :: RPN_COMM_is_io_pe, RPN_COMM_io_pe_size
  interface
    function RPN_COMM_io_pe_list(setno) result(list)
      implicit none
      integer, intent(IN) :: setno
      integer, dimension(:,:), pointer :: list
    end function RPN_COMM_io_pe_list
  end interface
  integer, intent(IN) :: pe_nx,pe_ny,pe_me,pe_mex,pe_mey
#else
  subroutine RPN_COMM_io_pe_test()
  use rpn_comm
  implicit none
#include <RPN_COMM_interfaces_int.inc>
#endif
  integer, external :: RPN_COMM_io_pe_test_callback
  integer setno,nio,me_io,setno2,setno3,status
  integer, dimension(1), target :: argv
  integer, dimension(:,:), pointer :: iopelist
  integer, dimension(1) :: tbcst
!
  if(pe_me == 0)  print *,'DEBUG: pe_nx,pe_ny',pe_nx,pe_ny
  nio = min(pe_nx,pe_ny)
  print 100,'RPN_COMM_io_pe test program, pe_nx,pe_ny,pe_me,pe_mex,pe_mey,nio=',pe_nx,pe_ny,pe_me,pe_mex,pe_mey,nio
  setno = RPN_COMM_create_io_set(nio+2,0)
  me_io = RPN_COMM_is_io_pe(setno)
  if(me_io .ne. -1) then
    print *,"I am a proud IO pe !"
  else
    print *,"I am a lazy  NON-IO pe !"
  endif
  print *,"set number, size of set='",setno,RPN_COMM_io_pe_size(setno)
  setno2 = RPN_COMM_create_io_set(nio,0)
  print *,"set number, size of set='",setno2,RPN_COMM_io_pe_size(setno2)
  setno = RPN_COMM_free_io_set(setno)
  print *,'DEBUG: freed IO set ',setno
  print *,"set number, size of set='",setno,RPN_COMM_io_pe_size(setno)
  setno = RPN_COMM_create_io_set(nio,0)
  print *,"set number, size of set='",setno,RPN_COMM_io_pe_size(setno)
  setno3 = RPN_COMM_create_io_set(nio-1,0)
  print *,"set number, size of set='",setno3,RPN_COMM_io_pe_size(setno3)
  print *,"set number, size of set='",4,RPN_COMM_io_pe_size(setno3)
  argv(1) = pe_me
  status = RPN_COMM_io_pe_call(setno3,RPN_COMM_io_pe_test_callback,C_LOC(argv(1)))
  print *,"after callback, status,argv=",status,argv(1)
  iopelist => RPN_COMM_io_pe_list(setno3)
  print *,"PE list x=",iopelist(:,1)
  print *,"PE list y=",iopelist(:,2)
  tbcst = pe_me
  if(RPN_COMM_is_io_pe(setno3) .ne. -1)then  ! part of the set only, bcst pe_me of highest rank in set
    call RPN_COMM_io_pe_bcast(tbcst,1,'MPI_INTEGER',RPN_COMM_io_pe_size(setno3)-1,setno3,status) 
    print *,'tbcst after broadcast',tbcst
  endif

100 format(A,10I5)
  return
end subroutine
!    status = callback(argv,setno,io_set(setno)%comm,io_set(setno)%me,   &
!                      io_set(setno)%x,io_set(setno)%y,io_set(setno)%npe)
function RPN_COMM_io_pe_test_callback(argv,setno,comm,me,x,y,npe) result(status)
  use ISO_C_BINDING
  implicit none
  integer, intent(IN) :: setno,comm,me,npe
  integer, dimension(npe), intent(IN) :: x,y
  type(C_PTR), intent(IN) :: argv
  integer :: status
  integer, dimension(:), pointer :: argvf
!
  call C_F_POINTER(argv,argvf,[1])
  print *,"IO PE CALLBACK set, me, argument = ",setno,me,argvf(1)
  print *,"CALLBACK x=",x
  print *,"CALLBACK y=",y
  argvf(1) = -(100 +argvf(1))
  status = -123
  return
end
!
!=========================   START of user callable functions   ===============================
!
function RPN_COMM_create_io_set(npes,method) result(setno)  !InTf!  ! create a set of IO PEs
  use RPN_COMM_io_pe_tables
  implicit none                                             !InTf!
  integer, intent(IN)  :: npes      !InTf!  ! number of IO PEs desired in set ( set size will be min(npes,pe_nx,pe_ny) )
  integer, intent(IN)  :: method    !InTf!  ! method 0 is the only one supported for the time being
  integer :: setno                  !InTf!  ! set number created ( -1 if creation failed )

  setno = create_ioset(npes,method)
  return
end function RPN_COMM_create_io_set  !InTf!  
!
function RPN_COMM_free_io_set(setno) result(freed)  !InTf!  ! delete a set of IO PEs created by RPN_COMM_create_io_set
  use RPN_COMM_io_pe_tables
  implicit none                       !InTf!
  integer, intent(IN) :: setno        !InTf!  ! set number as returned by RPN_COMM_create_io_set
  integer :: freed                    !InTf!  ! setno if operation succeeded, -1 if it failed
!
  freed = free_ioset(setno)
end function RPN_COMM_free_io_set     !InTf!  
!
function RPN_COMM_is_io_pe(setno) result(ordinal)   !InTf!  ! is this PE part of IO set setno ? 
  use RPN_COMM_io_pe_tables
  implicit none                      !InTf!
  integer, intent(IN) :: setno       !InTf!  ! set number as returned by RPN_COMM_create_io_set
  integer :: ordinal                 !InTf!  ! ordinal in IO PE set if a member, -1 if not
  ordinal = is_io_pe(setno)
end function RPN_COMM_is_io_pe        !InTf!  
!
function RPN_COMM_io_pe_comm(setno) result(communicator)  !InTf!   ! get communicator of IO set setno
  use RPN_COMM_io_pe_tables
  implicit none                       !InTf!
  integer, intent(IN) :: setno        !InTf!  ! set number as returned by RPN_COMM_create_io_set
  integer :: communicator             !InTf!  ! communicator for IO PE set if a member, null communicator otherwise
  communicator = io_pe_comm(setno)
end function RPN_COMM_io_pe_comm  !InTf!  
!
function RPN_COMM_io_pe_size(setno) result(population)   !InTf!  ! get size of IO set setno
  use RPN_COMM_io_pe_tables
  implicit none                       !InTf!
  integer, intent(IN) :: setno        !InTf!  ! set number as returned by RPN_COMM_create_io_set
  integer :: population               !InTf!  ! population of IO PE set, -1 if set is not valid
  population = io_pe_size(setno)
end function RPN_COMM_io_pe_size      !InTf!  
!
function RPN_COMM_io_pe_groups(setno) result(ngroups)   !InTf!  ! get number of groups in IO set setno
  use RPN_COMM_io_pe_tables
  implicit none                       !InTf!
  integer, intent(IN) :: setno        !InTf!  ! set number as returned by RPN_COMM_create_io_set
  integer :: ngroups                  !InTf!  ! ngroups of IO PE set, -1 if set is not valid
  ngroups = io_pe_groups(setno)
end function RPN_COMM_io_pe_groups    !InTf!  
!
function RPN_COMM_io_pe_list(setno) result(list)   !InTf!  ! get list of PEs in IO set setno
  use RPN_COMM_io_pe_tables
  implicit none                               !InTf!
  integer, intent(IN) :: setno                !InTf!  ! set number as returned by RPN_COMM_create_io_set
  integer, dimension(:,:), pointer :: list    !InTf!  ! list of IO PE set, null pointer if set is not valid
  list => io_pe_list(setno)                    ! list(:,1) x coordinates, list(:,2) y coordinates of IO PEs
end function RPN_COMM_io_pe_list              !InTf!  
!
! callback to be executed only on members of an IO PE set (usually called by all PEs on grid)
! if this PE is not a member of the IO set, nothing happens, status = 0
! if the IO set is invalid, nothing happens, status = -1
! if this PE is a member of the set, callback is called and its return value is put into status
!
! the calling sequence of callback is as follows
! status = callback(argv,setno,comm,ordinal,x,y,npe)
! type(C_PTR) :: argv 
! integer :: setno         ! IO PE set number
! integer :: comm          ! communicator used by this IO PE set
! integer :: ordinal       ! ordinal of this PE in the set
! integer, dimension(npe) :: x  ! x grid coordinates for IO PEs in this set
! integer, dimension(npe) :: y  ! y grid coordinates for IO PEs in this set
! integer :: npe           ! number of PEs in this set
!
function RPN_COMM_io_pe_call(setno,callback,argv) result(status) !InTf!
  use RPN_COMM_io_pe_tables
!!  import :: C_PTR      !InTf!
  implicit none                       !InTf!
  integer, intent(IN) :: setno        !InTf!  ! set number as returned by create_ioset
  type(C_PTR) :: argv                 !InTf!  ! "pickled" arguments to callback
  integer, external :: callback       !InTf!  ! user routine to be called
  integer :: status                   !InTf!
!
  status = io_pe_call(setno,callback,argv)
!
end function RPN_COMM_io_pe_call  !InTf!  
!
! broadcast within an IO PE set
!
subroutine RPN_COMM_io_pe_bcast(buffer,count,datatype,root,setno,ierr) ! cannot publish interface because of buffer
  use RPN_COMM_io_pe_tables
  implicit none
  integer, intent(IN) :: setno                ! set number as returned by RPN_COMM_create_io_set
  integer, intent(IN) :: count                ! number of elements of type datatype
  integer, intent(INOUT), dimension(*) :: buffer   ! data to be broadcast
  character (len=*), intent(IN) :: datatype   ! datatype RPN_COMM style
  integer, intent(IN) :: root                 ! root PE in set for broadcast (typically 0)
  integer, intent(OUT) :: ierr                ! status upon completion
  integer :: dtype
  integer, external :: RPN_COMM_datyp
!
  call MPI_bcast(buffer,count,RPN_COMM_datyp(datatype),root,io_pe_comm(setno),ierr)
end subroutine RPN_COMM_io_pe_bcast
!
!=========================    END of user callable functions    ===============================
!
#if defined SELF_TEST
program self_test
  implicit none
  include 'mpif.h'
  integer ierr,pe_me,pe_mex,pe_mey,npes,pe_nx,pe_ny,pe_indomm
  integer, external :: RPN_COMM_mype, RPN_COMM_comm

  call mpi_init(ierr)
  pe_nx = 1
  pe_ny = 1
  call mpi_comm_size(MPI_COMM_WORLD,npes,ierr)
  pe_indomm = RPN_COMM_comm('GRID')
  ierr = RPN_COMM_mype(pe_me,pe_mex,pe_mey)
  ierr = RPN_COMM_mype(pe_me,pe_mex,pe_mey)
  call mpi_allreduce(pe_mex,pe_nx,1,MPI_INTEGER,MPI_MAX,pe_indomm,ierr)
  pe_nx = pe_nx + 1
  call mpi_allreduce(pe_mey,pe_ny,1,MPI_INTEGER,MPI_MAX,pe_indomm,ierr)
  pe_ny = pe_ny + 1
  call RPN_COMM_io_pe_test(pe_nx,pe_ny,pe_me,pe_mex,pe_mey)
  call mpi_finalize(ierr)
  stop
  end
!
integer function RPN_COMM_comm(what)
  implicit none
  include 'mpif.h'
  character (len=*) :: what
  RPN_COMM_comm = MPI_COMM_WORLD
  return
end
integer function RPN_COMM_datyp(what)
  implicit none
  include 'mpif.h'
  character (len=*) :: what
  RPN_COMM_datyp = MPI_UNDEFINED
  if(trim(what) .ne. 'MPI_INTEGER') return
  RPN_COMM_datyp = MPI_INTEGER
end function RPN_COMM_datyp
integer function RPN_COMM_mype(pe_me,pe_mex,pe_mey)
  implicit none
  include 'mpif.h'
  integer, intent(OUT) :: pe_me,pe_mex,pe_mey
  integer :: i
  integer :: ierr,npes

  call mpi_comm_rank(MPI_COMM_WORLD,pe_me,ierr)
  call mpi_comm_size(MPI_COMM_WORLD,npes,ierr)
  RPN_COMM_mype = 0
  pe_mex = npes - 1
  pe_mey = 0
  do i=7,2,-1
    if(mod(npes,i) == 0 .and. npes .ne. i)then
      pe_mey = (pe_me) / i
      pe_mex = mod(pe_me, i)
      exit
    endif
  enddo
end function RPN_COMM_mype
#endif
