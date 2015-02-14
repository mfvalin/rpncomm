#if defined(SELF_TEST)
#define STAND_ALONE
#endif
!
! STAND_ALONE mode is mostly used for debugging
! in normal mode, these routines rely on internal module rpn_comm
!
module RPN_COMM_io_pe_tables   ! this module may also be used by data distribution routines
#if ! defined(STAND_ALONE)
  use rpn_comm
#endif
  type :: RPN_COMM_io_set
    integer, dimension(:), pointer :: x    ! x coordinate in grid of IO PEs in this set
    integer, dimension(:), pointer :: y    ! y coordinate in grid of IO PEs in this set
    integer :: ioset                       ! ID for this IO PE set
    integer :: comm                        ! MPI communicator  for this IO PE set
    integer :: me                          ! ordinal of this PE in above communicator
    integer :: npe                         ! number of OP PEs in this set
    integer :: ngroups                     ! number of PE groups in this set (size of a group may not exceed MIN(pe_nx,pe_ny) )
    integer :: groupsize                   ! size of groups (last group may be shorter)
  end type
  integer, save :: iosets=0
  integer, save :: scramblex = 0
  integer, save :: scrambley = 0
  integer, PARAMETER :: MAXIOSETS=16
  type(RPN_COMM_io_set), dimension(MAXIOSETS), save :: io_set
contains
!
  subroutine make_io_pe_list(x,y,npes,pe_nx,pe_ny,setno,method)
    integer, dimension(npes), intent(OUT) :: x
    integer, dimension(npes), intent(OUT) :: y
    integer, intent(IN) :: npes
    integer, intent(IN) :: pe_nx, pe_ny
    integer, intent(IN) :: setno, method
    real :: deltax, deltay
    integer :: i
!
    x = -1
    y = -1
    if(method .ne. 0) return      ! method 0 is the only one supported for the time being
    if(npes > min(pe_nx,pe_ny)) return
    deltax = (pe_nx -1.0) / (npes - 1)
    deltay = (pe_ny -1.0) / (npes - 1)
    do i = 1 , npes
      x(i) = nint( (i-1) * deltax )
!      if(mod(method,4)>=2) x(i) = mod(x(i)+setno-1,pe_nx)  ! spread along x 
      y(i) = nint( (i-1) * deltay )
!      if(mod(method,2)>=1) y(i) = mod(y(i)+setno-1,pe_ny)  ! spread along y
    enddo
  end subroutine make_io_pe_list
!
  function is_io_pe(setno) result(ordinal)  ! is this pe part of IO PE set setno. if yes return rank in set, else return -1
    implicit none
    integer, intent(IN) :: setno
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
    integer, intent(IN) :: setno
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
    integer, intent(IN) :: setno
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
  function io_pe_size(setno) result(population)  ! return set population if set is valid, else return -1
    implicit none
    integer, intent(IN) :: setno
    integer :: population

    population = -1
    if(setno < 1 .or. setno > iosets) return    ! setno out of bounds
    if(io_set(setno)%ioset .ne. setno) return   ! set no longer valid
!
    population = io_set(setno)%npe
  end function io_pe_size
!
  function free_ioset(setno) result(freed)
    implicit none
    include 'mpif.h'
    integer, intent(IN) :: setno
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
  function create_ioset(npes,method) result(setno)
    implicit none
    integer, intent(IN) :: npes
    integer, intent(IN) :: method
    integer :: setno
    integer :: ierr, my_color, i, newset

#if defined(STAND_ALONE)
    include 'mpif.h'
    integer :: pe_mex, pe_mey, pe_medomm, pe_me, pe_nx, pe_ny
    integer, external :: RPN_COMM_comm, RPN_COMM_mype
    pe_medomm = RPN_COMM_comm('GRID')          ! grid communicator
    ierr = RPN_COMM_mype(pe_me,pe_mex,pe_mey)  ! who am I, where am I ?
    call mpi_allreduce(pe_mex,pe_nx,1,MPI_INTEGER,MPI_MAX,pe_medomm,ierr)  ! roundabout way to get pe_nx
    pe_nx = pe_nx + 1
    call mpi_allreduce(pe_mey,pe_ny,1,MPI_INTEGER,MPI_MAX,pe_medomm,ierr)  ! roundabout way to get pe_ny
    pe_ny = pe_ny + 1
#endif
!
    setno = -1
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
    io_set(newset)%ngroups = 0
    io_set(setno)%groupsize = 0
    call make_io_pe_list(io_set(newset)%x,io_set(newset)%y,npes,pe_nx,pe_ny,newset,method)
    if(io_set(newset)%x(1) == -1) then        ! miserable failure, list of IO pes could not be created
      deallocate(io_set(newset)%x,io_set(newset)%y)
      nullify(io_set(newset)%x)
      nullify(io_set(newset)%y)
      setno = -1
      return
    endif
    io_set(newset)%ioset = newset
    io_set(newset)%npe = npes
    io_set(setno)%groupsize = min(pe_nx,pe_ny)   ! size of groups for this IO PE set (last group may be shorter)
    io_set(newset)%ngroups = (npes+io_set(setno)%groupsize-1)/(io_set(setno)%groupsize)   ! number of groups in this IO PE set
    my_color = MPI_UNDEFINED
    do i = 1 , npes
      if(pe_mex == io_set(newset)%x(i) .and. pe_mey == io_set(newset)%y(i)) then  ! I am an IO pe
        my_color = 1
        exit
      endif
    enddo
    call MPI_COMM_SPLIT(pe_medomm,my_color,pe_me,io_set(newset)%comm,ierr)   ! communicator for this set
    if(io_set(newset)%comm .ne. MPI_COMM_NULL) then                          ! get rank in communicator if part of set
      call mpi_comm_rank(io_set(newset)%comm,io_set(newset)%me,ierr)
    endif
    setno = newset
#if defined(SELF_TEST)
    if(pe_me == 0) then
      print *,'DEBUG: creating IO set ',newset,' npes=',npes
      print *,'DEBUG: npe, me',io_set(newset)%npe,io_set(newset)%me
      print *,'DEBUG: x=',io_set(newset)%x
      print *,'DEBUG: y=',io_set(newset)%y
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
  implicit none
  integer, intent(IN) :: pe_nx,pe_ny,pe_me,pe_mex,pe_mey
#else
  subroutine RPN_COMM_io_pe_test()
  use rpn_comm
  implicit none
#endif
  integer, external :: RPN_COMM_is_io_pe, RPN_COMM_create_io_set, RPN_COMM_free_io_set, RPN_COMM_io_pe_size
  integer setno,nio,me_io,setno2,setno3
!
  if(pe_me == 0)  print *,'DEBUG: pe_nx,pe_ny',pe_nx,pe_ny
  nio = min(pe_nx,pe_ny)
  print 100,'RPN_COMM_io_pe test program, pe_nx,pe_ny,pe_me,pe_mex,pe_mey,nio=',pe_nx,pe_ny,pe_me,pe_mex,pe_mey,nio
  setno = RPN_COMM_create_io_set(nio,0)
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
  setno3 = RPN_COMM_create_io_set(nio,0)
  print *,"set number, size of set='",setno3,RPN_COMM_io_pe_size(setno3)
  print *,"set number, size of set='",4,RPN_COMM_io_pe_size(setno3)
  

100 format(A,10I5)
  return
end subroutine
!
!=========================   START of user callable functions   ===============================
!
function RPN_COMM_create_io_set(npes,method) result(setno)  ! create a set of IO PEs
  use RPN_COMM_io_pe_tables
  implicit none
  integer, intent(IN)  :: npes    ! number of IO PEs desired in set ( set size will be min(npes,pe_nx,pe_ny) )
  integer, intent(IN)  :: method  ! method 0 is the only one supported for the time being
  integer :: setno                ! set number created ( -1 if creation failed )

  setno = create_ioset(npes,method)
  return
end function RPN_COMM_create_io_set
!
function RPN_COMM_free_io_set(setno) result(freed)  ! delete a set of IO PEs created by RPN_COMM_create_io_set
  use RPN_COMM_io_pe_tables
  implicit none
  integer, intent(IN) :: setno        ! set number as returned by RPN_COMM_create_io_set
  integer :: freed                    ! setno if operation succeeded, -1 if it failed
!
  freed = free_ioset(setno)
end function RPN_COMM_free_io_set
!
function RPN_COMM_is_io_pe(setno) result(ordinal)   ! is this PE part of IO set setno ? 
  use RPN_COMM_io_pe_tables
  implicit none
  integer, intent(IN) :: setno        ! set number as returned by RPN_COMM_create_io_set
  integer :: ordinal                  ! ordinal in IO PE set if a member, -1 if not
  ordinal = is_io_pe(setno)
end function RPN_COMM_is_io_pe
!
function RPN_COMM_io_pe_comm(setno) result(communicator)   ! get communicator of IO set setno
  use RPN_COMM_io_pe_tables
  implicit none
  integer, intent(IN) :: setno        ! set number as returned by RPN_COMM_create_io_set
  integer :: communicator             ! communicator for IO PE set if a member, null communicator otherwise
  communicator = io_pe_comm(setno)
end function RPN_COMM_io_pe_comm
!
function RPN_COMM_io_pe_size(setno) result(population)   ! get size of IO set setno
  use RPN_COMM_io_pe_tables
  implicit none
  integer, intent(IN) :: setno        ! set number as returned by RPN_COMM_create_io_set
  integer :: population               ! population of IO PE set, -1 if set is not valid
  population = io_pe_size(setno)
end function RPN_COMM_io_pe_size
!
function RPN_COMM_io_pe_list(setno) result(list)   ! get size of IO set setno
  use RPN_COMM_io_pe_tables
  implicit none
  integer, intent(IN) :: setno                ! set number as returned by RPN_COMM_create_io_set
  integer, dimension(:,:), pointer :: list    ! list of IO PE set, null pointer if set is not valid
  list = io_pe_list(setno)                    ! list(:,1) x coordinates, list(:,2) y coordinates of IO PEs
end function RPN_COMM_io_pe_list
!
!=========================    END of user callable functions    ===============================
!
#if defined SELF_TEST
program self_test
  implicit none
  include 'mpif.h'
  integer ierr,pe_me,pe_mex,pe_mey,npes,pe_nx,pe_ny,pe_medomm
  integer, external :: RPN_COMM_mype, RPN_COMM_comm

  call mpi_init(ierr)
  pe_nx = 1
  pe_ny = 1
  call mpi_comm_size(MPI_COMM_WORLD,npes,ierr)
  pe_medomm = RPN_COMM_comm('GRID')
  ierr = RPN_COMM_mype(pe_me,pe_mex,pe_mey)
  ierr = RPN_COMM_mype(pe_me,pe_mex,pe_mey)
  call mpi_allreduce(pe_mex,pe_nx,1,MPI_INTEGER,MPI_MAX,pe_medomm,ierr)
  pe_nx = pe_nx + 1
  call mpi_allreduce(pe_mey,pe_ny,1,MPI_INTEGER,MPI_MAX,pe_medomm,ierr)
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
