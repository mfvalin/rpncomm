#define IN_RPN_COMM_mpi_callbacks

module RPN_COMM_init_callbacks
  use ISO_C_BINDING
  save
  public
  integer, parameter :: MAXFNS = 8
  type(C_FUNPTR), dimension(MAXFNS) :: funcs    ! callback function pointers
  type(C_PTR), dimension(MAXFNS) :: args        ! argument address for callback function
  integer(C_INT) :: nf=0                        ! number of registered functions
  abstract interface                            ! abstract interface for callback functions
    subroutine func(arg)
      integer, dimension(*) :: arg
    end subroutine func
  end interface
contains
! register callbacks (up to MAXFNS) to be called AFTER mpi_init
  function at_mpi_init(fn,arg) bind(C,name='AtMpiInit') result(status)   !InTf!
!!  import :: C_FUNPTR, C_PTR, C_INT                                     !InTf!
    implicit none
    type(C_FUNPTR), intent(IN) :: fn                                     !InTf!
    type(C_PTR), intent(IN) :: arg                                       !InTf!
    integer(C_INT) :: status                                             !InTf!

    if(nf >= MAXFNS) then  ! OOPS, too many functions registered
      status = -1
      return
    endif
    nf = nf + 1
    funcs(nf) = fn
    args(nf) = arg
    status = MAXFNS - nf   ! return number of slots left to register functions
  end function at_mpi_init                                               !InTf!
end module RPN_COMM_init_callbacks

module RPN_COMM_finalize_callbacks
  use ISO_C_BINDING
  save
  public
  integer, parameter :: MAXFNS = 8
  type(C_FUNPTR), dimension(MAXFNS) :: funcs    ! callback function pointers
  type(C_PTR), dimension(MAXFNS) :: args        ! argument address for callback function
  integer(C_INT) :: nf=0                        ! number of registered functions
  abstract interface                            ! abstract interface for callback functions
    subroutine func(arg)
      integer, dimension(*) :: arg
    end subroutine func
  end interface
contains
! register callbacks (up to MAXFNS) to be called BEFORE mpi_finalize
  function at_mpi_finalize(fn,arg) bind(C,name='AtMpiFinalize') result(status)   !InTf!
!!  import :: C_FUNPTR, C_PTR, C_INT                                             !InTf!
    implicit none
    type(C_FUNPTR), intent(IN) :: fn                                             !InTf!
    type(C_PTR), intent(IN) :: arg                                               !InTf!
    integer(C_INT) :: status                                                     !InTf!

    if(nf >= MAXFNS) then  ! OOPS, too many functions registered
      status = -1
      return
    endif
    nf = nf + 1
    funcs(nf) = fn
    args(nf) = arg
    status = MAXFNS - nf   ! return number of slots left to register functions
  end function at_mpi_finalize                                                   !InTf!
end module RPN_COMM_finalize_callbacks

subroutine mpi_init(ierr)        ! overloading MPI library init function wrapper
  implicit none
  integer, intent(OUT) :: ierr
  call RPN_COMM_init_wrapper(ierr)
end subroutine mpi_init

subroutine RPN_COMM_init_wrapper(ierr)
  use RPN_COMM_init_callbacks
  implicit none
  integer, intent(OUT) :: ierr
  procedure(func), pointer :: my_func
  integer, dimension(:), pointer :: my_arg
  integer :: i
  integer, external :: c_len

  call PMPI_Init(ierr)           ! call the REAL MPI library init function
  do i = 1, nf                   ! then call all functions registered via at_mpi_init (in order)
    call c_f_procpointer(funcs(i),my_func)    ! pointer to function
    call c_f_pointer(args(i),my_arg,[1])      ! pointer to argument
    call my_func(my_arg)                      ! call function with argument
  enddo
  return
end subroutine RPN_COMM_init_wrapper

subroutine mpi_finalize(ierr)        ! overloading MPI library finalize function wrapper
  implicit none
  integer, intent(OUT) :: ierr
  call RPN_COMM_finalize_wrapper(ierr)
end subroutine mpi_finalize

subroutine RPN_COMM_finalize_wrapper(ierr)
  use RPN_COMM_finalize_callbacks
  implicit none
  integer, intent(OUT) :: ierr
  procedure(func), pointer :: my_func
  integer, dimension(:), pointer :: my_arg
  integer :: i

  do i = nf, 1, -1                  ! call all functions registered via at_mpi_finalize (in reverse order)
    call c_f_procpointer(funcs(i),my_func)    ! pointer to function
    call c_f_pointer(args(i),my_arg,[1])      ! pointer to argument
    call my_func(my_arg)                      ! call function with argument
  enddo
  call PMPI_Finalize(ierr)           ! then call the REAL MPI library finalize function
  return
end subroutine RPN_COMM_finalize_wrapper

#if defined(SELF_TEST)
program test_init_mpi   ! et le petit test qui va avec
  use ISO_C_BINDING
  implicit none
  include 'mpif.h'
  include 'RPN_COMM.inc'
  external :: my_sub, my_sub2
  integer :: ierr, i
  integer, dimension(10) :: my_arg, my_arg2, status1, status2

  do i = 1, 10
    my_arg(i) = 100 + i
    my_arg2(i) = 200 + i
    status1(i) = at_mpi_init(c_funloc(my_sub),c_loc(my_arg(i)))
    status2(i) = at_mpi_finalize(c_funloc(my_sub2),c_loc(my_arg2(i)))
  enddo

  call mpi_init(ierr)
  print 1,'ierr, status1 =',ierr, status1
  print 1,'ierr, status2 =',ierr, status2
  call mpi_finalize(ierr)
  stop
1 format(A,11I3)
end

subroutine my_sub(arg)
  integer, dimension(*) :: arg
  print *,'in my_sub, init, arg =',arg(1)
end subroutine my_sub

subroutine my_sub2(arg)
  integer, dimension(*) :: arg
  print *,'in my_sub2, finalize, arg =',arg(1)
end subroutine my_sub2
#endif
