program test_spread
  use ISO_C_BINDING
  implicit none
  include 'mpif.h'
  interface
    function RPN_COMM_spread_context(context,com,rootpe,pe,npts) result(status)
    use ISO_C_BINDING
    implicit none
      type(c_ptr), intent(OUT) :: context
      character (len=*), intent(IN) :: com             ! RPN_COMM communicator
      integer, intent(IN) :: npts                      ! number of data points
      integer, intent(IN) :: rootpe                    ! root PE for the spread operation
      integer, dimension(npts), intent(IN) :: pe       ! destination table, data point i will be sent to PE pe(i)
      integer :: status                                ! 0 if successful, non zero otherwise
    end function RPN_COMM_spread_context
    function RPN_COMM_spread(context,source,npts,ndata,dest) result(status)
    use ISO_C_BINDING
    implicit none
      type(c_ptr), intent(IN) :: context
      integer, intent(IN) :: npts, ndata
      real, dimension(npts,ndata), intent(IN) :: source  ! source array, used only on root PE
      real, dimension(:,:), pointer, intent(INOUT) :: dest
      integer :: status
    end function RPN_COMM_spread
  end interface

  integer, PARAMETER :: npts = 10
  integer :: npes, myrank
  integer, dimension(npts) :: pe
  real, dimension(npts,1) :: source
  real, dimension(npts,5) :: source2
  real, dimension(npts) :: rpe
  integer :: i, status, ierr, j
  type(c_ptr) :: context
  real, dimension(:,:), pointer, save :: dest => NULL()
  real, dimension(:,:), pointer, save :: dest2 => NULL()

  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD,npes,ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,myrank,ierr)
  do i=1,npts
    source(i,:)=i
    source2(i,:)=i
  enddo
  call random_number(rpe)
  pe = nint( rpe*(npes-1) )
!  pe(5:7) = -1
  pe( (/1,3,9/) ) = -1
  if(myrank==0) write(0,*)'PE=',pe
  status = RPN_COMM_spread_context(context,'GRID',0,pe,npts)
!  goto 111
  status = RPN_COMM_spread(context,source,npts,size(source,2),dest)
  if(associated(dest)) then
    write(0,*)'dest is now associated and has dimensions',shape(dest)
  endif
  write(0,*)'--------------------------------------------'
!  status = RPN_COMM_spread(context,source,npts,size(source,2),dest)
111 continue
  write(0,*)'--------------------------------------------'
    do j=npts,1,-1
!    write(0,*)'SOURCE data:',nint(source2(j,:))
    enddo
  status = RPN_COMM_spread(context,source2,npts,size(source2,2),dest2)
  if(associated(dest2)) then
    write(0,*)'dest2 is now associated and has dimensions',shape(dest2)
  endif
  if(status >= 0)  write(0,*)'dest2 received',status,' data point(s)'
  if(status <  0)  write(0,*)'ERROR while dest2 was receiving data'
  
  write(0,*)'--------------------------------------------'
!  status = RPN_COMM_spread(context,source,npts,size(source2,2),dest2)
  write(0,*)'============================================='
  call mpi_finalize()
  stop
end program test_spread
integer function RPN_COMM_comm(str)
include 'mpif.h'
character (len=*) :: str
RPN_COMM_comm = MPI_COMM_WORLD
return
end
subroutine xmpi_comm_size(comm,siz,ierr)
integer :: comm,siz,ierr
siz = 3
return
end
subroutine xmpi_comm_rank(comm,rank,ierr)
integer :: comm,rank,ierr
rank=0
return
end
subroutine xmpi_scatter()
return
end