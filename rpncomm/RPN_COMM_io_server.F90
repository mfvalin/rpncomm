      ! pmem is array pointing to shared memory segment
      ! 1006  : number of values in array + version_number * 1000
      ! (2)   : server kind
      ! (3)   : communicator for servers
      ! (4)   : server ordinal in pe_grid
      ! (5)   : size in MBytes of array pmem
      ! (6)   : number of Compute PEs in group
      ! (7)   : number of Server PEs in group
      ! open calling sequence, first element of array gives number of values
      subroutine RPN_COMM_io_server(pmem,parms) ! test/demo I/O server routine used to test rpn_comm_init
        use ISO_C_BINDING
        implicit none
        integer, intent(INOUT), dimension(*) :: pmem
        integer, intent(IN), dimension(*) :: parms
        integer :: pe_me, ierr
        type(C_PTR) :: ptr
        integer, dimension(:), pointer :: array
        integer, parameter :: COMMAND_BUFFER_SIZE = 4096
        type :: command_channel
          integer :: id, first, in, out, limit
          integer, dimension(COMMAND_BUFFER_SIZE) :: buf
        end type
        type(command_channel), dimension(:,:), pointer   :: cmem => NULL()   ! communication buffers between Compute and Server processes
        integer :: i,j

        ptr = C_LOC(pmem(1))
        call C_F_POINTER(ptr,array,[10])
        if(parms(1) .ne. 1006) then
          print *,'ERROR: version mismatch between rpn_comm_init and the IO server'
          print *,'       ( expected',1006,' got',parms(1),' )'
          return
        endif
        call MPI_comm_rank(parms(3),pe_me,ierr)
        print 100,'entering IO server type ',parms(2),', global server rank= ',parms(4)
        call C_F_POINTER(ptr,cmem,[parms(6),parms(7)])
        print 100,'Command buffers for server type ',parms(2)
        j = parms(2) + 1
        do i = 1, parms(6)
          print 101,'command buffer : CPE KIND GID FIRST IN OUT LIMIT =',i-1,j-1,cmem(i,j)%id,cmem(i,j)%first,cmem(i,j)%in,cmem(i,j)%out,cmem(i,j)%limit
        enddo
100     format(A,I2,A,I3)
101     format(A,10I5)
        return
      end subroutine RPN_COMM_io_server
