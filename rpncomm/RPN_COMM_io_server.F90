      ! pmem is array pointing to shared memory segment
      ! parms(1)   : 1006  = number of values in array + version_number * 1000
      ! parms(2)   : server kind
      ! parms(3)   : communicator for servers
      ! parms(4)   : server ordinal in pe_grid
      ! parms(5)   : size in MBytes of array pmem
      ! parms(6)   : number of Compute PEs in group
      ! parms(7)   : number of Server PEs in group
      ! open calling sequence, first element of array gives number of values
      subroutine RPN_COMM_io_server(pmem,parms) ! test/demo I/O server routine used to test rpn_comm_init
        use ISO_C_BINDING
        use rpn_comm_server_mod
        implicit none
        integer, intent(INOUT), dimension(*) :: pmem
        integer, intent(IN), dimension(*) :: parms
        integer :: pe_me, ierr
        type(C_PTR) :: ptr
        integer, dimension(:), pointer :: array
!         type(command_channel), dimension(:,:), pointer   :: mem => NULL()   ! communication buffers between Compute and Server processes
        type(cmd_channel), dimension(:,:), pointer   :: mem => NULL()   ! communication buffers between Compute and Server processes
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
        call C_F_POINTER(ptr,mem,[parms(7),parms(6)])
        print 100,'Command buffers for server type ',parms(2)
        j = parms(2) + 1
        do i = 1, parms(6)
          print 101,'Server : Client Me GroupID FIRST IN OUT LIMIT =',i-1,j-1,mem(i,j)%id  !  ,mem(i,j)%first,mem(i,j)%in,mem(i,j)%out,mem(i,j)%limit
        enddo
100     format(A,I2,A,I3)
101     format(A,10I5)
        return
      end subroutine RPN_COMM_io_server
