      ! pmem is array pointing to shared memory segment
      ! 104          : number of values in array + version_number * 100
      ! pe_server_id : server kind
      ! pe_grid      : communicator for servers
      ! pe_me_grid   : server ordinal in pe_grid
      ! mbytes       : size in MBytes of array pmem
      ! open calling sequence, first element of array gives number of values
      subroutine RPN_COMM_io_server(pmem,parms) ! test/demo I/O server routine used to test rpn_comm_init
        implicit none
        character(len=1), intent(INOUT), dimension(*) :: pmem
        integer, intent(IN), dimension(*) :: parms
        integer :: pe_me, ierr
        if(parms(1) .ne. 1004) then
          print *,'ERROR: version mismatch between rpn_comm_init and the IO server'
          print *,'       ( expected',1004,' got',parms(1),' )'
          return
        endif
        call MPI_comm_rank(parms(3),pe_me,ierr)
        print 100,'entering IO server type ',parms(2),', global server rank= ',parms(4)
100     format(A,I2,A,I3)
        return
      end subroutine RPN_COMM_io_server
