      ! pmem is array pointing to shared memory segment
      ! 4            : number of values in array
      ! pe_server_id : server kind
      ! pe_grid      : communicator for servers
      ! pe_me_grid   : server ordinal in pe_grid
      ! mbytes       : size in MBytes of array pmem
      ! open calling sequence, first element of array gives number of values
      subroutine RPN_COMM_io_server(pmem,parms) ! test/demo I/O server routine used to test rpn_comm_init
        implicit none
        character(len=1), intent(INOUT), dimension(*) :: pmem
        integer, intent(IN), dimension(*) :: parms
        print *,'entering IO server',parms(4)
        if(parms(1) .ne. 4) then
          print *,'ERROR: wrong version of rpn_comm_init calling this test IO server'
          return
        endif
        return
      end subroutine RPN_COMM_io_server
