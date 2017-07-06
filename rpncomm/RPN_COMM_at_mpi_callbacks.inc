       interface
        function at_mpi_init(func,arg) bind(C,name='AtMpiInit')                &
     &           result(status)
          import :: C_FUNPTR, C_PTR, C_INT
          type(C_FUNPTR), intent(IN) :: func
          type(C_PTR), intent(IN) :: arg
          integer(C_INT) :: status
        end function at_mpi_init
        function at_mpi_finalize(func,arg) bind(C,name='AtMpiFinalize')        &
     &           result(status)
          import :: C_FUNPTR, C_PTR, C_INT
          type(C_FUNPTR), intent(IN) :: func
          type(C_PTR), intent(IN) :: arg
          integer(C_INT) :: status
        end function at_mpi_finalize
      end interface
