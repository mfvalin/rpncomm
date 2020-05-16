  include 'RPN_COMM_mpi_symbols.inc'
  interface
    subroutine RPN_COMM_get_mpi_definitions(what, ierr)
      import :: RPN_COMM_mpi_definitions
      type(RPN_COMM_mpi_definitions), intent(INOUT) :: what
      integer, intent(OUT) :: ierr
    end subroutine RPN_COMM_get_mpi_definitions
  end interface
