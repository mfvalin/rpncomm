  include 'RPN_COMM_mpi_symbols.inc'
  include 'RPN_COMM_mpi_layout.inc'
  interface
    subroutine RPN_COMM_get_mpi_definitions(what, ierr)
      import :: RPN_COMM_mpi_definitions
      type(RPN_COMM_mpi_definitions), intent(INOUT) :: what
      integer, intent(OUT) :: ierr
    end subroutine RPN_COMM_get_mpi_definitions
    subroutine RPN_COMM_get_mpi_layout(what, ierr)
      import :: mpi_layout, C_INT
      type(mpi_layout), intent(INOUT) :: what
      integer(C_INT), intent(OUT) :: ierr
    end subroutine RPN_COMM_get_mpi_layout
  end interface
