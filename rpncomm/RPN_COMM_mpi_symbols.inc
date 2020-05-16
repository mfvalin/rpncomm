  !     user defined type containing frequently used MPI symbols
  integer, parameter :: mpi_symbols_version = 123456
  type :: RPN_COMM_mpi_definitions
    integer :: version = mpi_symbols_version
    !
    !     NULL "handles"
    !
    integer :: MPI_GROUP_NULL, MPI_COMM_NULL, MPI_DATATYPE_NULL
    integer :: MPI_REQUEST_NULL, MPI_OP_NULL, MPI_ERRHANDLER_NULL
    integer :: MPI_INFO_NULL, MPI_WIN_NULL
    !
    !     varia
    !
    integer :: MPI_STATUS_SIZE
    integer :: MPI_ANY_SOURCE, MPI_ANY_TAG
    integer :: MPI_SUCCESS, MPI_ERROR

    integer :: MPI_COMM_WORLD, MPI_COMM_SELF
    integer :: MPI_GROUP_EMPTY
    integer :: MPI_ERRORS_ARE_FATAL, MPI_ERRORS_RETURN
    !
    !     data types
    !
    integer :: MPI_BYTE, MPI_PACKED, MPI_UB, MPI_LB
    integer :: MPI_CHARACTER, MPI_LOGICAL
    integer :: MPI_INTEGER, MPI_INTEGER1, MPI_INTEGER2, MPI_INTEGER4
    integer :: MPI_INTEGER8, MPI_INTEGER16
    integer :: MPI_REAL, MPI_REAL2, MPI_REAL4, MPI_REAL8, MPI_REAL16
    integer :: MPI_DOUBLE_PRECISION
    integer :: MPI_COMPLEX, MPI_COMPLEX8, MPI_COMPLEX16, MPI_COMPLEX32
    integer :: MPI_DOUBLE_COMPLEX
    integer :: MPI_2REAL, MPI_2DOUBLE_PRECISION, MPI_2INTEGER
    integer :: MPI_2COMPLEX, MPI_2DOUBLE_COMPLEX
    integer :: MPI_LOGICAL1, MPI_LOGICAL2, MPI_LOGICAL4, MPI_LOGICAL8
    !
    !     operators
    !
    integer :: MPI_MAX, MPI_MIN, MPI_SUM, MPI_PROD, MPI_LAND
    integer :: MPI_BAND, MPI_LOR, MPI_BOR, MPI_LXOR, MPI_BXOR
    integer :: MPI_MAXLOC, MPI_MINLOC, MPI_REPLACE
  end type
