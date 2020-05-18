!/! RPN_COMM - Library of useful routines for C and FORTRAN programming
! ! Copyright (C) 1975-2020  Division de Recherche en Prevision Numerique
! !                          Environnement Canada
! !
! ! This library is free software; you can redistribute it and/or
! ! modify it under the terms of the GNU Lesser General Public
! ! License as published by the Free Software Foundation,
! ! version 2.1 of the License.
! !
! ! This library is distributed in the hope that it will be useful,
! ! but WITHOUT ANY WARRANTY; without even the implied warranty of
! ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! ! Lesser General Public License for more details.
! !
! ! You should have received a copy of the GNU Lesser General Public
! ! License along with this library; if not, write to the
! ! Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! ! Boston, MA 02111-1307, USA.
! !/
subroutine RPN_COMM_get_mpi_definitions(what, ierr)  ! get a copy of MPI definitions
  implicit none
  include 'mpif.h'
  include 'RPN_COMM_mpi_symbols.inc'
  type(RPN_COMM_mpi_definitions), intent(INOUT) :: what
  integer, intent(OUT) :: ierr

  ierr = MPI_ERROR
  if(what%version .ne. mpi_symbols_version) then
    write(0,*)'ERROR: (RPN_COMM_get_mpi_definitions) version mismatch, expected :',mpi_symbols_version,' got :',what%version
    what%version = -1
    return
  endif
  ierr = MPI_SUCCESS
  what%MPI_GROUP_NULL = MPI_GROUP_NULL
  what%MPI_COMM_NULL = MPI_COMM_NULL
  what%MPI_DATATYPE_NULL = MPI_DATATYPE_NULL
  what%MPI_REQUEST_NULL = MPI_REQUEST_NULL
  what%MPI_OP_NULL = MPI_OP_NULL
  what%MPI_ERRHANDLER_NULL = MPI_ERRHANDLER_NULL
  what%MPI_INFO_NULL = MPI_INFO_NULL
  what%MPI_WIN_NULL = MPI_WIN_NULL
  what%MPI_STATUS_SIZE = MPI_STATUS_SIZE
  what%MPI_ANY_SOURCE = MPI_ANY_SOURCE
  what%MPI_ANY_TAG = MPI_ANY_TAG
  what%MPI_SUCCESS = MPI_SUCCESS
  what%MPI_ERROR = MPI_ERROR
  what%MPI_COMM_WORLD = MPI_COMM_WORLD
  what%MPI_COMM_SELF = MPI_COMM_SELF
  what%MPI_GROUP_EMPTY = MPI_GROUP_EMPTY
  what%MPI_ERRORS_ARE_FATAL = MPI_ERRORS_ARE_FATAL
  what%MPI_ERRORS_RETURN = MPI_ERRORS_RETURN
  what%MPI_BYTE = MPI_BYTE
  what%MPI_PACKED = MPI_PACKED
  what%MPI_UB = MPI_UB
  what%MPI_LB = MPI_LB
  what%MPI_CHARACTER = MPI_CHARACTER
  what%MPI_LOGICAL = MPI_LOGICAL
  what%MPI_INTEGER = MPI_INTEGER
  what%MPI_INTEGER1 = MPI_INTEGER1
  what%MPI_INTEGER2 = MPI_INTEGER2
  what%MPI_INTEGER4 = MPI_INTEGER4
  what%MPI_INTEGER8 = MPI_INTEGER8
  what%MPI_INTEGER16 = MPI_INTEGER16
  what%MPI_REAL = MPI_REAL
  what%MPI_REAL2 = MPI_REAL2
  what%MPI_REAL4 = MPI_REAL4
  what%MPI_REAL8 = MPI_REAL8
  what%MPI_REAL16 = MPI_REAL16
  what%MPI_DOUBLE_PRECISION = MPI_DOUBLE_PRECISION
  what%MPI_COMPLEX = MPI_COMPLEX
  what%MPI_COMPLEX8 = MPI_COMPLEX8
  what%MPI_COMPLEX16 = MPI_COMPLEX16
  what%MPI_COMPLEX32 = MPI_COMPLEX32
  what%MPI_DOUBLE_COMPLEX = MPI_DOUBLE_COMPLEX
  what%MPI_2REAL = MPI_2REAL
  what%MPI_2DOUBLE_PRECISION = MPI_2DOUBLE_PRECISION
  what%MPI_2INTEGER = MPI_2INTEGER
  what%MPI_2COMPLEX = MPI_2COMPLEX
  what%MPI_2DOUBLE_COMPLEX = MPI_2DOUBLE_COMPLEX
  what%MPI_LOGICAL1 = MPI_LOGICAL1
  what%MPI_LOGICAL2 = MPI_LOGICAL2
  what%MPI_LOGICAL4 = MPI_LOGICAL4
  what%MPI_LOGICAL8 = MPI_LOGICAL8
  what%MPI_MAX = MPI_MAX
  what%MPI_MIN = MPI_MIN
  what%MPI_SUM = MPI_SUM
  what%MPI_PROD = MPI_PROD
  what%MPI_LAND = MPI_LAND
  what%MPI_BAND = MPI_BAND
  what%MPI_LOR = MPI_LOR
  what%MPI_BOR = MPI_BOR
  what%MPI_LXOR = MPI_LXOR
  what%MPI_BXOR = MPI_BXOR
  what%MPI_MAXLOC = MPI_MAXLOC
  what%MPI_MINLOC = MPI_MINLOC
  what%MPI_REPLACE = MPI_REPLACE

  return
end subroutine RPN_COMM_get_mpi_definitions

module RPN_COMM_mpi_layout
  use ISO_C_BINDING
  implicit none
  include 'RPN_COMM_mpi_layout.inc'
  type(mpi_layout), save :: internal_mpi_layout
end module RPN_COMM_mpi_layout

subroutine RPN_COMM_get_mpi_layout(what, ierr) ! get a copy on RPN_COMM internal mpi layout (communicators, ranks, sizes)
  use ISO_C_BINDING
  use RPN_COMM_mpi_layout
  implicit none
  include 'mpif.h'
  type(mpi_layout), intent(INOUT) :: what
  integer(C_INT), intent(OUT) :: ierr
  integer :: i

  ierr = MPI_ERROR
  if(what%version .ne. layout_version) then
    write(0,*)'ERROR: (RPN_COMM_get_mpi_layout) version mismatch, expected :',layout_version,' got :',what%version
    what%version = -1
    return
  endif
  ierr = MPI_SUCCESS

  what = internal_mpi_layout
  what%comm%sgrd%row       = MPI_COMM_NULL   ! row is not defined for supergrids
  what%rank%sgrd%row       = -1
  what%size%sgrd%row       = -1
  what%comm%sgrd%column    = MPI_COMM_NULL   ! and neither is column
  what%rank%sgrd%column    = -1
  what%size%sgrd%column    = -1
!   what%comm%appl%all       = internal_mpi_layout%comm%appl%all
!   what%comm%appl%same_node = internal_mpi_layout%comm%appl%same_node
!   what%comm%appl%same_numa = internal_mpi_layout%comm%appl%same_numa
! 
!   what%comm%sgrd%all       = internal_mpi_layout%comm%sgrd%all
!   what%comm%sgrd%compute   = internal_mpi_layout%comm%sgrd%compute
!   what%comm%sgrd%service   = internal_mpi_layout%comm%sgrd%service
!   what%comm%sgrd%same_node = internal_mpi_layout%comm%sgrd%same_node
!   what%comm%sgrd%same_numa = internal_mpi_layout%comm%sgrd%same_numa
!   what%comm%sgrd%node_peer = internal_mpi_layout%comm%sgrd%node_peer
!   what%comm%sgrd%numa_peer = internal_mpi_layout%comm%sgrd%numa_peer
!   what%comm%sgrd%grid_peer = internal_mpi_layout%comm%sgrd%grid_peer
!   what%comm%sgrd%row       = MPI_COMM_NULL
!   what%comm%sgrd%column    = MPI_COMM_NULL
! 
!   what%comm%grid%all       = internal_mpi_layout%comm%grid%all
!   what%comm%grid%compute   = internal_mpi_layout%comm%grid%compute
!   what%comm%grid%service   = internal_mpi_layout%comm%grid%service
!   what%comm%grid%same_node = internal_mpi_layout%comm%grid%same_node
!   what%comm%grid%same_numa = internal_mpi_layout%comm%grid%same_numa
!   what%comm%sgrd%node_peer = internal_mpi_layout%comm%grid%node_peer
!   what%comm%grid%numa_peer = internal_mpi_layout%comm%grid%numa_peer
!   what%comm%grid%grid_peer = internal_mpi_layout%comm%grid%grid_peer
!   what%comm%grid%row       = internal_mpi_layout%comm%grid%row
!   what%comm%grid%column    = internal_mpi_layout%comm%grid%column
  return
end subroutine RPN_COMM_get_mpi_layout
