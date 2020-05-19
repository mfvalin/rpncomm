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
!===================================================================
!
!  mpi_layout structure
!  table of elements
!  USAGE :
!       include 'RPN_COMM_mpi_definitions.inc'
!       type(mpi_layout) :: l
!       integer :: ierr
!       call RPN_COMM_get_mpi_layout(l, ierr)
!
!       internal RPN_COMM MPI information is now available, like:
!
!       l%comm%grid%service    ! communicator for service(IO) PEs in a grid
!       l%rank%sgrd%compute    ! rank in a supergrid compute PEs communicator
!       l%comm%grid%compute    ! communicator for compute PEs in a grid
!       l%comm%grid%column     ! grid column communicator
!       l%size%grid%row        ! size of grid row communicator
!
!  version                     ! version marker (consistency check between compilation and runtime)
!  host                        ! host id for this host
!  numa                        ! numa space on host for this PE
!       |
!       | WORLD, all PEs from all applications, compute and service
!       | % wrld % all         ! all PEs in all domains
!       |        % same_node   ! PEs in the same SMP node
!       |        % same_numa   ! PEs in the same NUMA space
!       |
!       | APPLICATION, containing M identical supergrids
!  comm | % appl % all         ! all PEs belonging to an application (model, domain)
!  rank |        % same_node   ! PEs in the same SMP node
!  size |        % same_numa   ! PEs in the same NUMA space
!       |
!       | SUPERGRID, M identical supergrids in a domain, containing N identical grids
!       | % sgrd % all         ! all PEs belonging to a supergrid (compute and service)
!       |        % compute     ! compute PEs belonging to a supergrid
!       |        % service     ! service (IO) PEs belonging to a supergrid
!       |        % same_node   ! PEs belonging to a supergrid in same SMP node
!       |        % same_numa   ! PEs belonging to a supergrid in same NUMA space
!       |        % node_peer   ! PEs belonging to a supergrid with same rank in the SMP nodes
!       |        % numa_peer   ! PEs belonging to a supergrid with same rank in NUMA spaces
!       |        % grid_peer   ! PEs from all supergrids in application with same rank in supergrid
!       |        % row         ! not used for supergrids
!       |        % column      ! not used for supergrids
!       |
!       | GRID, N identical grids in a supergrid
!       | % grid % all         ! all PEs belonging to a grid (compute and service)
!       |        % compute     ! compute PEs belonging to a grid
!       |        % service     ! service (IO) PEs belonging to a grid
!       |        % same_node   ! PEs belonging to a grid in same SMP node
!       |        % same_numa   ! PEs belonging to a grid in same NUMA space
!       |        % node_peer   ! PEs belonging to a grid with same rank in the SMP nodes
!       |        % numa_peer   ! PEs belonging to a grid with same rank in NUMA spaces
!       |        % grid_peer   ! PEs from all grids in application with same rank in grid
!       |        % row         ! compute PEs from a grid in the same row (same j/y ordinal)
!       |        % column      ! compute PEs from a grid in the same column (same i/x ordinal)
!       |
!       | % blck % all         ! all PEs belonging to a block (subgrid) (all are compute PEs)
!       |        % row         ! block PEs in the same row (same j/y ordinal)
!       |        % column      ! block PEs in the same column (same i/x ordinal)
!
! subgrids (blocks of PEs within a grid) , tentative
  type, BIND(C) :: subgrid
    integer(C_INT) :: all              ! all PEs
    integer(C_INT) :: row
    integer(C_INT) :: column
  end type
! grids and supergrids
  type, BIND(C) :: mpigrid
    integer(C_INT) :: all              ! all PEs, compute and service
    integer(C_INT) :: compute          ! compute PEs
    integer(C_INT) :: service          ! service (Input/Output) PEs
    integer(C_INT) :: color            ! ordinal of grid in supergrid or ordinal of supergrid in application
    integer(C_INT) :: same_node        ! PEs in same SMP node (compute and service)
    integer(C_INT) :: same_numa        ! PEs in same numa space (compute and service)
    integer(C_INT) :: node_peer        ! node to node peer (PEs with same rank in node)
    integer(C_INT) :: numa_peer        ! numa space to numa space peer (PEs with same rank in numa space)
    integer(C_INT) :: grid_peer        ! grid to grid / supergrid to supergrid peer (PEs with same rank in grid/supergrid)
    integer(C_INT) :: row              ! undefined in a supergrid, information for a grid row
    integer(C_INT) :: column           ! undefined in a supergrid, information for a grid column
  end type
! applications (domains) and world (less information than grids and supergrids)
  type, BIND(C) :: application
    integer(C_INT) :: all              ! all PEs, compute and service
    integer(C_INT) :: same_node        ! PEs in same SMP node (compute and service)
    integer(C_INT) :: same_numa        ! PEs in same numa space (compute and service)
    integer(C_INT) :: color            ! application id (color)
  end type
! the sets of information
  type, BIND(C) :: grid_hierarchy
    type(application) :: wrld          ! the entire world, conrains all applications
    type(application) :: appl          ! application (domain), contains (M identical supergrids)
    type(mpigrid)     :: sgrd          ! supergrid (contains N identical grids)
    type(mpigrid)     :: grid          ! basic computation grid
    type(subgrid)     :: blck          ! block of PEs inside a grid (subgrid)
  end type
! the three usual associates, communicators, rank, and size
  integer(C_INT), parameter :: layout_version = 10000 ! version 1.00.00
  type, BIND(C) :: mpi_layout
    integer(C_INT) :: version = layout_version
    integer(C_INT) :: host             ! SMP host id
    integer(C_INT) :: numa             ! numa space for this PE
    type(grid_hierarchy) :: comm       ! communicators
    type(grid_hierarchy) :: rank       ! rank in communicator
    type(grid_hierarchy) :: size       ! size of communicator
  end type

