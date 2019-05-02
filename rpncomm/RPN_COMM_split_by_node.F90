!/! RPN_COMM - Library of useful routines for C and FORTRAN programming
! ! Copyright (C) 2019  Division de Recherche en Prevision Numerique
! !                     Environnement Canada
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

!****P* rpn_comm/commnicators  (communicator management package)
!******
!InTf!
!****f* rpn_comm/RPN_COMM_split_by_node  split a communicator on a host basis
! SYNOPSIS
subroutine RPN_COMM_split_by_node(oldcomm, newcomm, rank, isiz, err)   !InTf!
! AUTHOR
!  M.Valin Recherche en Prevision Numerique 2015
! IGNORE
  use ISO_C_BINDING
  implicit none
! ARGUMENTS
  integer, intent(IN)  :: oldcomm   ! MPI communicator to split on a host basis        !InTf!
  integer, intent(OUT) :: newcomm   ! newcommunicator to be used py PEs on same host   !InTf!
  integer, intent(OUT) :: rank      ! rank in new communicator                         !InTf!
  integer, intent(OUT) :: isiz      ! size of new communicator                         !InTf!
  integer, intent(OUT) :: err       ! error code                                       !InTf!
!******
  include 'mpif.h'
  interface
    function gethostid() result(id) BIND(C,name='gethostid')
      import :: C_LONG
      integer(C_LONG) :: id
    end function gethostid
  end interface
  include 'RPN_COMM_constants.inc'

  integer, parameter :: MAX_CACHE=16
  integer :: myhost, myhost0, myhost1, tmpcomm, i
  integer, save :: ncached = 0
  integer, dimension(MAX_CACHE) :: cold, cnew

  err = RPN_COMM_ERROR      ! precondition for failure
  rank = -1
  isiz = 0
  newcomm = MPI_COMM_NULL

  do i = 1, ncached
    if(cold(i) == oldcomm) newcomm = cnew(i)  ! cached entry found
  enddo

  if(newcomm == MPI_COMM_NULL) then          ! nothing useful found in cache
    call mpi_comm_rank(oldcomm, rank, err)
    myhost  = gethostid()                    ! host id
    myhost0 = iand(myhost , Z'7FFFFFFF')     ! lower 31 bits
    myhost1 = iand( ishft(myhost, -31) , 1 ) ! upper bit

    call MPI_Comm_split(oldcomm , myhost0, rank, tmpcomm, err)     ! split oldcomm using the lower 31 bits of host id , weight=rank in base
    if(err .ne. MPI_SUCCESS) return
    call MPI_Comm_split(tmpcomm ,myhost1, rank, newcomm, err)     ! re split using the upper bit of host id , weight=rank in base
    if(err .ne. MPI_SUCCESS) return
  endif
  if(ncached < MAX_CACHE) then                ! add to cache if cache is not full
    ncached = ncached + 1
    cold(ncached) = oldcomm
    cnew(ncached) = newcomm
  endif

  call MPI_Comm_rank(newcomm, rank,err);                         ! rank of this PE on this SMP node
  if(err .ne. MPI_SUCCESS) return
  call MPI_Comm_size(newcomm, isiz, err);                        ! number of PEs on this SMP node
  if(err .ne. MPI_SUCCESS) return
  
  
end subroutine RPN_COMM_split_by_node  !InTf!
