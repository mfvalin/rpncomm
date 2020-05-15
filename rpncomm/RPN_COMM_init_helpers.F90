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
!InTf!
      subroutine RPN_COMM_mydomain (call_back, mydomain)             !InTf!
      use rpn_comm
      implicit none                                                  !InTf!
!
      external :: call_back                                          !InTf!
      integer, intent(OUT) :: mydomain                               !InTf!
!
!      include 'mpif.h'

      logical mpi_started
      integer ierr, err, err_all, pe_tot2, pe_me2, npe_per_domain,&
     &        offset, ndomains
!
      if(RPM_COMM_IS_INITIALIZED) then
         if (pe_me.eq.0) write (6,1002)
         stop
      endif
      call MPI_INITIALIZED(mpi_started,ierr)
      if (.not. mpi_started ) call MPI_init(ierr)

      call mpi_comm_size (MPI_COMM_WORLD,pe_tot2,ierr)
      call mpi_comm_rank (MPI_COMM_WORLD,pe_me2 ,ierr)

      call call_back (ndomains, offset, err)

      call mpi_allreduce &
     &   (err,err_all,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)

      if (err_all.lt.0) then
         if (pe_me2.eq.0) write (6,1001) 
         call rpn_comm_FINALIZE(ierr)
         stop
      endif

      npe_per_domain = pe_tot2 / ndomains
      mydomain       = pe_me2  / npe_per_domain + offset
!
 1001 format &
     &(/'===> rpn_comm_mydomain: FATAL error in call_back --- ABORTING ---'/)
 1002 format &
     &(/'===> rpn_comm_mydomain: RPN_COMM already initialized --- ABORTING ---'/)
      return
      end subroutine RPN_COMM_mydomain                               !InTf!
!===================================================================
!InTf!
      subroutine RPN_COMM_world_set(world_comm)                    !InTf!
      use rpn_comm
      implicit none                                                !InTf!
      integer, intent(IN) ::  world_comm                           !InTf!
!        world_comm=WORLD_COMM_MPI ! should rather be the other way around
      WORLD_COMM_MPI=world_comm
      WORLD_COMM_MPI_INIT=.true.
      return
      end subroutine RPN_COMM_world_set                            !InTf!
!===================================================================
