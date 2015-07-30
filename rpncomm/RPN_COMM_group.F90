!/! RPN_COMM - Library of useful routines for C and FORTRAN programming
! ! Copyright (C) 1975-2015  Division de Recherche en Prevision Numerique
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
      integer function RPN_COMM_group(com)         !InTf!
!      Luc Corbeil, 2000-11-21
!
!      lien entre chaine de caractere de groupe
!      GRID, EW et NS et leur numero assigne par
!      MPI.
!
      use rpn_comm
      implicit none                                 !InTf!
!      include mpif.h
!      include rpn_comm.h
      character(len=*), intent(IN) :: com           !InTf!
      character(len=32) comm
      integer ierr,world_group
      call rpn_comm_low2up(com,comm)

      if (comm(1:9) == RPN_COMM_GRIDPEERS) then
         RPN_COMM_group=pe_gr_grid_peers
         return
      endif
      if (comm(1:4) == RPN_COMM_GRID) then
         RPN_COMM_group=pe_gr_indomm
         return
      endif
      if (comm(1:4) == RPN_COMM_DOMM) then
         RPN_COMM_group=pe_gr_indomm
         return
      endif
      if (comm(1:5) == RPN_COMM_WORLD) then
         call MPI_COMM_GROUP(WORLD_COMM_MPI,world_group,ierr)
         RPN_COMM_group=world_group
         return
      endif
      if (comm(1:10) == RPN_COMM_ALLDOMAINS) then
         call MPI_COMM_GROUP(WORLD_COMM_MPI,world_group,ierr)
         RPN_COMM_group=world_group
         return
      endif
      if (comm(1:8) == RPN_COMM_UNIVERSE) then
         call MPI_COMM_GROUP(WORLD_COMM_MPI,world_group,ierr)
         RPN_COMM_group=world_group
         return
      endif
      if (comm(1:9) == RPN_COMM_ALLGRIDS) then
         RPN_COMM_group=pe_gr_a_domain
         return
      endif
      if (comm(1:9) == RPN_COMM_MULTIGRID) then
         RPN_COMM_group=pe_gr_indomms
         return
      endif
      if (comm(1:3) == RPN_COMM_ALL) then
         RPN_COMM_group=pe_gr_wcomm
         return
      endif
      if(comm(1:4) == RPN_COMM_DEFAULT) then
         RPN_COMM_group=pe_defgroup
         return
      endif
      if (comm(1:2) == RPN_COMM_EW) then
         RPN_COMM_group=pe_gr_myrow
         return
      endif
      if (comm(1:2) == RPN_COMM_NS) then
         RPN_COMM_group=pe_gr_mycol
         return
      endif
      if (comm(1:10) == RPN_COMM_BLOCMASTER) then
         RPN_COMM_group=pe_gr_blocmaster
         return
      endif
      if (comm(1:4) == RPN_COMM_BLOCK) then
         RPN_COMM_group=pe_gr_bloc
         return
      endif


      write(rpn_u,*) 'Unknown group ',com,', aborting'
        stop
        
      return
      end function RPN_COMM_group                  !InTf!
