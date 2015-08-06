
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
module rpncomm_com
  use rpn_comm
  implicit none
  integer, parameter :: MAX_COMM_TAB=128
  type(symtab), dimension(:), pointer, save :: com_tab => NULL()
  integer, save :: defcom_index=-1
  integer, save :: max_com_index=0
contains
!
! allocate and initialize com_tab, the communicator table
!
  subroutine init_com_tab
    implicit none
    integer :: i

    if( associated(com_tab) ) return  ! job already done

    allocate(com_tab(MAX_COMM_TAB))
    com_tab( 1) = symtab(pe_grid_peers,RPN_COMM_GRIDPEERS)
    com_tab( 2) = symtab(pe_indomm,RPN_COMM_GRID)
    com_tab( 3) = symtab(pe_indomm,RPN_COMM_DOMM)
    com_tab( 4) = symtab(WORLD_COMM_MPI,RPN_COMM_WORLD)
    com_tab( 5) = symtab(WORLD_COMM_MPI,RPN_COMM_ALLDOMAINS)
    com_tab( 6) = symtab(pe_a_domain,RPN_COMM_ALLGRIDS)
    com_tab( 7) = symtab(pe_indomms,RPN_COMM_MULTIGRID)
    com_tab( 8) = symtab(pe_wcomm,RPN_COMM_ALL)
    com_tab( 9) = symtab(pe_defcomm,RPN_COMM_DEFAULT)  ; defcom_index = 9  ! this item might get updated by rpn_comm_defo
    com_tab(10) = symtab(pe_myrow,RPN_COMM_EW)
    com_tab(11) = symtab(pe_mycol,RPN_COMM_NS)
    com_tab(12) = symtab(pe_blocmaster,RPN_COMM_BLOCMASTER)
    com_tab(13) = symtab(pe_bloc,RPN_COMM_BLOCK)
    com_tab(14) = symtab(MPI_COMM_WORLD,RPN_COMM_UNIVERSE)
    do i = 15,MAX_COMM_TAB
      com_tab(i) = symtab(MPI_COMM_NULL,"")
    enddo
  end subroutine init_com_tab
!
! get index of string com in com_tab
!
  function indx_com_tab(com) result(indx)
    implicit none
    character(len=*), intent(IN) :: com
    integer :: indx
    integer :: i

    if(.not. associated(com_tab)) call init_com_tab
    indx = MAX_COMM_TAB
    do i = 1,max_com_index
      if( trim(com_tab(i)%string) == trim(com) ) then
        indx = i
        return
      endif
    enddo
  end function indx_com_tab
end module rpncomm_com
!InTf!
      integer function RPN_COMM_comm(com)           !InTf!
!	Luc Corbeil, 2000-11-21
!
!	lien entre chaine de caractere de communicateur
!	GRID, EW et NS et leur numero assigne par
!	MPI.
!
      use rpncomm_com
      implicit none                                 !InTf!
!      include mpif.h
!        include rpn_comm.h
      character(len=*), intent(IN) :: com           !InTf!
      character(len=32) comm
      integer :: i

      call rpn_comm_low2up(com,comm)

      RPN_COMM_comm = com_tab(indx_com_tab(comm))%number

      return

      RPN_COMM_comm = MPI_COMM_NULL

      if (trim(comm) == RPN_COMM_GRIDPEERS) then
         RPN_COMM_comm=pe_grid_peers
         return
      endif
      if (trim(comm) == RPN_COMM_GRID) then
         RPN_COMM_comm=pe_indomm  ! alias pe_grid
         return
      endif
      if(trim(comm) == RPN_COMM_DOMM) then
         RPN_COMM_comm=pe_indomm  ! alias pe_grid
         return
      endif
      if (trim(comm) == RPN_COMM_WORLD) then
         RPN_COMM_comm=WORLD_COMM_MPI  ! alias pe_all_domains
         return
      endif
      if(trim(comm) == RPN_COMM_ALLDOMAINS) then
         RPN_COMM_comm=WORLD_COMM_MPI  ! alias pe_all_domains
         return
      endif
      if (trim(comm) == RPN_COMM_ALLGRIDS) then
         RPN_COMM_comm=pe_a_domain
         return
      endif
      if (trim(comm) == RPN_COMM_MULTIGRID) then
         RPN_COMM_comm=pe_indomms ! alias pe_multi_grid
         return
      endif
      if (trim(comm) == RPN_COMM_ALL) then
        RPN_COMM_comm=pe_wcomm  ! alias pe_grid
        return
      endif
      if(trim(comm) == RPN_COMM_DEFAULT) then
        RPN_COMM_comm=pe_defcomm
        return
      endif
      if (trim(comm) == RPN_COMM_EW) then
         RPN_COMM_comm=pe_myrow
         return
      endif
      if (trim(comm) == RPN_COMM_NS) then
         RPN_COMM_comm=pe_mycol
         return
      endif
      if (trim(comm) == RPN_COMM_BLOCMASTER) then
         RPN_COMM_comm=pe_blocmaster
         return
      endif
      if (trim(comm) == RPN_COMM_BLOCK) then
         RPN_COMM_comm=pe_bloc
         return
      endif
      if (trim(comm) == RPN_COMM_UNIVERSE) then
         RPN_COMM_comm=MPI_COMM_WORLD
         return
      endif

      write(rpn_u,*) 'Unknown communicator ',comm,', aborting'
      stop
      return
      end function RPN_COMM_comm                                  !InTf!
!InTf!
      integer function RPN_COMM_custom_comm(com,name,mode)        !InTf!
      use rpn_comm
      implicit none                                               !InTf!
!     lookup or create a communicator with a rpn_comm style name
      character(len=*), intent(IN) :: name                        !InTf!
      integer, intent(IN) :: com                                  !InTf!
      integer, intent(IN) :: mode                                 !InTf!
!
      integer :: i
      character (len=32) :: name2
!
      integer, parameter :: MAX_NAMES=128
      type(SYMTAB), save, dimension(:), pointer :: names => NULL()
      integer, save :: entries=0
!
      if(.not. associated(names)) allocate(names(MAX_NAMES))

      RPN_COMM_custom_comm=MPI_COMM_NULL
      name2 = trim(name)
!
      if(mode==RPN_COMM_GET) then                 ! look for rpn_comm communicator named "name"
         do i = 1 , entries
            if(trim(name2)==trim(names(i)%string)) then
               RPN_COMM_custom_comm = names(i)%number
               return
            endif
         enddo
      else if(mode==RPN_COMM_SET) then             ! add "name" and com to the rpn_comm communicators
         if(entries<MAX_NAMES) then
            entries = entries + 1
            names(entries)%string = trim(name2)
            names(entries)%number = com
            RPN_COMM_custom_comm=com
         else
            write(rpn_u,*) 'ERROR: communicator table full'
         endif
      else if(mode==RPN_COMM_DEL) then              ! delete "name" and com from rpn_comm communicators
      else
         write(rpn_u,*) 'ERROR: RPN_COMM_custom_comm illegal mode'
      endif
      return
      end function RPN_COMM_custom_comm                      !InTf!
!
!       fill an entity of type rpncomm_communicator from type string
!       ctyp_c : character version of communicator
!       ctyp   : new item of type rpncomm_communicator
!InTf!
        subroutine RPN_COMM_i_comm(ctyp_c,ctyp)            !InTf!
        use rpn_comm
!!      import :: rpncomm_communicator                       !InTf!
        implicit none
        type(rpncomm_communicator), intent(OUT) :: ctyp      !InTf!
        character(len=*), intent(IN) :: ctyp_c               !InTf!
        integer, external :: RPN_COMM_comm

        ctyp%p = C_LOC(WORLD_COMM_MPI)            ! signature
        ctyp%t2 = RPN_COMM_comm(ctyp_c)           ! communicator value
        ctyp%t1 = -ctyp%t2                        ! - communicator value
        end subroutine RPN_COMM_i_comm                    !InTf!

