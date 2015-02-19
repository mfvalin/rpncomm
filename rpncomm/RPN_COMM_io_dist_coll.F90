!/! RMNLIB - Library of useful routines for C and FORTRAN programming
! ! Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
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
subroutine RPN_COMM_shuf_dist(setno,  &
                              global,gni,gnj,gnk,  &
                              local,mini,maxi,minj,maxj,lnk,  &
                              liste_i,liste_o,  &
                              start_x,count_x,start_y,count_y,  &
                              periodx,periody)
  use RPN_COMM_io_pe_tables
  implicit none
  integer, intent(IN) :: setno,gni,gnj,gnk,mini,maxi,minj,maxj,lnk
  integer, intent(IN), dimension(gni,gnj,gnk) :: global
  integer, intent(OUT), dimension(mini:maxi,minj:maxj,lnk) :: local
  integer, intent(IN), dimension(pe_nx) :: start_x,count_x
  integer, intent(IN), dimension(pe_ny) :: start_y,count_y
  integer, intent(IN), dimension(gnk) :: liste_i
  integer, intent(OUT), dimension(lnk) :: liste_o
  logical, intent(IN) :: periodx,periody
  integer :: i

  local = 0
  do i=1,gnk
    liste_o(liste_i(i)) = .false.
  enddo
  do i=1,gnk
    call RPN_COMM_shuf_dist_1(setno)
  enddo
end subroutine RPN_COMM_shuf_dist
!
subroutine RPN_COMM_shuf_dist_1(setno,  &
                              global,gni,gnj,gk,  &
                              local,mini,maxi,minj,maxj,lnk,  &
                              liste_o,  &
                              start_x,count_x,start_y,count_y,  &
                              periodx,periody)
  use rpn_comm
  use RPN_COMM_io_pe_tables
  implicit none
  integer, intent(IN) :: setno,gni,gnj,gk,mini,maxi,minj,maxj,lnk
  integer, intent(IN), dimension(gni,gnj) :: global
  integer, intent(OUT), dimension(mini:maxi,minj:maxj,lnk) :: local
  integer, intent(IN), dimension(pe_nx) :: start_x,count_x
  integer, intent(IN), dimension(pe_ny) :: start_y,count_y
  integer, intent(OUT), dimension(lnk) :: liste_o
  logical, intent(IN) :: periodx,periody

  liste_o(gk) = .false.
  local = 0
end subroutine RPN_COMM_shuf_dist_1
!
subroutine RPN_COMM_shuf_coll()
  use rpn_comm
  use RPN_COMM_io_pe_tables
  implicit none
end subroutine RPN_COMM_shuf_coll
!
subroutine RPN_COMM_shuf_coll_1()
  use rpn_comm
  use RPN_COMM_io_pe_tables
  implicit none
end subroutine RPN_COMM_shuf_coll_1
!
