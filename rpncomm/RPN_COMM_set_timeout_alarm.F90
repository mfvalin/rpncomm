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
      function RPN_COMM_set_timeout_alarm(seconds) result(seconds_since)  !InTf!
      use ISO_C_BINDING
      implicit none
      integer, intent(IN) :: seconds  !InTf!
      integer :: seconds_since  !InTf!

      interface
      function c_alarm(seconds) result(seconds_since) BIND(C,name='alarm')
        use ISO_C_BINDING
        implicit none
        integer(C_INT), intent(IN), value :: seconds
        integer(C_INT) :: seconds_since
      end function c_alarm
      end interface

      seconds_since = c_alarm(seconds)
!      print *,'alarm set to ',seconds,' seconds'
      return
      end function RPN_COMM_set_timeout_alarm                             !InTf!
