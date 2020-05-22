!/* RPN_COMM - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2012  Division de Recherche en Prevision Numerique
! *                          Environnement Canada
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
! */
!	FORTRAN to C bridging routines
!
!	RPN_COMM_chdir : bridge to f_RPN_COMM_chdir, null terminates string before calling C routine
!
! may 2020 : use BIND(C) to access chdir, gethostid, and environment variables
!            this obsoletes RPN_COMM_getenv.c, RPN_COMM_chdir.c
!InTf!
      integer function RPN_COMM_chdir(string) !InTf!
      use ISO_C_BINDING
      implicit none                           !InTf!
      character (len=*) :: string             !InTf!
      interface
        function f_chdir(path) result(status) BIND(C,name='chdir')
          import :: C_INT, C_CHAR
          character(C_CHAR), dimension(*), intent(IN) :: path
          integer(C_INT) :: status
        end function f_chdir
      end interface
!       integer,external :: f_RPN_COMM_chdir
!       RPN_COMM_chdir=f_RPN_COMM_chdir(trim(string)//achar(0))
      RPN_COMM_chdir=f_chdir(trim(string)//achar(0))
      return
      end function RPN_COMM_chdir             !InTf!
!
!       RPN_COMM_env_var : get value of environment variable
!       returns blanks if variable does not exist
!
!InTf!
      subroutine RPN_COMM_env_var(varname,value)  !InTf!
      implicit none                               !InTf!
      character (len=*), intent(IN) :: varname    !InTf!
      character (len=*), intent(OUT) :: value     !InTf!

      integer status,length
      value = " "
!       call RPN_COMM_getenv(trim(varname)//achar(0),value,len(value))
      CALL GET_ENVIRONMENT_VARIABLE(varname, value, length, status)
      if(status .ne. 0) value = " "
      return
      end subroutine RPN_COMM_env_var              !InTf!
!
!       get the 32 bit host identifier
!
!InTf!
      integer function RPN_COMM_hostid()       !InTf!
      use ISO_C_BINDING
      implicit none                            !InTf!
      interface
	function f_gethostid() result(id) BIND(C,name='gethostid')
	  import :: C_LONG
	  integer(C_LONG) :: id
	end function f_gethostid
      end interface
      RPN_COMM_hostid=f_gethostid()
      return
      end function RPN_COMM_hostid !InTf!
