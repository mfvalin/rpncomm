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
        integer function RPN_COMM_datyp_indx(data_int)            !InTf!
        use rpn_comm
        implicit none                                             !InTf!
!        include 'mpif.h'
!        include 'rpn_comm.h'
        character(len=*), intent(IN) :: data_int                  !InTf!
        character(len=32) :: datatype
        integer :: i

        call rpn_comm_low2up(data_int,datatype)

        do i = 1 , size(type_tab)
          if(type_tab(i)%string == datatype) then
            RPN_COMM_datyp_indx = i
            return
          endif
        enddo

        RPN_COMM_datyp_indx = -1  ! error return
        return
        end function RPN_COMM_datyp_indx                           !InTf!
!InTf!
        integer function RPN_COMM_datyp(data_int)                 !InTf!
!   Luc Corbeil, 2000-11-20
!   lien entre datatype et MPI_datatype
        use rpn_comm
        implicit none                                             !InTf!
!        include 'mpif.h'
!        include 'rpn_comm.h'
        character(len=*), intent(IN) :: data_int                  !InTf!
        character(len=32) :: datatype
        integer :: i
        integer, external :: RPN_COMM_datyp_indx

        RPN_COMM_datyp = MPI_DATATYPE_NULL

        i = RPN_COMM_datyp_indx(data_int)
        if(i >= 0) RPN_COMM_datyp = type_tab(i)%number
        return

#if defined(DEPRECATED)
        call rpn_comm_low2up(data_int,datatype)

        RPN_COMM_datyp = -999999  ! precondition to error return

        do i = 1 , size(type_tab)
          if(type_tab(i)%string == datatype) then
            RPN_COMM_datyp = type_tab(i)%number
            return
          endif
        enddo
#endif
#if defined(DEPRECATED)
        goto 777
        if (datatype(1:11).eq.'MPI_INTEGER') then
           RPN_COMM_datyp=MPI_INTEGER
           return
        endif
        if (datatype(1:8).eq.'MPI_REAL') then
           RPN_COMM_datyp=MPI_REAL
           return
        endif
        if (datatype(1:13).eq.'MPI_CHARACTER') then
           RPN_COMM_datyp=MPI_CHARACTER
           return
        endif
        if (datatype(1:8).eq.'MPI_BYTE') then
           RPN_COMM_datyp=MPI_BYTE
           return
        endif
        if (datatype(1:12).eq.'MPI_INTEGER2') then
           RPN_COMM_datyp=MPI_INTEGER2
           return
        endif
        if (datatype(1:18).eq.'MPI_DOUBLE_COMPLEX') then
           RPN_COMM_datyp=MPI_DOUBLE_COMPLEX
           return
        endif
        if (datatype(1:20).eq.'MPI_DOUBLE_PRECISION') then 
           RPN_COMM_datyp=MPI_DOUBLE_PRECISION
           return
        endif
        if (datatype(1:9).eq.'MPI_REAL4') then
           RPN_COMM_datyp=MPI_REAL4
           return
        endif
        if (datatype(1:9).eq.'MPI_REAL8') then
           RPN_COMM_datyp=MPI_REAL8
           return
        endif
        if (datatype(1:11).eq.'MPI_COMPLEX') then
           RPN_COMM_datyp=MPI_COMPLEX
            return
        endif
        if (datatype(1:11).eq.'MPI_LOGICAL') then
           RPN_COMM_datyp=MPI_LOGICAL
           return
        endif
777     write(rpn_u,*) 'ERROR: Unknown datatype ',datatype
#endif
        return
        end function RPN_COMM_datyp                             !InTf!
!
!       fill an entity of type rpncomm_datatype from rpn_comm type string
!       dtyp_c : character version of data type
!       dtyp   : new item of type rpncomm_datatype
        subroutine RPN_COMM_i_datyp(dtyp_c,dtyp)           !InTf!
        use rpn_comm
!!      import :: rpncomm_datatype                           !InTf!
        implicit none
        type(rpncomm_datatype), intent(OUT) :: dtyp          !InTf!
        character(len=*), intent(IN) :: dtyp_c               !InTf!
        integer, external :: RPN_COMM_datyp_indx

        dtyp%p = C_LOC(WORLD_COMM_MPI)            ! signature
        dtyp%t1 = RPN_COMM_datyp_indx(dtyp_c)     ! index of datatype from internal table
        dtyp%t2 = type_tab(dtyp%t1)%number        ! datatype value
        end subroutine RPN_COMM_i_datyp                    !InTf!
!
!       is dtyp a valid item of type rpncomm_datatype ?
!
        function RPN_COMM_i_valid_datyp(dtyp) result(valid)    !InTf!
        use rpn_comm
!!      import :: rpncomm_datatype                           !InTf!
        implicit none
        type(rpncomm_datatype), intent(IN) :: dtyp           !InTf!
        logical :: valid                                     !InTf!
        integer, external :: RPN_COMM_datyp_indx
        type(C_PTR) :: temp

        temp = C_LOC(WORLD_COMM_MPI)                      ! correct signature pointer ?
        valid = c_associated( dtyp%p , temp )
        if(.not. valid) return

        valid = (dtyp%t1 > 0) .and. (dtyp%t1 <= size(type_tab))   ! plausible t1 ?
        if(.not. valid) return

        valid = type_tab(dtyp%t1)%number == dtyp%t2   ! consistent t1 and t2 ?

        end function RPN_COMM_i_valid_datyp                    !InTf!
