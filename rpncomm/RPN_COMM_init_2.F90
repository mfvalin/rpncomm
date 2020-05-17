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
!InTf!
!!INTEGER FUNCTION RPN_COMM_init_multi_level_io(Userinit,Pelocal,Petotal,Pex,Pey,MultiGrids,Grids,ApplID,Io)
      INTEGER FUNCTION RPN_COMM_init_multi_level_io &                !InTf!
      (Userinit,Pelocal,Petotal,Pex,Pey,MultiGrids,Grids,ApplID,Io)  !InTf!
      use rpn_comm
      implicit none                                                  !InTf!
      external :: Userinit                                           !InTf!
      integer, intent(out)   :: Pelocal,Petotal                      !InTf!
      integer, intent(inout) :: Pex,Pey                              !InTf!
      integer, intent(in)    :: MultiGrids                           !InTf!
      integer, intent(in)    :: Grids                                !InTf!
      integer, intent(in)    :: ApplID                               !InTf!
      integer, intent(in)    :: Io                                   !InTf!
!arguments
!  I	Userinit	User routine that will be called by PE 0 to
!		get the processor grid topology if it is not supplied
!		(Pex .eq. 0 ) or (Pey .eq. 0)
!  O	Pelocal	PE rank (number) of local PE
!  O	Petotal	Number of PEs in job
!  I/O	Pex	Number of PEs along the X axis. If Pex=0 upon entry
!		it will be set to the proper value upon exit
!  I/O	Pey	Number of PEs along the Y axis. If Pey=0 upon entry
!		it will be set to the proper value upon exit
!  I    MultiGrids  number of simultaneous MultiGrids (usually 1)
!  I    Grids   number of grids in MultiGrids
!  I    ApplID  application ID, positive number
!  I/O  Io      Io mode (there will be mod(io,10) IO PEs per io/10 compute PEs)
!               Io = 182 : 2 IO PEs for 18 compute PEs (2 out of 20 PEs are IO PEs)
!
!notes
!	processor topology common /pe/ will be filled here
!	positions are calculated from 0 (ORIGIN 0)
!!
!       this is also intended to be a cleanup/refactoring of RPN_COMM_init_multi_level
!
      integer ierr, i, j, count, npe, reste, nslots, key, status
      logical mpi_started
      integer gr1, gr2
      INTEGER newgroup,rowcomm
      integer, allocatable, dimension(:) :: proc_indomm
      integer ndom, lndom, nproc, procmax,procmin
      type(domm), allocatable, dimension(:) :: locdom
      logical ok, allok
      logical RPN_COMM_grank
      integer RPN_COMM_petopo, RPN_COMM_version
      character *4096 SYSTEM_COMMAND
      character *32 access_mode
      integer ncolors,my_color,directory_file_num,iun
      integer version_marker, version_max, version_min
      integer pe_my_location(8)
      external RPN_COMM_unbind_process
      integer, external :: RPN_COMM_get_a_free_unit, RPN_COMM_set_timeout_alarm, fnom
!
!      if(RPM_COMM_IS_INITIALIZED) then ! ignore with warning message or abort ?
!      endif
      pe_indomm=-1
      pe_indomms=-1
      pe_dommtot=-1
      pe_a_domain=-1
      pe_me_a_domain=-1
      pe_all_domains=-1
      pe_me_all_domains=-1
      pe_multi_grid=-1
      pe_me_multi_grid=-1
      pe_grid=-1
      pe_me_grid=-1
      RPM_COMM_IS_INITIALIZED=.true.
      if( .not. WORLD_COMM_MPI_INIT ) then ! world not set before, use MPI_COMM_WORLD
           WORLD_COMM_MPI_INIT=.true.
           WORLD_COMM_MPI=MPI_COMM_WORLD
      endif
      RPN_COMM_init_multi_level_io = 0
      ok = .true.

      call MPI_INITIALIZED(mpi_started,ierr)
      status = RPN_COMM_set_timeout_alarm(60)
      if (.not. mpi_started ) call MPI_init(ierr)
      status = RPN_COMM_set_timeout_alarm(0)
      pe_wcomm=WORLD_COMM_MPI                     ! UNIVERSE at this point
      call MPI_COMM_RANK(pe_wcomm,pe_me,ierr)     ! rank in UNIVERSE
      call MPI_COMM_SIZE(pe_wcomm,pe_tot,ierr)    ! size of UNIVERSE

      if(pe_me == 0)then                   
        call RPN_COMM_env_var("RPN_COMM_MONITOR",SYSTEM_COMMAND)
        if(SYSTEM_COMMAND .ne. " ") then
          iun = RPN_COMM_get_a_free_unit()
          open(iun,file=trim(SYSTEM_COMMAND),form='FORMATTED')
          write(iun,'(A)')'mpi_init successful'
          close(iun)
        endif
      endif
      call RPN_COMM_unbind_process ! unbind processes if applicable (FULL_UNBIND environment variable, linux pnly)

!     check that all Processes use the same version of rpn_comm
      version_marker=RPN_COMM_version()
      call mpi_allreduce(version_marker, version_max, 1, MPI_INTEGER,&
     &                   MPI_BOR, WORLD_COMM_MPI, ierr)
      if(version_max .ne. version_marker)then
        write(rpn_u,*) 'ERROR: RPN_COMM version mismatch, ABORTING execution'
        call mpi_finalize(ierr)
        stop
      endif
!     all domains
      pe_all_domains = pe_wcomm
      pe_me_all_domains = pe_me
      pe_tot_all_domains = pe_tot
      call MPI_COMM_GROUP(pe_all_domains,pe_gr_all_domains,ierr)
!
      my_color = 0
      diag_mode=1
      SYSTEM_COMMAND=" "
      call RPN_COMM_env_var("RPN_COMM_DIAG",SYSTEM_COMMAND)
      if( SYSTEM_COMMAND .ne. " " ) read(SYSTEM_COMMAND,*) diag_mode
!
!     even if environment variable RPN_COMM_DOM is set, ignore it, deprecated feature
!     domain/appl split is now performed using applid
!
!     --------------------------------------------------------------------------
!     if appl ID is nonzero and positive, split by application
!
      if(applid > 0) then
        my_color = applid
        RPN_COMM_init_multi_level_io = my_color
        call MPI_COMM_SPLIT(WORLD_COMM_MPI,my_color,pe_me,pe_wcomm,ierr)  ! split using appl color
        call MPI_COMM_RANK(pe_wcomm,pe_me,ierr)               ! my rank
        call MPI_COMM_SIZE(pe_wcomm,pe_tot,ierr)              ! size of my appl
      endif
!     --------------------------------------------------------------------------
!     my domain (appl)
      my_colors(1) = my_color
      RPN_COMM_init_multi_level_io = my_color
      pe_a_domain = pe_wcomm
      pe_me_a_domain = pe_me
      pe_tot_a_domain = pe_tot
      if(pe_me_a_domain .eq.0 .and. diag_mode .ge. 1) then
        write(rpn_u,1000)my_color,pe_tot,MultiGrids,Grids
1000  	format('domain=',I1,I4,' PEs,',I3,' SuperGrids with',I3,' Grids each')
      endif
      call MPI_COMM_GROUP(pe_a_domain,pe_gr_a_domain,ierr)
!
!     verify that pe_tot_a_domain is a multiple of MultiGrids and Grids
!
      ok = .true.
      i=((pe_tot_a_domain/MultiGrids)/Grids)  ! number of PEs in a grid
      j=pe_tot_a_domain - i*MultiGrids*Grids  ! remainder
      if(j .ne. 0) ok = .false.               ! OOPS
      call mpi_allreduce(ok, allok, 1, MPI_LOGICAL,MPI_LAND,WORLD_COMM_MPI,ierr)
      if(.not.allok .and. pe_me_a_domain .eq.0) then
           write(rpn_u,*)'ERROR: number of PEs in domain is not a multiple of MultiGrids*Grids, ABORTING execution'
           call RPN_COMM_finalize(ierr)
           stop
      endif
!     --------------------------------------------------------------------------
!     if multigrids in program, (re)split domain into multigrids
!
      my_color = 0
      if(MultiGrids .gt. 1) then
        my_color=pe_me/(pe_tot/MultiGrids)
        RPN_COMM_init_multi_level_io = my_color
        call MPI_COMM_SPLIT(pe_a_domain,my_color,pe_me,pe_wcomm,ierr)
        call MPI_COMM_RANK(pe_wcomm,pe_me,ierr)               ! my rank
        call MPI_COMM_SIZE(pe_wcomm,pe_tot,ierr)              ! size of my subdomain
      endif
      my_colors(2) = my_color                                 ! my multigrid number
      pe_multi_grid = pe_wcomm                                ! multigrid communicator
      pe_me_multi_grid = pe_me                                ! my rank in multigrid
      pe_tot_multi_grid = pe_tot                              ! multigrid population
      call MPI_COMM_GROUP(pe_multi_grid,pe_gr_multi_grid,ierr)
      pe_indomms = pe_multi_grid                              ! multigrid communicator
      call MPI_COMM_GROUP(pe_indomms,pe_gr_indomms,ierr)
!
!     make communicator for PEs on same host in this multigrid
!       call RPN_COMM_split_by_node(pe_multi_grid,
!                                   pe_multigrid_host, pe_node_peers,
!                                   rank_on_host,      rank_in_peers, n_on_node, ierr)
!
!     --------------------------------------------------------------------------
!     if there are subgrids, (re)split multigrid into grids 
!     (compute and IO processors at this point)
!
      my_color = 0
      if(Grids .gt. 1) then
        Pelocal = pe_me                                       ! rank in my domain
        Petotal = pe_tot                                      ! number of pes in domain
        my_color=pe_me/(pe_tot/Grids)
        RPN_COMM_init_multi_level_io = my_color
        call MPI_COMM_SPLIT(pe_multi_grid,my_color,&
     &                      pe_me_multi_grid,pe_wcomm,ierr)
        call MPI_COMM_RANK(pe_wcomm,pe_me,ierr)               ! my rank
        call MPI_COMM_SIZE(pe_wcomm,pe_tot,ierr)              ! size of my subdomain
      endif
!     --------------------------------------------------------------------------
!     TODO : take care of IO processors  (at the grid level)
!            "grid" will have to be split into compute and IO processes
!            will need to look at RPN_COMM_IO_CONFIG environment variable
!            may have to redefine pe_wcomm to be compute PEs
!            do we instead want IO at the supergrid level ?
!            grids will have to cooperate at the supergrid level
!            my_colors for compute PEs only ?
!     options:
!          return to the caller with special code
!          call user supplied subroutine
!          call own IO server code
!     --------------------------------------------------------------------------
      my_colors(3) = my_color                                 ! my grid number in my multigrid
      pe_grid = pe_wcomm                                      ! grid communicator
      pe_me_grid = pe_me                                      ! my rank in grid
      pe_tot_grid = pe_tot                                    ! grid population
!      write(rpn_u,*)'pe_tot_grid=',pe_tot_grid
      call MPI_COMM_GROUP(pe_grid,pe_gr_grid,ierr)

      Pelocal = pe_me                                          ! my rank in my grid
      Petotal = pe_tot                                         ! number of pes in my grid

      pe_indomm = pe_grid
      pe_gr_indomm = pe_gr_grid
      pe_pe0 = 0
      pe_dommtot = pe_tot
      call MPI_COMM_GROUP(pe_wcomm,pe_gr_wcomm,ierr)
!     --------------------------------------------------------------------------
!     communicator for PEs on same host in this grid
      call mpi_comm_split_type(pe_grid, MPI_COMM_TYPE_SHARED, pe_me_grid, MPI_INFO_NULL, pe_grid_host, ierr)
!       my_color = abs(f_gethostid())  ! coloring by host identifier
!       call MPI_COMM_SPLIT(pe_grid,my_color,pe_me_grid,pe_grid_host,ierr)        ! same (sub)grid, same host node communicator
      call MPI_COMM_RANK(pe_grid_host,pe_me_grid_host,ierr)    ! my rank on this host
      call MPI_COMM_SIZE(pe_grid_host,pe_tot_grid_host,ierr)   ! population of this host
      call MPI_COMM_GROUP(pe_grid_host,pe_gr_grid_host,ierr)   ! group communicator
!     --------------------------------------------------------------------------
!     create peer to peer communicators between grids within a multigrid
!     do we want to include IO processors ? (probably not)
!
      my_color = pe_me_grid
      call MPI_COMM_SPLIT(pe_multi_grid,my_color,pe_me_multi_grid,pe_grid_peers,ierr)
      call MPI_COMM_SIZE(pe_grid_peers,pe_tot_peer,ierr)  ! This must be equal to the number of grids
      call MPI_COMM_RANK(pe_grid_peers,pe_me_peer,ierr)   ! in a multigrid (or else...)
      call MPI_COMM_GROUP(pe_grid_peers,pe_gr_grid_peers,ierr)
!     --------------------------------------------------------------------------
!     Grid initialization (compute PEs) 
!     get PEs along X and Y
!     call Userinit if appropriate
!
      if(pe_me .eq. pe_pe0) then
	  if ( Pex.eq.0 .or. Pey.eq.0  ) then ! get processor topology from Userinit
	      WORLD_pe(1)=pe_tot
	      WORLD_pe(2)=1
	      call Userinit(WORLD_pe(1),WORLD_pe(2))
	      if(WORLD_pe(1)*WORLD_pe(2).ne.pe_dommtot) then
		ok = .false.
		write(rpn_u,*) 'RPN_COMM_init_io: Inconsistency between'
		write(rpn_u,*) 'userinit Subroutine and total number'
		write(rpn_u,*) 'of PE: please check'
	      endif
	      if(diag_mode.ge.1) then
		write(rpn_u,*)'Requested topology = ',WORLD_pe(1),' by ',WORLD_pe(2)
		write(rpn_u,*)'Grid will use ',pe_dommtot,' processes'
	      endif
          else ! ( Pex.eq.0 .or. Pey.eq.0  )
	      write(rpn_u,*) 'RPN_COMM_init_io: Forced topology'
	      WORLD_pe(1) = Pex
	      WORLD_pe(2) = Pey
	      if(WORLD_pe(1)*WORLD_pe(2).ne.pe_dommtot) then
		ok = .false.
		write(rpn_u,*) 'RPN_COMM_init_io: Inconsistency between Pex'
		write(rpn_u,*) 'and Pey args and total number of PE: '
		write(rpn_u,*) 'please check'
	      endif
	      if(diag_mode.ge.1) then
		write(rpn_u,*)'Requested topology = ',WORLD_pe(1),' by ',WORLD_pe(2)
	      endif
	  endif ! ( Pex.eq.0 .or. Pey.eq.0  )
!
	  if(WORLD_pe(1)*WORLD_pe(2) .gt. pe_dommtot) then
	    write(rpn_u,*)' ERROR: not enough PEs for decomposition '
	    write(rpn_u,*)' REQUESTED=',WORLD_pe(1)*WORLD_pe(2),' AVAILABLE=',pe_dommtot
	    ok = .false.
	  endif
      endif  ! (pe_me .eq. pe_pe0)
!
      call mpi_allreduce(ok, allok, 1, MPI_LOGICAL,MPI_LAND,WORLD_COMM_MPI,ierr)
      if(.not.allok) then
           if(.not. ok) write(rpn_u,*)'ERROR DETECTED'
           call RPN_COMM_finalize(ierr)
           stop
      endif
!
!	send WORLD topology to all PEs. That will allow all PEs
!	to compute other PE topology parameters locally.
!       for doing this, we need to define some basic domains
!       communicators.

	call MPI_COMM_rank(pe_indomm,pe_medomm,ierr)
	pe_defcomm = pe_indomm
	pe_defgroup = pe_gr_indomm
!
!       broadcast number of PEs along X and Y axis
!       broadcast PE block sizes (deltai and deltaj)
!
        WORLD_pe(3) = deltai
        WORLD_pe(4) = deltaj
	call MPI_BCAST(WORLD_pe,size(WORLD_pe),MPI_INTEGER,0,&
     &	               pe_indomm,ierr)
	pe_nx  = WORLD_pe(1)
	pe_ny  = WORLD_pe(2)
	deltai = WORLD_pe(3)
	deltaj = WORLD_pe(4)
!
	if ( Pex.eq.0 .or. Pey.eq.0  ) then ! return processor topology
	  Pex = WORLD_pe(1)
	  Pey = WORLD_pe(2)
	endif
!
!	pe_pe0 is not equal to 0 if there are more than one domain
!	computational grid
!
	count = pe_pe0
!
!	fill tables containing the position along the X axis (pe_xtab)
!	and along the Y axis (pe_ytab) for all processors
!        write(rpn_u,*)'calling petopo'
	ierr = RPN_COMM_petopo(WORLD_pe(1),WORLD_pe(2))
!        write(rpn_u,*)'after petopo'

      pe_my_location(1) = pe_mex
      pe_my_location(2) = pe_mey
      pe_my_location(3) = pe_me_grid
      pe_my_location(4) = my_colors(3)
      pe_my_location(5) = pe_me_multi_grid
      pe_my_location(6) = my_colors(2)
      pe_my_location(7) = pe_me_a_domain
      pe_my_location(8) = my_colors(1)
      allocate(pe_location(8,0:pe_tot-1))
      call MPI_allgather(&
     &     pe_my_location,8,MPI_INTEGER,&
     &     pe_location,   8,MPI_INTEGER,WORLD_COMM_MPI,&
     &     ierr)
      if( pe_me_all_domains .eq. 0 .and. diag_mode .ge.3) then
        write(rpn_u,*)'                         FULL PE MAP'
        write(rpn_u,*)&
     & '    mex     mey   me(g)    grid  me(sg)   sgrid   me(d)  domain'
        do j=0,pe_tot_all_domains-1
           write(rpn_u,1001)(pe_location(i,j),i=1,8)
1001       format(8I8)
        enddo
      endif

      BLOC_EXIST   =.false.
      BLOC_SIZEX   = 1
      BLOC_SIZEY   = 1
      BLOC_mybloc  = 0
      BLOC_myblocx = 0
      BLOC_myblocy = 0
      BLOC_me      = pe_me
      BLOC_comm_world = pe_indomm
      BLOC_comm_row = pe_myrow
      BLOC_comm_col = pe_mycol
      BLOC_corner = pe_pe0
      BLOC_master = 0
      if(pe_me.eq.pe_pe0) then
         BLOC_master=1
      endif
      pe_bloc = pe_indomm

      call MPI_Group_incl(pe_gr_indomm, 1, 0, pe_gr_blocmaster, ierr) 
      call MPI_Comm_create(pe_indomm,pe_gr_blocmaster, &
     &            pe_blocmaster, ierr)

!       for each communicator, store the corresponding group
!
      call MPI_COMM_GROUP(pe_myrow,pe_gr_myrow,ierr)
      call MPI_COMM_GROUP(pe_mycol,pe_gr_mycol,ierr)
      call MPI_COMM_GROUP(pe_bloc,pe_gr_bloc,ierr)
!      write(rpn_u,*)'peer communicator =',pe_grid_peers
!      write(rpn_u,*)'peer grid size =',pe_tot_peer
!      write(rpn_u,*)'peer rank =',pe_me_peer
!
!      ! what follows has to be added to rpn_comm or another module
!
!      integer, pointer, dimension(:,:) :: grid_id_table
!      integer, dimension(3) :: id_table
!      integer :: pe_pe0s, pe_me_pe0s  ! communicator for PE0 set and rank within
!
!      ! end of addition to rpn_comm or another module
!
!      my_color = min(1,pe_me_grid)  ! i am a PE 0 or not
!      call MPI_COMM_SPLIT(pe_all_domains,my_color,&
!     &                    pe_me_all_domains,pe_pe0s,ierr)
!      call MPI_COMM_rank(pe_pe0s,pe_me_pe0s,ierr) ! communicator nicknamed RPN_COMM_PE0='PE_00'
!      if(pe_me_grid /= 0) then  ! make sure it is invalid if not a PE0
!        pe_pe0s = -1
!        pe_me_pe0 = -1
!      endif
!
!      allocate(grid_id_table(3,pe_tot_all_domains))
!      id_table(1) = my_colors(1)*1024*1024 
!      id_table(1) = id_table(1) + my_colors(2)*1024
!      id_table(1) = id_table(1) + my_colors(3)   ! this will become a "grid id"
!      id_table(2) = pe_me_grid   ! local rank in grid
!      id_table(3) = pe_me_all_domains  ! global rank in "universe"
!      call MPI_allgather(id_table,3,MPI_INTEGER,&
!     &                    grid_id_table,3,MPI_INTEGER,&
!     &                    pe_all_domains,ierr)   ! progagate grid_id/local_rank/global_rank table
      return
!      contains

      end FUNCTION RPN_COMM_init_multi_level_io                        !InTf!
