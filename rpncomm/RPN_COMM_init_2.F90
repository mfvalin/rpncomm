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
!
! multiple levels initialization function application/supergrid/grid/block
!
!!INTEGER FUNCTION RPN_COMM_init_all_levels(Userinit,Pelocal,Petotal,Pex,Pey,MultiGrids,Grids,AppID,Io)
      INTEGER FUNCTION RPN_COMM_init_all_levels &                    !InTf!
      (Userinit,Pelocal,Petotal,Pex,Pey,MultiGrids,Grids,AppID,Io)   !InTf!
      use rpn_comm
      use RPN_COMM_mpi_layout
      implicit none                                                  !InTf!
      external :: Userinit                                           !InTf!
      integer, intent(out)   :: Pelocal,Petotal                      !InTf!
      integer, intent(inout) :: Pex,Pey                              !InTf!
      integer, intent(in)    :: MultiGrids                           !InTf!
      integer, intent(in)    :: Grids                                !InTf!
      character(len=*)       :: AppID                                !InTf!
      integer, intent(in)    :: Io                                   !InTf!
!arguments
!  I	Userinit	User routine that will be called by PE 0 to
!		get the processor grid topology if it is not supplied
!		(Pex .eq. 0 ) or (Pey .eq. 0)
!  O	Pelocal	PE rank (number) of local PE in its GRID
!  O	Petotal	Number of PEs in this PEs GRID
!  I/O	Pex	Number of PEs along the X axis. If Pex=0 upon entry
!		it will be set to the proper value upon exit
!  I/O	Pey	Number of PEs along the Y axis. If Pey=0 upon entry
!		it will be set to the proper value upon exit
!  I    MultiGrids  number of simultaneous MultiGrids (usually 1)
!  I    Grids   number of grids in MultiGrids
!  I    AppID   application ID, character string, 5 chars max
!  I/O  Io      Io mode (there will be mod(io,10) service (IO) PEs per io/10 compute PEs)
!               Io = 182 : 2 service PEs for 18 compute PEs (2 out of 20 PEs are IO PEs)
!
!notes
!	processor topology common /pe/ will be filled here
!	positions are calculated from 0 (ORIGIN 0)
!!
!       this is also intended to be a cleanup/refactoring of RPN_COMM_init_multi_level
!
      interface
	function gethostid() result(id) BIND(C,name='gethostid')
	  import :: C_LONG
	  integer(C_LONG) :: id
	end function gethostid
      end interface
      integer ierr, i, j, count, npe, reste, nslots, key, status, core
      integer :: pe_type    ! 0 = compute, 1 = service
      integer :: service
      integer :: compute
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
      integer pe_my_location(10)
      external RPN_COMM_unbind_process
      integer, external :: RPN_COMM_get_a_free_unit, RPN_COMM_set_timeout_alarm
      integer :: ApplID
      character(len=5) :: appid5
!
!     build application (domain) "color" from first 5 (at most) characters of AppID
      ApplID = 0
      do i = 1, min(5 , len(AppID))    ! 6 bits per character ASCII 32-96, case insensitive
        ApplID = ApplID * 64 + and(63, 32 + ichar(AppID(i:i)))  ! i th character
      enddo
      call RPN_COMM_reset_mpi_layout   ! initialize NEW style layout structure
      call RPN_COMM_get_core_and_numa(core, ml%numa)  ! get numa space for this PE
      ml%host = gethostid()
      compute = Io / 10         ! number of compute PEs in a PE block (compute + service PEs)
      service = mod(Io, 10)     ! number of service(IO) PEs in a PE block
!     initialize OLD style variables
      pe_indomm=-1
      pe_indomms=-1
      pe_dommtot=-1
      pe_all_domains=MPI_COMM_NULL
      pe_me_all_domains=-1
      pe_a_domain=MPI_COMM_NULL
      pe_me_a_domain=-1
      pe_multi_grid=MPI_COMM_NULL
      pe_me_multi_grid=-1
      pe_grid=MPI_COMM_NULL
      pe_me_grid=-1
!
      RPM_COMM_IS_INITIALIZED=.true.
      if( .not. WORLD_COMM_MPI_INIT ) then ! world not set before, use MPI_COMM_WORLD
           WORLD_COMM_MPI_INIT=.true.
           WORLD_COMM_MPI=MPI_COMM_WORLD
      endif
!
      call MPI_INITIALIZED(mpi_started,ierr)      ! is the MPI library already initialized ?
      status = RPN_COMM_set_timeout_alarm(60)     ! maximum of 60 seconds for MPI_init
      if (.not. mpi_started ) call MPI_init(ierr)
      status = RPN_COMM_set_timeout_alarm(0)      ! timeout reset to infinity (no timeout)
!
!     --------------------------------------------------------------------------
!     all applications (domains)
      pe_wcomm=WORLD_COMM_MPI                     ! the UNIVERSE
      call MPI_COMM_RANK(pe_wcomm,pe_me,ierr)     ! rank in UNIVERSE
      call MPI_COMM_SIZE(pe_wcomm,pe_tot,ierr)    ! size of UNIVERSE
!     NEW style
      ml%comm%wrld%all = pe_wcomm
      ml%rank%wrld%all = pe_me
      ml%size%wrld%all = pe_tot
      call MPI_COMM_SPLIT_TYPE(pe_wcomm, MPI_COMM_TYPE_SHARED, pe_me, MPI_INFO_NULL, ml%comm%wrld%same_node, ierr) ! same node
      call MPI_COMM_RANK(ml%comm%wrld%same_node, ml%rank%wrld%same_node, ierr)
      call MPI_COMM_SIZE(ml%comm%wrld%same_node, ml%size%wrld%same_node, ierr)
      call MPI_COMM_SPLIT(ml%comm%wrld%same_node, ml%numa, pe_me, ml%comm%wrld%same_numa, ierr)                    ! same numa space
      call MPI_COMM_RANK(ml%comm%wrld%same_numa, ml%rank%wrld%same_numa, ierr)
      call MPI_COMM_SIZE(ml%comm%wrld%same_numa, ml%size%wrld%same_numa, ierr)
!     OLD style
      pe_all_domains = pe_wcomm
      pe_me_all_domains = pe_me
      pe_tot_all_domains = pe_tot
      call MPI_COMM_GROUP(pe_all_domains,pe_gr_all_domains,ierr)
!
!     --------------------------------------------------------------------------
      if(pe_me == 0)then                           ! if requested produce a "status" file
        call RPN_COMM_env_var("RPN_COMM_MONITOR",SYSTEM_COMMAND)  ! export RPN_COMM_MONITOR=filename
        if(SYSTEM_COMMAND .ne. " ") then
          iun = RPN_COMM_get_a_free_unit()
          open(iun,file=trim(SYSTEM_COMMAND),form='FORMATTED')
          write(iun,'(A)')'mpi_init successful'
          close(iun)
        endif
      endif
      call RPN_COMM_unbind_process ! unbind processes if requested (FULL_UNBIND environment variable, linux only)
!
!     --------------------------------------------------------------------------
!     set diagnostic mode. if >= 2, some diagnostics are printed (pe topo)
      diag_mode=1
      SYSTEM_COMMAND=" "
      call RPN_COMM_env_var("RPN_COMM_DIAG",SYSTEM_COMMAND)  ! export RPN_COMM_DIAG=n
      if( SYSTEM_COMMAND .ne. " " ) read(SYSTEM_COMMAND,*) diag_mode
!
!     --------------------------------------------------------------------------
!     check that all Processes use the same version of rpn_comm
      version_marker=RPN_COMM_version()
      call mpi_allreduce(version_marker, version_max, 1, MPI_INTEGER, MPI_BOR, WORLD_COMM_MPI, ierr)
      if(version_max .ne. version_marker)then
        write(rpn_u,*) 'ERROR: RPN_COMM version mismatch, PLS recompile, ABORTING execution'
        call mpi_finalize(ierr)
        stop
      endif
!
!     even if environment variable RPN_COMM_DOM is set, ignore it (deprecated feature)
!
!     --------------------------------------------------------------------------
!     split WORLD_COMM_MPI by application(domain) using applid as "color"
!
      my_color = ApplID
      call MPI_COMM_SPLIT(WORLD_COMM_MPI, my_color, pe_me, pe_wcomm, ierr)
      call MPI_COMM_RANK(pe_wcomm,pe_me,ierr)               ! my rank
      call MPI_COMM_SIZE(pe_wcomm,pe_tot,ierr)              ! size of my appl
      if(pe_me .eq. 0 .and. diag_mode .ge. 1) then
        write(rpn_u,1000)trim(AppID),pe_tot,MultiGrids,Grids
1000    format('domain=',A,I6,' PEs,',I3,' SuperGrids with',I3,' Grids each')
      endif
      ok = .true.                             
      i=((pe_tot/MultiGrids)/Grids)           ! number of PEs in a grid
      j=pe_tot - i*MultiGrids*Grids           ! remainder
      ok = (j .eq. 0)                         ! check that pe_tot is a multiple of MultiGrids and Grids
      if(.not. ok .and. pe_me .eq. 0) &
         write(rpn_u,*)'ERROR: number of PEs in domain is not a multiple of MultiGrids*Grids, ABORTING execution'
      call mpi_allreduce(ok, allok, 1, MPI_LOGICAL,MPI_LAND,WORLD_COMM_MPI,ierr)
      if(.not.allok .and. pe_me .eq. 0) then
           call RPN_COMM_finalize(ierr)
           stop
      endif
      RPN_COMM_init_all_levels = my_color
!     NEW style
      ml%comm%appl%all = pe_wcomm
      ml%rank%appl%all = pe_me
      ml%size%appl%all = pe_tot
      call MPI_COMM_SPLIT_TYPE(pe_wcomm, MPI_COMM_TYPE_SHARED, pe_me, MPI_INFO_NULL, ml%comm%appl%same_node, ierr) ! same node
      call MPI_COMM_RANK(ml%comm%appl%same_node, ml%rank%appl%same_node, ierr)
      call MPI_COMM_SIZE(ml%comm%appl%same_node, ml%size%appl%same_node, ierr)
      call MPI_COMM_SPLIT(ml%comm%appl%same_node, ml%numa, pe_me, ml%comm%appl%same_numa, ierr)                    ! same numa space
      call MPI_COMM_RANK(ml%comm%appl%same_numa, ml%rank%appl%same_numa, ierr)
      call MPI_COMM_SIZE(ml%comm%appl%same_numa, ml%size%appl%same_numa, ierr)
      ml%colors(1)     = my_color
!     OLD style
      my_colors(1) = my_color
      pe_a_domain = pe_wcomm
      pe_me_a_domain = pe_me
      pe_tot_a_domain = pe_tot
      call MPI_COMM_GROUP(pe_a_domain,pe_gr_a_domain,ierr)
if(pe_me == 0) print *,'application split done'
!
!     --------------------------------------------------------------------------
!     if more than 1 supergrid in application, split application into supergrids
!     (compute and service PEs at this point)
!
      my_color = 0
      pe_wcomm = ml%comm%appl%all
      if(MultiGrids .gt. 1) then
        my_color=ml%rank%appl%all / (ml%size%appl%all / MultiGrids)
        RPN_COMM_init_all_levels = my_color
        call MPI_COMM_SPLIT(ml%comm%appl%all, my_color, ml%rank%appl%all, pe_wcomm, ierr)
      endif
      call MPI_COMM_RANK(pe_wcomm,pe_me,ierr)                 ! rank in supergrid
      call MPI_COMM_SIZE(pe_wcomm,pe_tot,ierr)                ! size of supergrid
!     NEW style
      ml%comm%sgrd%all = pe_wcomm                             ! supergrid communicator
      ml%rank%sgrd%all = pe_me                                ! ordinal in supergrid communicator
      ml%size%sgrd%all = pe_tot                               ! population in supergrid communicator
      ml%colors(2)     = my_color
      ml%comm%sgrd%row = MPI_COMM_NULL                        ! row and column are not defined for supergrids
      ml%rank%sgrd%row = -1
      ml%size%sgrd%row = -1
      ml%comm%sgrd%column = MPI_COMM_NULL
      ml%rank%sgrd%column = -1
      ml%size%sgrd%column = -1
!
!     supergrid to supergrid peers in application (PEs with same rank in supergrid)
!
      call MPI_COMM_SPLIT(ml%comm%appl%all, pe_me, ml%rank%appl%all, ml%comm%sgrd%grid_peer, ierr)  ! supergrid peers
      call MPI_COMM_RANK(ml%comm%sgrd%grid_peer, ml%rank%sgrd%grid_peer, ierr)
      call MPI_COMM_SIZE(ml%comm%sgrd%grid_peer, ml%size%sgrd%grid_peer, ierr)
!
!     split supergrid communicator into same node and node peers
!
      call MPI_COMM_SPLIT_TYPE(pe_wcomm, MPI_COMM_TYPE_SHARED, pe_me, MPI_INFO_NULL, ml%comm%sgrd%same_node, ierr)  ! same node
      call MPI_COMM_RANK(ml%comm%sgrd%same_node, ml%rank%sgrd%same_node, ierr)
      call MPI_COMM_SIZE(ml%comm%sgrd%same_node, ml%size%sgrd%same_node, ierr)
      call MPI_COMM_SPLIT(pe_wcomm, ml%rank%sgrd%same_node, ml%rank%sgrd%all, ml%comm%sgrd%node_peer, ierr)         ! node peers
      call MPI_COMM_RANK(ml%comm%sgrd%node_peer, ml%rank%sgrd%node_peer, ierr)
      call MPI_COMM_SIZE(ml%comm%sgrd%node_peer, ml%size%sgrd%node_peer, ierr)
!
!     split same node into same NUMA space and supergrid into NUMA space peers
!
      call MPI_COMM_SPLIT(ml%comm%sgrd%same_node, ml%numa, pe_me, ml%comm%sgrd%same_numa, ierr)             ! same numa space
      call MPI_COMM_RANK(ml%comm%sgrd%same_numa, ml%rank%sgrd%same_numa, ierr)
      call MPI_COMM_SIZE(ml%comm%sgrd%same_numa, ml%size%sgrd%same_numa, ierr)
      call MPI_COMM_SPLIT(pe_wcomm, ml%rank%sgrd%same_numa, ml%rank%sgrd%all, ml%comm%sgrd%numa_peer, ierr) ! numa peers
      call MPI_COMM_RANK(ml%comm%sgrd%numa_peer, ml%rank%sgrd%numa_peer, ierr)
      call MPI_COMM_SIZE(ml%comm%sgrd%numa_peer, ml%size%sgrd%numa_peer, ierr)
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
!
!     split supergrid into compute and service (IO) PEs
!
      pe_type = 0   ! compute by default, service PEs out of every (service+compute) are service PEs
      if( mod(pe_me,compute+service) < service ) pe_type = 1    ! determine if compute or service PE
! print *,'PE type =',pe_type
!
      call MPI_COMM_SPLIT(ml%comm%sgrd%all, pe_type, pe_me, pe_wcomm, ierr)
      call MPI_COMM_RANK(pe_wcomm,pe_me,ierr)                   ! rank
      call MPI_COMM_SIZE(pe_wcomm,pe_tot,ierr)                  ! size
      if(pe_type == 0) then                                     ! compute PEs
!     	NEW style
	ml%comm%sgrd%compute = pe_wcomm
	ml%rank%sgrd%compute = pe_me
	ml%size%sgrd%compute = pe_tot
	ml%comm%sgrd%service = MPI_COMM_NULL
	ml%rank%sgrd%service = -1
	ml%size%sgrd%service = -1
      else                                                      ! service PEs
	ml%comm%sgrd%compute = MPI_COMM_NULL
	ml%rank%sgrd%compute = -1
	ml%size%sgrd%compute = -1
	ml%comm%sgrd%service = pe_wcomm
	ml%rank%sgrd%service = pe_me
	ml%size%sgrd%service = pe_tot
      endif
!     OLD style, compute PEs only
      my_colors(2) = my_color                                 ! my multigrid number
      pe_multi_grid = pe_wcomm                                ! multigrid communicator
      pe_me_multi_grid = pe_me                                ! my rank in multigrid
      pe_tot_multi_grid = pe_tot                              ! multigrid population
      call MPI_COMM_GROUP(pe_multi_grid,pe_gr_multi_grid,ierr)
      pe_indomms = pe_multi_grid                              ! multigrid communicator
      call MPI_COMM_GROUP(pe_indomms,pe_gr_indomms,ierr)
! print *,'supergrid split done'
!     --------------------------------------------------------------------------
!     if more than 1 grid per supergrid, split supergrid into grids 
!     (compute and service PEs at this point)
!
      my_color = 0
      pe_wcomm = ml%comm%sgrd%all
      if(Grids .gt. 1) then
        my_color=ml%rank%sgrd%all / (ml%size%sgrd%all / Grids)
        RPN_COMM_init_all_levels = my_color
        call MPI_COMM_SPLIT(ml%comm%sgrd%all, my_color, ml%rank%sgrd%all, pe_wcomm, ierr)
      endif
      call MPI_COMM_RANK(pe_wcomm,pe_me,ierr)               ! my rank
      call MPI_COMM_SIZE(pe_wcomm,pe_tot,ierr)              ! size of my subdomain
!     NEW style
      ml%comm%grid%all = pe_wcomm
      ml%rank%grid%all = pe_me
      ml%size%grid%all = pe_tot
      ml%colors(3)     = my_color
!
!     grid to grid peers in supergrid (PEs with same rank in grid)
!
! print *,'splitting for grid peers'
      call MPI_COMM_SPLIT(ml%comm%sgrd%all, pe_me, ml%rank%sgrd%all, ml%comm%grid%grid_peer, ierr)  ! grid peers
      call MPI_COMM_RANK(ml%comm%grid%grid_peer, ml%rank%grid%grid_peer, ierr)
      call MPI_COMM_SIZE(ml%comm%grid%grid_peer, ml%size%grid%grid_peer, ierr)
!     OLD style
      pe_grid_peers = ml%comm%grid%grid_peer
      call MPI_COMM_SIZE(pe_grid_peers,pe_tot_peer,ierr)  ! This must be equal to the number of grids
      call MPI_COMM_RANK(pe_grid_peers,pe_me_peer,ierr)   ! in a supergrid (or else...)
      call MPI_COMM_GROUP(pe_grid_peers,pe_gr_grid_peers,ierr)
!
!     split grid into nodes and node peers
!
! print *,'splitting for same node'
      call MPI_COMM_SPLIT_TYPE(pe_wcomm, MPI_COMM_TYPE_SHARED, pe_me, MPI_INFO_NULL, ml%comm%grid%same_node, ierr)  ! same node
      call MPI_COMM_RANK(ml%comm%grid%same_node, ml%rank%grid%same_node, ierr)
      call MPI_COMM_SIZE(ml%comm%grid%same_node, ml%size%grid%same_node, ierr)
! print *,'splitting for node peers'
      call MPI_COMM_SPLIT(pe_wcomm, ml%rank%grid%same_node, ml%rank%sgrd%all, ml%comm%grid%node_peer, ierr) ! node peers
      call MPI_COMM_RANK(ml%comm%grid%node_peer, ml%rank%grid%node_peer, ierr)
      call MPI_COMM_SIZE(ml%comm%grid%node_peer, ml%size%grid%node_peer, ierr)
!     OLD style
      pe_grid_host = ml%comm%grid%same_node
      call MPI_COMM_RANK(pe_grid_host,pe_me_grid_host,ierr)    ! my rank on this host
      call MPI_COMM_SIZE(pe_grid_host,pe_tot_grid_host,ierr)   ! population of this host
      call MPI_COMM_GROUP(pe_grid_host,pe_gr_grid_host,ierr)   ! group communicator
!
!     split grid into NUMA spaces and NUMA peers
!
! print *,'splitting for same numa'
      call MPI_COMM_SPLIT(ml%comm%grid%same_node, ml%numa, pe_me, ml%comm%grid%same_numa, ierr)                     ! same numa
      call MPI_COMM_RANK(ml%comm%grid%same_numa, ml%rank%grid%same_numa, ierr)
      call MPI_COMM_SIZE(ml%comm%grid%same_numa, ml%size%grid%same_numa, ierr)
! print *,'splitting for numa peers'
      call MPI_COMM_SPLIT(ml%comm%grid%all, ml%rank%grid%same_numa, ml%rank%grid%all, ml%comm%grid%numa_peer, ierr) ! numa peers
      call MPI_COMM_RANK(ml%comm%grid%numa_peer, ml%rank%grid%numa_peer, ierr)
      call MPI_COMM_SIZE(ml%comm%grid%numa_peer, ml%size%grid%numa_peer, ierr)
!
!     split grid into compute and service (IO) PEs
!
!     pe_type already set above when dealing with supergrids
! print *,'splitting for service'
      call MPI_COMM_SPLIT(ml%comm%grid%all, pe_type, pe_me, pe_wcomm, ierr)
      call MPI_COMM_RANK(pe_wcomm,pe_me,ierr)                   ! rank
      call MPI_COMM_SIZE(pe_wcomm,pe_tot,ierr)                  ! size
      if(pe_type == 0) then                                     ! compute PEs
!     	NEW style
	ml%comm%grid%compute = pe_wcomm
	ml%rank%grid%compute = pe_me
	ml%size%grid%compute = pe_tot
	ml%comm%grid%service = MPI_COMM_NULL
	ml%rank%grid%service = -1
	ml%size%grid%service = -1
      else                                                      ! service PEs
	ml%comm%grid%compute = MPI_COMM_NULL
	ml%rank%grid%compute = -1
	ml%size%grid%compute = -1
	ml%comm%grid%service = pe_wcomm
	ml%rank%grid%service = pe_me
	ml%size%grid%service = pe_tot
      endif
!     OLD style, compute PEs only
      my_colors(3) = my_color                                 ! my grid number in my multigrid
      pe_grid = pe_wcomm                                      ! grid communicator
      pe_me_grid = pe_me                                      ! my rank in grid
      pe_tot_grid = pe_tot                                    ! grid population
      call MPI_COMM_GROUP(pe_grid,pe_gr_grid,ierr)

      Pelocal = pe_me                                          ! my rank in my grid
      Petotal = pe_tot                                         ! number of pes in my grid

      pe_indomm = pe_grid
      pe_gr_indomm = pe_gr_grid
      pe_pe0 = 0
      pe_dommtot = pe_tot
      call MPI_COMM_GROUP(pe_wcomm,pe_gr_wcomm,ierr)
! if(pe_me == 0) print *,'grid split done'
!     --------------------------------------------------------------------------
!     Grid initialization (compute PEs) 
!     get PEs along X and Y
!     call Userinit if appropriate
!
      ok = .true.
      if(pe_me .eq. pe_pe0) then
	  if ( Pex.eq.0 .or. Pey.eq.0  ) then ! get processor topology from Userinit
	      WORLD_pe(1)=pe_tot
	      WORLD_pe(2)=1
	      call Userinit(WORLD_pe(1),WORLD_pe(2))
	      if(WORLD_pe(1)*WORLD_pe(2).ne.pe_dommtot) then
		ok = .false.
		write(rpn_u,*) 'ERROR: (RPN_COMM_init_io) Inconsistency between'
		write(rpn_u,*) '       userinit Subroutine and total number of PEs'
		write(rpn_u,*) '       please double check topology'
	      endif
	      if(diag_mode.ge.1) then
		write(rpn_u,*)'INFO: Requested topology = ',WORLD_pe(1),' by ',WORLD_pe(2)
		write(rpn_u,*)'      Grid will use ',pe_dommtot,' processes'
	      endif
          else ! ( Pex.ne.0 .and. Pey.ne.0  )
	      write(rpn_u,*) 'INFO: (RPN_COMM_init_io) Forced topology :',Pex,' by',Pey
	      WORLD_pe(1) = Pex
	      WORLD_pe(2) = Pey
	      if(WORLD_pe(1)*WORLD_pe(2).ne.pe_dommtot) then
		ok = .false.
		write(rpn_u,*) 'ERROR: (RPN_COMM_init_io) Inconsistency between Pex, Pey and total number of PEs'
		write(rpn_u,*) '       please double check topology'
	      endif
	      if(diag_mode.ge.1) then
		write(rpn_u,*)'Requested topology =',WORLD_pe(1),' by ',WORLD_pe(2)
	      endif
	  endif ! ( Pex.eq.0 .or. Pey.eq.0  )
!
	  if(WORLD_pe(1)*WORLD_pe(2) .gt. pe_dommtot) then
	    write(rpn_u,*)' ERROR: not enough PEs for requested decomposition '
	    write(rpn_u,*)'        REQUESTED=',WORLD_pe(1)*WORLD_pe(2)
	    write(rpn_u,*)'        AVAILABLE=',pe_dommtot
	    ok = .false.
	  endif
      endif  ! (pe_me .eq. pe_pe0)
!
      call mpi_allreduce(ok, allok, 1, MPI_LOGICAL, MPI_LAND, WORLD_COMM_MPI, ierr)
      if(.not.allok) then
           if(.not. ok .and. pe_me .eq. pe_pe0 ) write(rpn_u,*)'ERROR: problem in grid initialization'
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
      call MPI_BCAST(WORLD_pe, size(WORLD_pe), MPI_INTEGER, 0, pe_indomm, ierr)
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
!     --------------------------------------------------------------------------
!     PE topology
      ierr = RPN_COMM_petopo(WORLD_pe(1),WORLD_pe(2))
      ml%comm%grid%row    = pe_myrow
      call MPI_COMM_RANK(ml%comm%grid%compute, ml%rank%grid%row, ierr)
      call MPI_COMM_SIZE(ml%comm%grid%compute, ml%size%grid%row, ierr)
      ml%comm%grid%column = pe_mycol
      call MPI_COMM_RANK(ml%comm%grid%compute, ml%rank%grid%column, ierr)
      call MPI_COMM_SIZE(ml%comm%grid%compute, ml%size%grid%column, ierr)

      pe_my_location(1) = pe_mex
      pe_my_location(2) = pe_mey
      pe_my_location(3) = pe_me_grid
      pe_my_location(4) = my_colors(3)
      pe_my_location(5) = pe_me_multi_grid
      pe_my_location(6) = my_colors(2)
      pe_my_location(7) = pe_me_a_domain
      pe_my_location(8) = my_colors(1)
      pe_my_location(9) = ml%host
      pe_my_location(10)= ml%numa

      pe_tot = ml%size%wrld%all
print *,'pe_tot =',pe_tot
! diag_mode = 3
      allocate(pe_location(10,0:pe_tot-1))
      call MPI_allgather( &
     &     pe_my_location,10,MPI_INTEGER, &
     &     pe_location,   10,MPI_INTEGER, &
     &     WORLD_COMM_MPI, ierr)
      if( pe_me_all_domains .eq. 0 .and. diag_mode .ge.3) then
        write(rpn_u,*)'                         FULL PE MAP'
        write(rpn_u,*)'    mex     mey   me(g)    grid  me(sg)   sgrid   me(d)  domain host   numa'
        do j=0,pe_tot_all_domains-1
           appid5 = '     '
           reste = pe_location(8,j)
           do i = 5, 1, -1
             appid5(i:i) = char(32 + mod(reste,64))
             reste = reste / 64
           enddo
           write(rpn_u,1001)pe_location(1:8,j),appid5,pe_location(9:10,j)
1001       format(7I8,I12,3X,A5,Z9,I3)
        enddo
      endif
      deallocate(pe_location)
!     --------------------------------------------------------------------------
!     PE blocks, initialized to 1 x 1 blocks
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
      call MPI_Comm_create(pe_indomm,pe_gr_blocmaster, pe_blocmaster, ierr)

!       for each communicator, store the corresponding group
!
      call MPI_COMM_GROUP(pe_myrow,pe_gr_myrow,ierr)
      call MPI_COMM_GROUP(pe_mycol,pe_gr_mycol,ierr)
      call MPI_COMM_GROUP(pe_bloc,pe_gr_bloc,ierr)
      return
      end FUNCTION RPN_COMM_init_all_levels                        !InTf!
