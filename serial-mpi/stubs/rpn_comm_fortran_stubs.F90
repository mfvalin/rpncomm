      subroutine MPI_abort
      write(6,*) 'MPI_abort is called, ABORT'
      call ABORT
      return
      end
!
      subroutine MPI_allgather (a, cnt, b, c, cnt2, d, e, ierr )
        implicit none
        include 'mpif.h'
        integer MPI_STUBS_length
        integer a(*),cnt,b,c(*),cnt2,d,e,ierr,i
        do i=1,min(cnt,cnt2)*MPI_STUBS_length(b)
           c(i)=a(i)
        enddo
        ierr=MPI_SUCCESS
        return
      end
!
      subroutine MPI_allgatherv
      write(6,*) 'MPI_allgatherv not authorized in serial-mpi mode, ABORT'
      call ABORT
      return
      end
!
      subroutine MPI_allreduce(send,recv,ni,l,m,n,ierr)
      implicit none
!
      include 'mpif.h'
!
      integer MPI_STUBS_length
      integer i,l,m,n,ierr
      integer send(*),recv(*),ni
      do i=1,ni*MPI_STUBS_length(l)
        recv(i)=send(i)
      enddo
      ierr=MPI_SUCCESS
      return
      end
!
      subroutine MPI_alltoall
      write(6,*)'MPI_alltoall not authorized in serial-mpi mode, ABORT'
      call ABORT
      return
      end
!
      subroutine MPI_alltoallv
      write(6,*)'MPI_alltoallv not authorized in serial-mpi mode, ABORT'
      call ABORT
      return
      end
!
      subroutine MPI_barrier(i,ierr)
      implicit none
      include 'mpif.h'
      integer i,ierr
      ierr=MPI_SUCCESS
      return
      end
!
      subroutine MPI_bcast(i,j,k,l,m,ierr)
      implicit none
!
      include 'mpif.h'
!
      integer i,j,k,l,m,ierr
      ierr=MPI_SUCCESS
      return
      end
!
      subroutine mpi_comm_create(pe_indomm,pe_gr_blocmaster, pe_blocmaster, ierr)
      implicit none
      include 'mpif.h'
      integer pe_indomm,pe_gr_blocmaster,  pe_blocmaster, ierr
      pe_blocmaster = 0
      ierr = MPI_SUCCESS
      return
      end
!
      subroutine MPI_comm_group(i,group,ierr)
      implicit none
      integer i,group,ierr
!
      include 'mpif.h'
!
      group=-1
      ierr=MPI_SUCCESS
      return
      end
!
      subroutine MPI_COMM_RANK(comm,pe_me,ierr)
      implicit none
!
      include 'mpif.h'
!
      integer comm,pe_me,ierr
      ierr=MPI_SUCCESS
      pe_me=0
      return
      end
!
      subroutine MPI_COMM_SIZE(comm,pe_tot,ierr)
      implicit none
!
      include 'mpif.h'
!
      integer comm,pe_tot,ierr
      ierr=MPI_SUCCESS
      pe_tot=1
      return
      end
!
      subroutine MPI_comm_split(i,j,k,newcomm,ierr)
      implicit none
      integer i,j,k,newcomm,ierr
!
      include 'mpif.h'
!
      newcomm=-1
      ierr=MPI_SUCCESS
      return
      end
!
      subroutine MPI_finalize
      return
      end
!
      subroutine MPI_gather(a, cnt, b, c, cnt2, d, e, f,ierr )
      implicit none
      include 'mpif.h'
      integer MPI_STUBS_length
      integer a(*),cnt,b,c(*),cnt2,d,e,f,ierr,i
      do i=1,min(cnt,cnt2)*MPI_STUBS_length(b)
         c(i)=a(i)
      enddo
      ierr=MPI_SUCCESS
      return
      end
!
      subroutine MPI_gatherv(a, cnt, b, c, cnt2s, displ, d, e, f,ierr )
      implicit none
      include 'mpif.h'
      integer MPI_STUBS_length
      integer a(*),cnt,b,c(*),cnt2s(*),cnt2,d,e,f,ierr,i,displ(*)
        cnt2=cnt2s(1)
      do i=1,min(cnt,cnt2)*MPI_STUBS_length(b)
         c(i+displ(1))=a(i+displ(1))
      enddo
      ierr=MPI_SUCCESS
      return
      end
!
      subroutine MPI_get_count
      write(6,*) 'MPI_get_count not authorized in serial-mpi mode, ABORT'
      call ABORT
      return
      end
!
      subroutine MPI_GROUP_incl(pe_gr_wcomm, pe_dommtot,proc_indomm,pe_gr_indomm,ierr)
      implicit none
      include 'mpif.h'
      integer pe_gr_wcomm, pe_dommtot,proc_indomm,pe_gr_indomm,ierr
      pe_gr_indomm = 0
      ierr = MPI_SUCCESS
      return
      end
!
      subroutine MPI_GROUP_RANK(comm,pe_me,ierr)
      implicit none
!
      include 'mpif.h'
!
      integer comm,pe_me,ierr
      ierr=MPI_SUCCESS
      pe_me=0
      return
      end
!
      subroutine MPI_GROUP_SIZE(comm,pe_tot,ierr)
      implicit none
!
      include 'mpif.h'
!
      integer comm,pe_tot,ierr
      ierr=MPI_SUCCESS
      pe_tot=1
      return
      end
!
      subroutine MPI_init(ierr)
      implicit none
      integer pe_tot
      save pe_tot
!
      include 'mpif.h'
!
      integer ierr
      data pe_tot / -1/
      ierr=0
!     if(pe_tot .ne. -1) call ABORT
      pe_tot=0
      return
      end
!
      subroutine MPI_INITIALIZED(mpi_started,ierr)
      implicit none
!
      include 'mpif.h'
!
      logical mpi_started
      integer ierr
!
      ierr=MPI_SUCCESS
      mpi_started=.false.
      return
      end
!
      subroutine MPI_irecv
      write(6,*) 'MPI_irecv not authorized in serial-mpi mode, ABORT'
      call ABORT
      return
      end
!
      subroutine MPI_isend
      write(6,*) 'MPI_isend not authorized in serial-mpi mode, ABORT'
      call ABORT
      return
      end
!
      subroutine MPI_recv
      write(6,*) 'MPI_recv not authorized in serial-mpi mode, ABORT'
      call ABORT
      return
      end
!
      subroutine MPI_reduce(send,recv,ni,l,m,n,o,ierr)
      implicit none
!
      include 'mpif.h'
!
      integer MPI_STUBS_length
      integer i,l,m,n,o,ierr
      integer send(*),recv(*),ni
      do i=1,ni*MPI_STUBS_length(l)
        recv(i)=send(i)
      enddo
      ierr=MPI_SUCCESS
      return
      end
!
      subroutine MPI_scatterv
      write(6,*) 'MPI_scatterv not authorized in serial-mpi mode, ABORT'
      call ABORT
      return
      end
!
      subroutine MPI_send
      write(6,*) 'MPI_send not authorized in serial-mpi mode, ABORT'
      call ABORT
      return
      end
!
      subroutine MPI_sendrecv
      write(6,*)  'MPI_sendrecv not authorized in serial-mpi mode, ABORT'
      call ABORT
      return
      end
!
      subroutine MPI_ssend
      write(6,*) 'MPI_ssend not authorized in serial-mpi mode, ABORT'
      call ABORT
      return
      end
!
      subroutine MPI_comm_free
      write(6,*) 'MPI_comm_free not authorized in serial-mpi mode, ABORT'
      call ABORT
      return
      end
!
      subroutine MPI_type_get_extent(dtyp_m,extent,ierr)  ! replaces MPI_type_extent
      integer,intent(IN) :: dtyp_m
      !integer, intent(out) :: ierr
      integer :: ierr
      integer :: extent
      write(6,*) 'MPI_type_get_extent not authorized in serial-mpi mode, ABORT'
      call ABORT
      ierr = -1
      return
      end
      subroutine MPI_type_extent(dtyp_m,extent,ierr)  ! old, deprecated entry
      integer,intent(IN) :: dtyp_m
      !integer, intent(out) :: ierr
      integer :: ierr
      integer :: extent
      write(6,*) 'MPI_type_extent not authorized in serial-mpi mode, ABORT'
      call ABORT
      ierr = -1
      return
      end
!
      integer function MPI_STUBS_length(itype)
      implicit none
      include 'mpif.h'
      integer itype

      MPI_STUBS_length=0

      if(itype .eq. MPI_DOUBLE_PRECISION) MPI_STUBS_length=2
      if(itype .eq. MPI_2DOUBLE_PRECISION) MPI_STUBS_length=2
      if(itype .eq. MPI_REAL8) MPI_STUBS_length=2
      if(itype .eq. MPI_INTEGER8) MPI_STUBS_length=2

      if(itype .eq. MPI_LOGICAL) MPI_STUBS_length=1

      if(itype .eq. MPI_INTEGER) MPI_STUBS_length=1
      if(itype .eq. MPI_INTEGER4) MPI_STUBS_length=1
      if(itype .eq. MPI_2INTEGER) MPI_STUBS_length=1
      if(itype .eq. MPI_REAL) MPI_STUBS_length=1
      if(itype .eq. MPI_REAL4) MPI_STUBS_length=1
      if(itype .eq. MPI_2REAL) MPI_STUBS_length=1

        if(MPI_STUBS_length.eq.0)then
        write(6,*)'MPI_STUBS_length ERROR: unrecognized type'
        call abort
      endif
      return
      end
!
      subroutine MPI_wait
      return
      end
!
      subroutine MPI_waitall
      return
      end
!
      REAL*8 function  MPI_wtick()
      MPI_wtick = 1.0E-9
      return
      end
      REAL*8 function  PMPI_wtick()
      PMPI_wtick = 1.0E-9
      return
      end
!
      REAL*8 function  MPI_wtime()
      real *8, save :: dummy_time=1.0E-9
      MPI_wtime=dummy_time
      dummy_time=dummy_time+1.0E-9
      return
      end
      REAL*8 function  PMPI_wtime()
      real *8, save :: dummy_time=1.0E-9
      PMPI_wtime=dummy_time
      dummy_time=dummy_time+1.0E-9
      return
      end

subroutine mpi_accumulate
  write(6,*) 'ERROR: mpi_accumulate not authorized in serial-mpi mode, ABORT'
  call ABORT
  return
end

subroutine mpi_alloc_mem
  write(6,*) 'ERROR: mpi_alloc_mem not authorized in serial-mpi mode, ABORT'
  call ABORT
  return
end

subroutine mpi_fetch_and_op
  write(6,*) 'ERROR: mpi_fetch_and_op not authorized in serial-mpi mode, ABORT'
  call ABORT
  return
end

subroutine mpi_free_mem
  write(6,*) 'ERROR: mpi_free_mem not authorized in serial-mpi mode, ABORT'
  call ABORT
  return
end

subroutine mpi_get
  write(6,*) 'ERROR: mpi_get not authorized in serial-mpi mode, ABORT'
  call ABORT
  return
end

subroutine mpi_group_free
  write(6,*) 'ERROR: mpi_group_free not authorized in serial-mpi mode, ABORT'
  call ABORT
  return
end

subroutine mpi_put
  write(6,*) 'ERROR: mpi_put not authorized in serial-mpi mode, ABORT'
  call ABORT
  return
end

subroutine mpi_scatter
  write(6,*) 'ERROR: mpi_scatter not authorized in serial-mpi mode, ABORT'
  call ABORT
  return
end

subroutine mpi_win_allocate
  write(6,*) 'ERROR: mpi_win_allocate not authorized in serial-mpi mode, ABORT'
  call ABORT
  return
end

subroutine mpi_win_complete
  write(6,*) 'ERROR: mpi_win_complete not authorized in serial-mpi mode, ABORT'
  call ABORT
  return
end

subroutine mpi_win_create
  write(6,*) 'ERROR: mpi_win_create not authorized in serial-mpi mode, ABORT'
  call ABORT
  return
end

subroutine mpi_win_fence
  write(6,*) 'ERROR: mpi_win_fence not authorized in serial-mpi mode, ABORT'
  call ABORT
  return
end

subroutine mpi_win_free
  write(6,*) 'ERROR: mpi_win_free not authorized in serial-mpi mode, ABORT'
  call ABORT
  return
end

subroutine mpi_win_lock
  write(6,*) 'ERROR: mpi_win_lock not authorized in serial-mpi mode, ABORT'
  call ABORT
  return
end

subroutine mpi_win_post
  write(6,*) 'ERROR: mpi_win_post not authorized in serial-mpi mode, ABORT'
  call ABORT
  return
end

subroutine mpi_win_start
  write(6,*) 'ERROR: mpi_win_start not authorized in serial-mpi mode, ABORT'
  call ABORT
  return
end

subroutine mpi_win_unlock
  write(6,*) 'ERROR: mpi_win_unlock not authorized in serial-mpi mode, ABORT'
  call ABORT
  return
end

subroutine mpi_win_wait
  write(6,*) 'ERROR: mpi_win_wait not authorized in serial-mpi mode, ABORT'
  call ABORT
  return
end

