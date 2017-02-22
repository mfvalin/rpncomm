module RPN_COMM_io  ! file replicator across nodes , async file copy, test and helper routines
  use iso_c_binding
  implicit none

#include "RPN_COMM_types.inc"
#define IN_RPN_COMM_io
#include "RPN_COMM_interfaces_int.inc"

  logical, save :: rpn_comm_io_debug=.false.

  contains

    integer function rpn_comm_unlink(name)
!
! delete a file
!
    implicit none
    character (len=*), intent(IN) ::name
    character (len=1), dimension(len_trim(name)+1), target :: temp
    
    temp = transfer( trim(name)//achar(0) , temp )
    rpn_comm_unlink=c_rpn_comm_unlink(c_loc(temp))
  end function rpn_comm_unlink

  integer function rpn_comm_open(name, mode)
!
! open a file for reading or writing
! mode = 0 : read
! mode = 1 : write
!
    implicit none
    character (len=*), intent(IN) ::name
    integer, intent(IN) :: mode
    character (len=1), dimension(len_trim(name)+1), target :: temp
    
    temp = transfer( trim(name)//achar(0) , temp )
    rpn_comm_open=c_rpn_comm_open(c_loc(temp),mode)
  end function rpn_comm_open

  integer function RPN_COMM_file_copy_test()  ! test asynchronous copying in background and file replication across nodes
    implicit none
    include 'mpif.h'
    integer :: fd1, fd2, status, id, ierr, rank
    integer, external :: RPN_COMM_file_copy_start, RPN_COMM_bcst_to_local_file

    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)                         ! my rank in world
    rpn_comm_io_debug=.true.

    RPN_COMM_file_copy_test = -1
    status = rpn_comm_copy(-1,-1,1)  ! set verbose debug mode

    if(rank /= 0) goto 222  ! first part of test for PE 0 only, asynchronous copying in background 

    print *,"START OF COPY TEST"
    fd1 = rpn_comm_open("existing_input_file",0) ! open for read
    if(fd1 < 0) then
      print *,"ERROR: cannot open existing_input_file for read or open error"
      return
    endif
    fd2 = rpn_comm_open("new_output_file",0) ! open for read
    if(fd2 >0 ) then
      print *,"WARNING: new_output_file exists"
    endif
    fd2 = rpn_comm_open("new_output_file",1) ! open for write
    if(fd2 <0 ) then
      print *,"ERROR: failed to open new_output_file for write"
      return
    endif
!
    status = rpn_comm_copy(fd1,fd2,0)  ! synchronous copy test
    if(status /= 0) then
      print *,"ERROR: synchronous copy failed"
      return
    endif
!
    id = rpn_comm_file_copy_start("existing_input_file","new_output_file_2")  ! asynchronous copy test
    print *,"INFO: asynchronous copy id =",id
    status = rpn_comm_wait(id)
    if(status /= 0) then
      print *,"ERROR: asynchronous copy failed"
      return
    endif
    print *,"INFO: cmp new_output_file new_output_file_2"
    call system("cmp new_output_file new_output_file_2 && echo ' INFO: files compare OK'")
!
    status = rpn_comm_unlink("new_output_file")  ! get rid of created files
    if(status /= 0) then
      print *,"WARNING: failed to delete new_output_file"
    endif
    status = rpn_comm_unlink("new_output_file_2")
    if(status /= 0) then
      print *,"WARNING: failed to delete new_output_file_2"
    endif
    print *,"END OF COPY TEST"

222 call mpi_barrier(MPI_COMM_WORLD,ierr)  ! all PEs wait for PE 0 to start replication test

    if(rank == 0) print *,"START OF REPLICATION TEST target_dir/existing_input_file"
    status = RPN_COMM_bcst_to_local_file("target_dir/existing_input_file","/dev/shm/LocalCopy",MPI_COMM_WORLD)
    if(rank == 0) print *,"INFO: replication status =",status
    if(rank == 0) print *,"END OF REPLICATION TEST"
    call mpi_finalize(ierr)
    return
  end function RPN_COMM_file_copy_test

end module RPN_COMM_io

integer function RPN_COMM_file_copy_start(name1,name2)  ! start asynchronous file copy  !InTf!
  use RPN_COMM_io
  implicit none
  character (len=*), intent(IN) :: name1, name2   ! copy file name1 into file name2     !InTf!
  integer :: fd1, fd2

  RPN_COMM_file_copy_start = 0  ! O.K. status
  fd1 = rpn_comm_open(name1,0)  ! open input file, get C file descriptor fd1
  if(fd1 < 0) goto 999
  fd2 = rpn_comm_open(name2,1)  ! open output file, get C file descriptor fd2
  if(fd2 < 0) goto 999

  RPN_COMM_file_copy_start = rpn_comm_copy(fd1,fd2,1)   ! copy fd1 onto fd2 asynchronously
  
  return

999 continue   ! error exception return
  RPN_COMM_file_copy_start = -1  ! error status
  return
end function RPN_COMM_file_copy_start                                                   !InTf!

integer function RPN_COMM_file_bcst(name,com)                                           !InTf!
!
! replicate a file across all the hosts in a communication domain
! PE 0 of communicator reads the file in chunks and broadcasts the chunks and their length
! one PE per host then writes the chunk into the local file ( except PE 0)
!
! this function is normally called after PE 0 has read the file from shared storage and written it to 
! fast local storage on its host if lname = name
!
! this is deemed more efficient than having all processes reading the same file on shared storage
! there will be a few processes on the same host reading the same local file subsequently
!
  implicit none
  character (len=*), intent(IN) :: name   ! name of the file to replicate               !InTf!
  integer, intent(IN) :: com              ! communicator                                !InTf!
  integer, external :: RPN_COMM_bcst_to_local_file

  RPN_COMM_file_bcst = RPN_COMM_bcst_to_local_file(name,name,com)

  return
end function RPN_COMM_file_bcst                                                         !InTf!

integer function RPN_COMM_bcst_to_local_file(name,lname,com)                                    !InTf!
!
! replicate a file across all the hosts in a communication domain
! PE 0 of communicator reads the file in chunks and broadcasts the chunks and their length
! one PE per host then writes the chunk into the local file ( except PE 0 if local name is the same as source file name)
!
! this is deemed more efficient than having all processes reading the same file on shared storage
! there will be a few processes on the same host reading the same local file subsequently
!
  use RPN_COMM_io
  implicit none
  character (len=*), intent(IN) :: name   ! name of the file to replicate ("slow" shared storage)         !InTf!
  character (len=*), intent(IN) :: lname  ! name of the replica (normally written to fast local storage)  !InTf!
  integer, intent(IN) :: com              ! communicator                                                  !InTf!

  interface
    integer(C_INT) function c_gethostid()BIND(C,name='gethostid')  ! interface to libc gethostid
      import :: C_INT
    end function c_gethostid
  end interface

  include 'mpif.h'

  integer(C_LONG_LONG), parameter :: NW8 = 1024*1024
  integer, dimension(NW8), target :: buf  ! 4 MegaByte buffer
  integer :: fdi, fdo, status, rank, rank_on_host, ierr, my_color, same_host, errors, sum_errors, com0, i
  integer(C_LONG_LONG) :: nwr, nww, nwf

  call mpi_comm_rank(com,rank,ierr)                         ! my rank in global communicator
  my_color = abs(c_gethostid())                             ! coloring by host identifier
  call MPI_COMM_SPLIT(com,my_color,rank,same_host,ierr)     ! same host communicator
  call mpi_comm_rank(same_host,rank_on_host,ierr)           ! my rank on host
  call MPI_COMM_SPLIT(com,rank_on_host,rank,com0,ierr) ! split global communicator by rank on host
  ! com0 will be used by rank_on_host=0 PEs to do the job

  fdi = 0   ! must initialize because some PEs may not call open
  fdo = 0   ! must initialize because some PEs may not call open
  if(rank == 0) then
    fdi = rpn_comm_open(name,0)   ! PE 0, open for read (file MUST exist)
    if(rpn_comm_io_debug) print 111,"rank=",rank," read  fdi=",fdi," file=",trim(name)," hostid=",my_color
    if(trim(lname) /= trim(name)) then
      fdo = rpn_comm_open(lname,1)   ! PE 0, open for write (if necessary)
      if(rpn_comm_io_debug) print 111,"rank=",rank," read  fdo=",fdo," file=",trim(lname)," hostid=",my_color
    endif
  else if(rank_on_host ==0) then ! open for write, one process per host, but not if same host as PE 0
    fdo = rpn_comm_open(name,1) 
    if(rpn_comm_io_debug) print 111,"rank=",rank," write fdo=",fdo," file=",trim(lname)," hostid=",my_color
  else                           ! nothing to do on this PE, neither global rank 0 nor local rank 0
    if(rpn_comm_io_debug) print 111,"rank=",rank," noop  fd=",-1," file=",trim(name)," hostid=",my_color
  endif

  errors = 0
  if(fdi < 0 .or. fdo < 0) errors = errors + 1
  call mpi_allreduce(errors,sum_errors,1,MPI_INTEGER,MPI_SUM,com,ierr)  !  get global eror count
  if(sum_errors > 0) goto 999  ! an open failed somewhere

  if(rank == 0) nwf = rpn_comm_file_size(fdi)             ! file size in Bytes

  errors = 0
  if(rank_on_host == 0) then                              ! reader/writer(s) association
    call mpi_bcast(nwf,1,MPI_INTEGER,0,same_host,ierr)    ! broadcast size to writers
    do i = 1,nwf,NW8*4                                    ! while something left to write
      nwr = min(NW8*4,nwf+1-i)                            ! size of block to read / broadcast
      if(rank == 0) nwr = rpn_comm_read(fdi,c_loc(buf(1)),nwr)   ! PE ranked at 0 reads 
      call mpi_bcast(buf,nwr,MPI_BYTE,0,com0,ierr)               ! and broadcasts
      if(rank == 0) then                ! global rank 0 only writes if lname different from name
        if(trim(lname) /= trim(name)) then
          nww = rpn_comm_write(fdo,c_loc(buf(1)),nwr)
          if(rpn_comm_io_debug) print 111,"DEBUG: rank no",rank," writing"
        else
          nww = nwr   ! not writing, let's not trigger errors
        endif
      else
        if(rpn_comm_io_debug) print 111,"DEBUG: rank no",rank," writing"
        nww = rpn_comm_write(fdo,c_loc(buf(1)),nwr)  ! one PE per node writes
      endif
      if(nww /= nwr) errors = errors + 1                             ! write error, OOPS
    enddo
  else
    if(rpn_comm_io_debug) print 111,"DEBUG: rank no",rank," doing nothing"
  endif
  call mpi_allreduce(errors,sum_errors,1,MPI_INTEGER,MPI_SUM,com,ierr)  ! I/O error anywhere ?
  if(sum_errors > 0) goto 999  ! yes, OUCH !!

  status = 0
  if(fdi > 0) status = status + rpn_comm_close(fdi)
  if(fdo > 0) status = status + rpn_comm_close(fdo)
  RPN_COMM_bcst_to_local_file = status   ! will be zero if no errors 
  return

111 format(A,I4,A,I4,A,A,A,Z9)
999 continue   ! error exception
  if(rank == 0) print 111,"ERROR: RPN_COMM_bcst_to_local_file errors =", sum_errors," rank =",rank
  RPN_COMM_bcst_to_local_file = sum_errors
  return

end function RPN_COMM_bcst_to_local_file                                                              !InTf!
