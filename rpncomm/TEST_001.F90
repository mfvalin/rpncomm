subroutine rpn_comm_test_001
  use ISO_C_BINDING
  implicit none
  include 'RPN_COMM.inc'
  external :: RPN_COMM_grid_redist_test
  integer :: RPN_COMM_grid_redist_test
  external :: RPN_COMM_xch_halo_test
  integer :: RPN_COMM_xch_halo_test
  integer, external :: RPN_COMM_xch_halo_flip_test
  external TestUserInit, get_a_free_unit
  integer :: get_a_free_unit
  integer :: RPN_COMM_dist_test
  external RPN_COMM_dist_test
  integer :: Pelocal,Petotal,Pex,Pey,ierr,iun,test_to_perform
  integer :: nparams, i, ier, status
  integer, dimension(100) :: params
  character(len=256) :: RPN_COMM_TEST_CFG
  integer, external :: rpn_comm_2dgrid_test

  Pex = 0
  Pey = 0
!       UserInit supplied by TEST_helpers.f
  call RPN_COMM_init(TestUserInit,Pelocal,Petotal,Pex,Pey)
  print *,' Pelocal,Petotal,Pex,Pey =',Pelocal,Petotal,Pex,Pey

  call get_environment_variable("RPN_COMM_TEST_CFG",RPN_COMM_TEST_CFG,i,ier)
  if(ier==0) then
    read(RPN_COMM_TEST_CFG,FMT=*)test_to_perform,nparams,params(1:nparams)
  else
    iun=get_a_free_unit()
    open(UNIT=iun,FILE='TEST_001.cfg',STATUS='OLD')
    read(UNIT=iun,FMT=*)test_to_perform,nparams,params(1:nparams)
    close(UNIT=iun)
  endif
  if(IAND(test_to_perform,1)==1)then
    ierr=RPN_COMM_dist_test(Petotal)
  endif
  if(IAND(test_to_perform,2)==2)then
!          print *,'start grid_redist test'
    ierr=RPN_COMM_grid_redist_test(nparams,params)
  endif
  if(IAND(test_to_perform,4)==4)then
    print *,'start halo exchange test'
    ierr=RPN_COMM_xch_halo_test(nparams,params)
  endif
  if(IAND(test_to_perform,8)==8)then
    print *,'start haloflip exchange test'
    ierr=RPN_COMM_xch_halo_flip_test(nparams,params)
  endif
  if(IAND(test_to_perform,16)==16)then
    print *,'start fast distribution test'
    call RPN_COMM_fast_dist_test(nparams,params)
  endif
  if(IAND(test_to_perform,32)==32)then
    print *,'start shuffle distribution test'
    call RPN_COMM_io_dist_coll_test(nparams,params)
  endif
  if(IAND(test_to_perform,64)==64)then
    print *,'start 2D grid definition test'
    status = rpn_comm_2dgrid_test(nparams,params)
  endif
  if(IAND(test_to_perform,128)==128)then
    print *,'start one sided communications test'
    call RPN_COMM_i_win_test(nparams,params)
  endif
  if(IAND(test_to_perform,256)==256)then
    print *,'start shared memory segment test'
    call RPN_COMM_shm_test(nparams,params)
  endif
  if(IAND(test_to_perform,512)==512)then
    print *,'start one sided request/reply test'
    call RPN_COMM_i_win_req_test(nparams,params)
  endif
  call RPN_COMM_finalize(ierr)
  stop
end
subroutine TestUserInit(NX,NY) ! try to get NX,NY from file TEST.cfg if it exists
  external :: get_a_free_unit
  integer :: get_a_free_unit
  integer :: iun,ier,i
  character(len=128) :: RPN_COMM_TEST_SHAPE
  call get_environment_variable("RPN_COMM_TEST_SHAPE",RPN_COMM_TEST_SHAPE,i,ier)
  if(ier == 0) then
    read(RPN_COMM_TEST_SHAPE,*)NX,NY
    return
  endif
  iun=get_a_free_unit()
  if(iun<0)return
!        print *,'attempting to read TEST.cfg'
!        print *,'nx , ny =',nx,ny
  open(UNIT=iun,FILE='TEST.cfg',STATUS='OLD',ACTION='READ',IOSTAT=ier)
!        print *,'open iostat=',ier
  if(ier .ne. 0) then
!          print *,'attempting auto distribution, nx , ny =',nx,ny
    do i=7,2,-1
!            print *,'nx i mod(nx,i) =',nx, i, mod(nx,i)
      if(mod(NX,i) == 0 .and. nx .ne. i)then
        NY = NX / i
        NX = NX / NY
        return
      endif
    enddo
    return
  endif
  read(UNIT=iun,IOSTAT=ier,FMT=*)NX,NY
!        print *,'read iostat=',ier
  close(UNIT=iun)
  return
  end
  integer function get_a_free_unit()
  implicit none
  integer :: i
  character (len=16) :: access_mode
    get_a_free_unit=-1
    do i = 1 , 99  ! find an available unit number
      inquire(UNIT=i,ACCESS=access_mode)
      if(trim(access_mode) == 'UNDEFINED')then ! found
        get_a_free_unit = i
        exit
      endif
    enddo
  return
end
subroutine RPN_COMM_shm_test(nparams,params)
  use rpn_comm
  implicit none
  include 'RPN_COMM_interfaces.inc'
  integer, intent(IN) :: nparams
  integer, intent(IN), dimension(nparams) :: params

  type(C_PTR) :: my_ptr
  integer, parameter :: ni=80
  integer, parameter :: nj=40
  integer, parameter :: nk=62
  integer :: siz, i, j, k, ierr, ihost, errors
  integer, dimension(:,:,:,:), pointer :: array

  siz = ni*nj*nk*4*pe_tot_grid_host
  my_ptr = rpn_comm_shmget(pe_grid_host , siz)
  if( .not. c_associated(my_ptr) ) then
    print *,'ERROR: shared segment creation failed'
  else
    if(pe_me == 0) print *,'INFO: shared segment creation successful, all PEs attached'
    call c_f_pointer(my_ptr,array,[ni,nj,pe_tot_grid_host,nk])
  endif
  array(:,:,pe_me_grid_host+1,:) = -1
  do k=1,nk
  do j=1,nj
  do i=1,ni
    array(i,j,pe_me_grid_host+1,k) = pe_me_grid_host*1000000+i*10000+j*100+k
  enddo
  enddo
  enddo
  call mpi_barrier(pe_grid_host,ierr)
  if(pe_me_grid_host == 0) then
    errors = 0
    do ihost=0,pe_tot_grid_host-1
    do k=1,nk
    do j=1,nj
    do i=1,ni
      if( array(i,j,ihost+1,k) .ne. ihost*1000000+i*10000+j*100+k ) errors = errors + 1
    enddo
    enddo
    enddo
    enddo
    print *,'INFO: number of errors =',errors
  endif
end subroutine RPN_COMM_shm_test
