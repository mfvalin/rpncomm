
function RPN_COMM_debug_str(name, string) result(status)   ! from :name=value extract value. return "" if :name= not found
  implicit none
  character (len=*), intent(IN) :: name
  character (len=*), intent(OUT) :: string
  integer :: status

  character(len=1024), save :: master=""
  integer, save :: lmaster=-1
  integer :: pos1, pos2, lstring

  if (lmaster == -1) then
    call get_environment_variable("RPN_COMM_DEBUG",master,lmaster,status)
    if(status .ne. 0) lmaster = 0
  endif

  string = ""
  status = -1
  if(lmaster == 0) return   ! environment variable was not defined

  pos1 = index(master,':'//trim(name)//'=')
  if(pos1 == 0) return   ! name not found

  do while(master(pos1:pos1) .ne. '=')  ! skip until after =
    pos1 = pos1 + 1
  enddo
  pos1 = pos1 + 1

  pos2 = 0
  do while(master(pos1:pos1) .ne. ':' .and. pos1 <= lmaster .and. pos2 < lstring) ! copy until end of master or next : or lstring characters
    pos2 = pos2 + 1
    string(pos2:pos2) = master(pos1:pos1)
    pos1 = pos1 + 1
  enddo
  if(pos2 <= lstring) status = 0

  return
end function RPN_COMM_debug_str

function RPN_COMM_debug_i(name,flag) result(status)
  implicit none
  character (len=*), intent(IN) :: name
  integer, intent(OUT) :: flag
  integer :: status

  integer, external :: RPN_COMM_debug_str
  character (len=32) :: string

  flag = 0
  status = RPN_COMM_debug_str(name,string)
  read(string,*)flag

  return
end function RPN_COMM_debug_i

function RPN_COMM_debug_r(name,flag) result(status)
  implicit none
  character (len=*), intent(IN) :: name
  real*4, intent(OUT) :: flag
  integer :: status

  integer, external :: RPN_COMM_debug_str
  character (len=32) :: string

  flag = 0.0
  status = RPN_COMM_debug_str(name,string)
  read(string,*)flag

  return
end function RPN_COMM_debug_r

function RPN_COMM_debug_l(name,flag) result(status)
  implicit none
  character (len=*), intent(IN) :: name
  logical, intent(OUT) :: flag
  integer :: status

  integer, external :: RPN_COMM_debug_str
  character (len=32) :: string

  flag = .false.
  status = RPN_COMM_debug_str(name,flag)
  read(string,*)flag

  return
end function RPN_COMM_debug_l

function RPN_COMM_debug_s(name,flag) result(status)
  implicit none
  character (len=*), intent(IN) :: name
  character (len=*), intent(OUT) :: flag
  integer :: status

  integer, external :: RPN_COMM_debug_str

  flag = ' '
  status = RPN_COMM_debug_str(name,flag)

  return
end function RPN_COMM_debug_s
  