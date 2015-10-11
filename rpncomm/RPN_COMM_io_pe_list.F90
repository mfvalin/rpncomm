!
! helper routine used by create_ioset
! for now only method 0 is supported for IO PE dispersion
!
  subroutine RPN_COMM_make_io_pe_list(x,y,npes,pe_nx,pe_ny,method)  !InTf!
    implicit none
    integer, dimension(npes), intent(OUT) :: x  !InTf!   ! x coordinates of PEs in set
    integer, dimension(npes), intent(OUT) :: y  !InTf!   ! y coordinates of PEs in set
    integer, intent(IN) :: npes                 !InTf!   ! number of PEs in set
    integer, intent(IN) :: pe_nx, pe_ny         !InTf!   ! number of PEs along x and y in PE grid
    integer, intent(IN) :: method               !InTf!   ! fill method
    integer :: i
    integer :: deltax, deltay, pe_nxy
    integer, save :: scramblex = 0
    integer, save :: scrambley = 0
    integer, dimension(16), save :: primes = [ &
         5,      7,      11,     13,   &
        17,     19,      23,     29,   &
        31,     37,      41,     43,   &
        47,     53,      59,     61  ]
!
    if(scramblex == 0) then    ! initialize
      scramblex = 1
      scrambley = 1
      if(pe_nx > pe_ny) then   ! scramblex = lowest number that is prime with respect to pe_nx
        do i = 1 , size(primes)
          scramblex = primes(i)
          if(mod(pe_nx,primes(i)) .ne. 0) exit
        enddo
      else                     ! scrambley = lowest number that is prime with respect to pe_ny
        do i = 1 , size(primes)
          scrambley = primes(i)
          if(mod(pe_ny,primes(i)) .ne. 0) exit
        enddo
      endif
!      print *,"DEBUG: scramblex, scrambley =",scramblex, scrambley
    endif
    x = -1
    y = -1
    if(method .ne. 0) return      ! method 0 is the only one supported for the time being
    deltax = 1
    deltay = 1
    pe_nxy = min(pe_nx,pe_ny)
    if(method == 0) then
      deltax = scramblex
      deltay = scrambley
    endif
    if( npes > pe_nxy * pe_nxy ) return
    do i = 0 , npes-1
      if(npes > pe_nxy) then
        if(pe_nx > pe_ny) then
          x(i+1) = mod( i * deltax , pe_nx )
          y(i+1) = mod( mod( i , pe_ny) + i / pe_nx , pe_ny)
        else
          x(i+1) = mod( mod( i , pe_nx ) + i / pe_ny , pe_nx)
          y(i+1) = mod( i * deltay , pe_ny)
        endif
      else
        x(i+1) = mod( i * deltax , pe_nx )
        y(i+1) = mod( i * deltay , pe_ny )
      endif
    enddo
  end subroutine RPN_COMM_make_io_pe_list  !InTf!
!
!
! check that this PE set is valid (no duplicates) and print IO PE map
! also check that no column has 2 IO PEs in a group and neither has any row
! used to validate what create_ioset has produced
!
  function RPN_COMM_check_ioset(newset, x ,y, npes, pe_nx, pe_ny, pe_me, diag) result(status)  !InTf!
    implicit none
    integer, intent(IN) :: newset, npes        !InTf!   ! set number, number of PEs in set
    integer, intent(IN) :: pe_nx, pe_ny        !InTf!   ! number of PEs along x and y
    integer, intent(IN) :: pe_me               !InTf!   ! ordinal in grid of this PE
    logical, intent(IN) :: diag                !InTf!   ! if .true., PE 0 prints the IO PE map for this set
    integer, intent(IN), dimension(npes) :: x  !InTf!   ! x coordinates of IO PEs
    integer, intent(IN), dimension(npes) :: y  !InTf!   ! y coordinates of IO PEs
    integer :: status                          !InTf!   ! 0 if set is OK, -1 otherwise
!
    integer*1, dimension(0:pe_nx-1,0:pe_ny-1) :: safe
    integer*1, dimension(0:pe_ny-1) :: row      ! there are pe_ny rows
    integer*1, dimension(0:pe_nx-1) :: col      ! there are pe_nx columns
    integer :: i, j, groupsize, low, high
!
    status = -1
    safe = 0
    groupsize = min(pe_nx, pe_ny)
    do i = 1, npes, groupsize  ! loop over groups
      row = 0
      col = 0
      do j = i , min(i+groupsize-1,npes)      ! for PE at coordinates ( x(j), y(j) )
        row(y(j)) = row(y(j)) + 1             ! add one to count for this row
        col(x(j)) = col(x(j)) + 1             ! add one to count for this column
      enddo
      if(any(row > 1) .or. any(col > 1) ) then
        print *,"ERROR: problem creating IO set, there are 2 or more PEs on row or column"
        print *,"ERROR: x = ",x
        print *,"ERROR: row = ",row
        print *,"ERROR: y = ",y
        print *,"ERROR: col = ",col
        print *,"ERROR: group, group limits = ",(i-1)/groupsize+1,i,min(i+groupsize-1,npes)
        status = -1
        return
      endif
    enddo
    do i = 1 , npes
      if(safe(x(i),y(i)) .ne. 0) then  ! OOPS, 2 PEs in this set are the same
        print *,"ERROR: problem creating IO set, there are duplicates"
        status = -1
        return
      else
        safe(x(i),y(i)) = 1 + mod(  ( (i-1) / min(pe_nx, pe_ny) ) , 9)  ! group number, 9 if group number > 9
      endif
    enddo
    if(pe_me == 0 .and. diag)then
      print 101,"===== IO PE map for set",newset," (",(npes+min(pe_nx, pe_ny)-1)/min(pe_nx, pe_ny)," groups) ====="
      print 102,"INFO: x =",x
      print 102,"INFO: y =",y
      do j = pe_ny-1 , 0 , -1
        print 100,j,safe(:,j)
      enddo
    endif
100   format(1X,I4,1x,128I1)
101   format(A,I3,A,I3,A)
102   format(A,15I4)
    status = 0
    return
  end function RPN_COMM_check_ioset  !InTf!
