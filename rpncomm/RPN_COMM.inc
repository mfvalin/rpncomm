      interface
        function RPN_COMM_spread_context(context,com,rootpe,pe,npts)    &
     &                                 result(status)
        use ISO_C_BINDING
        implicit none
          type(c_ptr), intent(OUT) :: context
          character (len=*), intent(IN) :: com             ! RPN_COMM communicator
          integer, intent(IN) :: npts                      ! number of data points
          integer, intent(IN) :: rootpe                    ! root PE for the spread operation
          integer, dimension(npts), intent(IN) :: pe       ! destination table, data point i will be sent to PE pe(i)
          integer :: status                                ! 0 if successful, non zero otherwise
        end function RPN_COMM_spread_context

        function RPN_COMM_spread(context,source,npts,ndata,dest)        &
     &                           result(status)
        use ISO_C_BINDING
        implicit none
          type(c_ptr), intent(IN) :: context
          integer, intent(IN) :: npts, ndata
          real, dimension(npts,ndata), intent(IN) :: source  ! source array, used only on root PE
          real, dimension(:,:), pointer, intent(INOUT) :: dest
          integer :: status
        end function RPN_COMM_spread
      end interface
