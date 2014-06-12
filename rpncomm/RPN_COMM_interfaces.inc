      interface
                                                                        
      subroutine RPN_COMM_carac(npex,npey,me,medomm,mex,mey,sizex,sizey,&
     &ismaster, mymaster, mybloc, myblocx,myblocy,blocme,domname)       
      implicit none                                                     
      integer, intent(out) :: npex,npey,me,mex,mey,sizex,sizey,ismaster 
      integer, intent(out) :: mymaster, mybloc, myblocx,myblocy,blocme  
      integer, intent(out) :: medomm                                    
      character(len=*), intent(out) :: domname                          
      end subroutine RPN_COMM_carac                                     
                                                                        
      integer function RPN_COMM_comm(com)                               
      implicit none                                                     
      character(len=*), intent(IN) :: com                               
      end function RPN_COMM_comm                                        
                                                                        
      integer function RPN_COMM_custom_comm(com,name,mode)              
      implicit none                                                     
      character(len=*), intent(IN) :: name                              
      integer, intent(IN) :: com                                        
      integer, intent(IN) :: mode                                       
      end function RPN_COMM_custom_comm                                 
                                                                        
        integer function RPN_COMM_datyp(data_int)                       
        implicit none                                                   
        character(len=*), intent(IN) :: data_int                        
        end function RPN_COMM_datyp                                     
                                                                        
      SUBROUTINE RPN_COMM_finalize(ierr)                                
      implicit none                                                     
      integer, intent(OUT) ::  ierr                                     
      end SUBROUTINE RPN_COMM_finalize                                  
                                                                        
      subroutine RPN_COMM_mydomain (call_back, mydomain)                
      implicit none                                                     
      external :: call_back                                             
      integer, intent(OUT) :: mydomain                                  
      end subroutine RPN_COMM_mydomain                                  
                                                                        
      subroutine RPN_COMM_world_set(world_comm)                         
      implicit none                                                     
      integer, intent(IN) ::  world_comm                                
      end subroutine RPN_COMM_world_set                                 
                                                                        
      SUBROUTINE RPN_COMM_init(Userinit,Pelocal,Petotal,Pex,Pey)        
      implicit none                                                     
      integer, intent(out)   :: Pelocal,Petotal                         
      integer, intent(inout) :: Pex,Pey                                 
      external Userinit                                                 
      end SUBROUTINE RPN_COMM_init                                      
                                                                        
      INTEGER FUNCTION RPN_COMM_init_multigrid(Userinit,Pelocal,Petotal,&
     &Pex,Pey,MultiGrids)                                               
      implicit none                                                     
      external :: Userinit                                              
      integer, intent(out)   :: Pelocal,Petotal                         
      integer, intent(inout) :: Pex,Pey                                 
      integer, intent(in)    :: MultiGrids                              
      end FUNCTION RPN_COMM_init_multigrid                              
                                                                        
      INTEGER FUNCTION RPN_COMM_init_multi_level(Userinit,Pelocal,Petota&
     &l,Pex,Pey,MultiGrids,Grids)                                       
      implicit none                                                     
      external :: Userinit                                              
      integer, intent(out)   :: Pelocal,Petotal                         
      integer, intent(inout) :: Pex,Pey                                 
      integer, intent(in)    :: MultiGrids                              
      integer, intent(in)    :: Grids                                   
      end FUNCTION RPN_COMM_init_multi_level                            
                                                                        
      integer function RPN_COMM_mype(Me,Me_x,Me_y)                      
      implicit none                                                     
      integer, intent(OUT) ::  Me,Me_x,Me_y                             
      end  function RPN_COMM_mype                                       
                                                                        
      SUBROUTINE RPN_COMM_rank( com, rank ,ierr )                       
      implicit none                                                     
      character(len=*),intent(in) :: com                                
      integer, intent(out) :: rank                                      
      integer, intent(out) :: ierr                                      
      end SUBROUTINE RPN_COMM_rank                                      
                                                                        
      integer function RPN_COMM_topo_2(nxg,minx,maxx,nxl,nxlmax,halo,nx0&
     &,alongx,fill,relax,abort)                                         
      implicit none                                                     
      integer, intent(in) :: nxg,halo,relax                             
      logical, intent(in) :: alongx,fill,abort                          
      integer, intent(out):: minx,maxx,nxl,nxlmax,nx0                   
        end  function RPN_COMM_topo_2                                   
                                                                        
      integer function RPN_COMM_topo(nxg,minx,maxx,nxl,nxlmax,halo,nx0,a&
     &longx,fill)                                                       
      implicit none                                                     
      integer, intent(in) :: nxg,halo                                   
      logical, intent(in) :: alongx,fill                                
      integer, intent(out):: minx,maxx,nxl,nxlmax,nx0                   
        end function RPN_COMM_topo                                      
                                                                        
        SUBROUTINE RPN_COMM_barrier(com,ierr)                           
        implicit none                                                   
        integer, intent(OUT) :: ierr                                    
        character(len=*), intent(IN) ::  com                            
      end SUBROUTINE RPN_COMM_barrier                                   
                                                                        
        integer function RPN_COMM_barrier_data(level,values,nvalues)    
        integer , intent(IN) :: level                                   
        integer , intent(IN) :: nvalues                                 
        real *8, dimension(nvalues), intent(OUT) :: values              
        end function RPN_COMM_barrier_data                              
                                                                        
      integer function RPN_COMM_bloc(nblocx,nblocy)                     
      implicit none                                                     
      integer, intent(IN) :: nblocx, nblocy                             
      end function RPN_COMM_bloc                                        
                                                                        
      integer function RPN_COMM_bloc_find(nblocx,nblocy,set)            
      implicit none                                                     
      integer, intent(IN) :: nblocx, nblocy                             
      logical, intent(IN) :: set                                        
      end function RPN_COMM_bloc_find                                   
                                                                        
      integer function RPN_COMM_bloc_create(nblocx,nblocy)              
      implicit none                                                     
      integer, intent(IN) :: nblocx, nblocy                             
      end function RPN_COMM_bloc_create                                 
                                                                        
      function RPN_COMM_spread(context,source,npts,ndata,dest) result(st&
     &atus)                                                             
        use ISO_C_BINDING                                               
        implicit none                                                   
        type(c_ptr), intent(IN) :: context                              
        integer, intent(IN) :: npts, ndata                              
        real, dimension(npts,ndata), intent(IN) :: source               
        real, dimension(:,:), pointer, intent(INOUT) :: dest            
        integer :: status                                               
      end function RPN_COMM_spread                                      
                                                                        
      function RPN_COMM_spread_context(context,com,rootpe,pe,npts) resul&
     &t(status)                                                         
        use ISO_C_BINDING                                               
        implicit none                                                   
        type(c_ptr), intent(OUT) :: context                             
        character (len=*), intent(IN) :: com                            
        integer, intent(IN) :: npts                                     
        integer, intent(IN) :: rootpe                                   
        integer, dimension(npts), intent(IN) :: pe                      
        integer :: status                                               
      end function RPN_COMM_spread_context                              
                                                                        
        function rpn_comm_shmget(comm,size) result(where)               
        use ISO_C_BINDING                                               
        implicit none                                                   
          integer(C_INT), intent(IN) :: comm                            
          integer(C_INT), intent(IN) :: size                            
          type(C_PTR) :: where                                          
        end function rpn_comm_shmget                                    
      end interface
