      interface
                                                                        
      integer function RPN_COMM_chdir(string)                           
      implicit none                                                     
      character (len=*) :: string                                       
      end function RPN_COMM_chdir                                       
                                                                        
      subroutine RPN_COMM_env_var(varname,value)                        
      implicit none                                                     
      character (len=*), intent(IN) :: varname                          
      character (len=*), intent(OUT) :: value                           
      end subroutine RPN_COMM_env_var                                   
                                                                        
      integer function RPN_COMM_hostid()                                
      implicit none                                                     
      end function RPN_COMM_hostid                                      
                                                                        
        real*8 function RPN_COMM_wtime()                                
        implicit none                                                   
        end function RPN_COMM_wtime                                     
                                                                        
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
                                                                        
      SUBROUTINE RPN_COMM_bcastc(buffer,count,datatype,root,com,ierr)   
      implicit none                                                     
      integer, intent(IN) :: count,root                                 
      integer, intent(OUT) :: ierr                                      
      character(len=*), intent(INOUT) :: buffer                         
      character(len=*), intent(IN) :: datatype,com                      
      end SUBROUTINE RPN_COMM_bcastc                                    
                                                                        
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
                                                                        
      subroutine RPN_COMM_bloctopo(blocme,blocmex,blocmey,blocsizex,bloc&
     &sizey)                                                            
      implicit none                                                     
      integer, intent(out) :: blocme,blocmex,blocmey                    
      integer, intent(out) :: blocsizex, blocsizey                      
      end subroutine RPN_COMM_bloctopo                                  
                                                                        
      subroutine RPN_COMM_carac(npex,npey,me,medomm,mex,mey,sizex,sizey,&
     &ismaster, mymaster, mybloc, myblocx,myblocy,blocme,domname)       
      implicit none                                                     
      integer, intent(out) :: npex,npey,me,mex,mey,sizex,sizey,ismaster 
      integer, intent(out) :: mymaster, mybloc, myblocx,myblocy,blocme  
      integer, intent(out) :: medomm                                    
      character(len=*), intent(out) :: domname                          
      end subroutine RPN_COMM_carac                                     
                                                                        
      integer function RPN_COMM_colors(comm)                            
      implicit none                                                     
      character(len=*) :: comm                                          
      end function RPN_COMM_colors                                      
                                                                        
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
                                                                        
          subroutine RPN_COMM_const(IER)                                
          implicit none                                                 
          integer, intent(OUT) :: IER                                   
          end   subroutine RPN_COMM_const                               
                                                                        
        integer function RPN_COMM_datyp(data_int)                       
        implicit none                                                   
        character(len=*), intent(IN) :: data_int                        
        end function RPN_COMM_datyp                                     
                                                                        
      	SUBROUTINE RPN_COMM_defo(com)                                    
      	implicit none                                                    
      	character(len=*), intent(IN) ::  com                             
      	end SUBROUTINE RPN_COMM_defo                                     
                                                                        
      SUBROUTINE RPN_COMM_finalize(ierr)                                
      implicit none                                                     
      integer, intent(OUT) ::  ierr                                     
      end SUBROUTINE RPN_COMM_finalize                                  
                                                                        
        logical function RPN_COMM_grank(com)                            
        implicit none                                                   
        character(len=*), intent(IN) :: com                             
        end function RPN_COMM_grank                                     
                                                                        
        logical function RPN_COMM_ngrank(group)                         
        implicit none                                                   
        integer, intent(IN) :: group                                    
        end function RPN_COMM_ngrank                                    
                                                                        
      integer function RPN_COMM_group(com)                              
      implicit none                                                     
      character(len=*), intent(IN) :: com                               
      end function RPN_COMM_group                                       
                                                                        
      subroutine RPN_COMM_haloflip(g,minx,maxx,miny,maxy,ni,nj,nk,halox,&
     &haloy,gni)                                                        
        use rpn_comm                                                    
        implicit none                                                   
        integer, intent(IN) :: minx,maxx,miny,maxy,ni,nj,nk,halox,haloy,&
     &gni                                                               
        integer, intent(INOUT) :: g(minx:maxx,miny:maxy,nk)             
      end subroutine RPN_COMM_haloflip                                  
                                                                        
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
                                                                        
      integer function RPN_COMM_limit_2(my_id, npe, gmin, gmax,lmini,lma&
     &xi,count, offset,relax)                                           
      implicit none                                                     
      integer, intent(IN) ::  my_id, npe, gmin, gmax, relax             
      integer, intent(OUT) :: lmini,lmaxi                               
      integer, intent(OUT) :: count(npe),offset(npe)                    
      end function RPN_COMM_limit_2                                     
                                                                        
      integer function RPN_COMM_limit(my_id, npe, gmin, gmax, lmini,lmax&
     &i,count, offset)                                                  
      implicit none                                                     
      integer, intent(IN) ::  my_id, npe, gmin, gmax                    
      integer, intent(OUT) :: lmini,lmaxi                               
      integer, intent(OUT) :: count(npe),offset(npe)                    
      end function RPN_COMM_limit                                       
                                                                        
        SUBROUTINE RPN_COMM_UP2LOW(S1,S2)                               
        implicit none                                                   
        character (len=*), intent(IN) :: S1                             
        character (len=*), intent(OUT) :: S2                            
        END SUBROUTINE RPN_COMM_UP2LOW                                  
                                                                        
      	integer function RPN_COMM_mype(Me,Me_x,Me_y)                     
      	implicit none                                                    
      	integer, intent(OUT) ::  Me,Me_x,Me_y                            
      	end  function RPN_COMM_mype                                      
                                                                        
        integer function RPN_COMM_oper(op)                              
        implicit none                                                   
        character(len=*), intent(IN) :: op                              
        end function RPN_COMM_oper                                      
                                                                        
      	integer function RPN_COMM_option_L(name_S,value_L)               
      	implicit none                                                    
      	character(len=*), intent(IN) :: name_S                           
      	logical, intent(IN) :: value_L                                   
      	end function RPN_COMM_option_L                                   
                                                                        
      	integer function RPN_COMM_option(name_S,value)                   
      	implicit none                                                    
      	character(len=*), intent(IN) :: name_S                           
      	integer, intent(IN) :: value                                     
      	end function RPN_COMM_option                                     
                                                                        
      	SUBROUTINE RPN_COMM_optn(op_type,op_ival,op_rval,op_cval)        
      	implicit none                                                    
      	character(len=*), intent(IN) :: op_type, op_cval                 
      	integer, intent(IN) :: op_ival                                   
      	real *4, intent(IN) :: op_rval                                   
      	end SUBROUTINE RPN_COMM_optn                                     
                                                                        
      subroutine RPN_COMM_set_petopo(sizx,sizy)                         
      implicit none                                                     
      integer, intent(IN) :: sizx,sizy                                  
      end subroutine RPN_COMM_set_petopo                                
                                                                        
      integer function RPN_COMM_petopo(pex,pey)                         
      implicit none                                                     
      integer, intent(IN) :: pex,pey                                    
      end function RPN_COMM_petopo                                      
                                                                        
      integer function RPN_COMM_get_pe(x,y,grd,sgrd,communicator)       
      implicit none                                                     
      integer, intent(IN) :: x,y,grd,sgrd                               
      character *(*), intent(IN) :: communicator                        
      end function RPN_COMM_get_pe                                      
                                                                        
      	subroutine RPN_COMM_qadl(ipe,ni,nj,i0,in,j0,jn)                  
      	implicit none                                                    
          integer, intent(OUT) :: i0,in,j0,jn                           
          integer, intent(IN) :: ipe,ni,nj                              
      	end subroutine RPN_COMM_qadl                                     
                                                                        
      SUBROUTINE RPN_COMM_rank( com, rank ,ierr )                       
      implicit none                                                     
      character(len=*),intent(in) :: com                                
      integer, intent(out) :: rank                                      
      integer, intent(out) :: ierr                                      
      end SUBROUTINE RPN_COMM_rank                                      
                                                                        
        SUBROUTINE RPN_COMM_size( com, size ,ierr )                     
      implicit none                                                     
      integer, intent(OUT) :: ierr,size                                 
      character(len=*), intent(IN) ::  com                              
      end SUBROUTINE RPN_COMM_size                                      
                                                                        
      integer function RPN_COMM_split(mex,nx,nxg,minx,maxx,nxl,nxlmax,ha&
     &lox,nx0,fill)                                                     
      	implicit none                                                    
      	integer, intent(IN) :: nx, mex, nxg, halox                       
      	integer, intent(OUT) :: minx,maxx,nxl,nxlmax,nx0                 
      	logical, intent(IN) :: fill                                      
      	end function RPN_COMM_split                                      
                                                                        
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
                                                                        
        subroutine RPN_COMM_status_size(status_size)                    
        implicit none                                                   
        integer, intent(OUT) :: status_size                             
      end subroutine RPN_COMM_status_size                               
                                                                        
      SUBROUTINE RPN_COMM_tmg_wrt(Iun)                                  
      implicit none                                                     
      integer,intent(IN) :: Iun                                         
      end SUBROUTINE RPN_COMM_tmg_wrt                                   
                                                                        
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
                                                                        
      integer function RPN_COMM_topo_xy(nig,njg,mini,maxi,minj,maxj,nil,&
     &njl,haloi,haloj,peri,perj)                                        
      implicit none                                                     
      integer, intent(IN) :: nig,njg,haloi,haloj                        
      integer, intent(OUT) :: nil,njl,mini,maxi,minj,maxj               
      logical, intent(IN) ::  peri,perj                                 
      end function RPN_COMM_topo_xy                                     
                                                                        
        subroutine RPN_COMM_unit(unit)                                  
        implicit none                                                   
        integer, intent(IN) ::  unit                                    
        end subroutine RPN_COMM_unit                                    
                                                                        
        function rpn_comm_shmget(comm,size) result(where)               
        use ISO_C_BINDING                                               
        implicit none                                                   
          integer(C_INT), intent(IN) :: comm                            
          integer(C_INT), intent(IN) :: size                            
          type(C_PTR) :: where                                          
        end function rpn_comm_shmget                                    
      end interface
