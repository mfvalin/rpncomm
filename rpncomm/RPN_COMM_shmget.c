#include <stdio.h>
#include <stdlib.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <mpi.h>
#include <unistd.h>


/* these Fortran functions need an explicit Fortran interface using ISO_C_BINDING */
/* because they return a C pointer and have value type arguments */

#ifdef MUST_NEVER_BE_TRUE
        integer, parameter :: BY_HOST = 0
        integer, parameter :: BY_NUMA = 1
        integer, parameter :: BY_SOCKET = 1
!InTf!
        function RPN_Comm_shm_bcast(tag,offset,count,root,comm) result(status) bind(C,name='F_RPN_Comm_shm_bcast')    !InTf!
        import C_INT                                          !InTf!
        implicit none                                         !InTf!
          integer(C_INT), intent(IN), value :: tag            !InTf!   ! shared memory area tag (from shmget)
          integer(C_INT), intent(IN), value :: offset         !InTf!   ! offset into shared memory area
          integer(C_INT), intent(IN), value :: count          !InTf!   ! number of elements to broadcast
          integer(C_INT), intent(IN), value :: root           !InTf!   ! rank of broadcast root PE in comm
          integer(C_INT), intent(IN), value :: comm           !InTf!   ! Fortran communicator (coherency with tag will be checked)
          integer(C_INT) :: status                            !InTf!   ! tag of shared memory area
        end function RPN_Comm_shm_bcast                       !InTf!
!InTf!
        function RPN_Comm_shmget(comm,size) result(tag) bind(C,name='F_RPN_Comm_shmget')    !InTf!
        import C_INT                                          !InTf!
        implicit none                                         !InTf!
          integer(C_INT), intent(IN), value :: comm           !InTf!   ! Fortran communicator
          integer(C_INT), intent(IN), value :: size           !InTf!   ! size in bytes of shared memory area
          integer(C_INT) :: tag                               !InTf!   ! tag of shared memory area
        end function RPN_Comm_shmget                          !InTf!
!InTf!
        function RPN_Comm_shmget_numa(comm,size) result(tag) bind(C,name='F_RPN_Comm_shmget_numa')    !InTf!
        import C_INT                                          !InTf!
        implicit none                                         !InTf!
          integer(C_INT), intent(IN), value :: comm           !InTf!   ! Fortran communicator
          integer(C_INT), intent(IN), value :: size           !InTf!   ! size in bytes of shared memory area
          integer(C_INT) :: tag                               !InTf!   ! tag of shared memory area
        end function RPN_Comm_shmget_numa                     !InTf!
!InTf!
        function RPN_Comm_shm_comm(tag) result(comm) bind(C,name='F_RPN_Comm_shm_comm')   !InTf!
        import C_INT                                          !InTf!
        implicit none                                         !InTf!
          integer(C_INT), intent(IN), value :: tag            !InTf!   ! shared memory area tag (from shmget)
          integer(C_INT) :: comm                              !InTf!   ! Fortran communicator associated with tag
        end function RPN_Comm_shm_comm                        !InTf!
!InTf!
        function RPN_Comm_shm_rank(tag) result(rank) bind(C,name='RPN_Comm_shm_rank')   !InTf!
        import C_INT                                          !InTf!
        implicit none                                         !InTf!
          integer(C_INT), intent(IN), value :: tag            !InTf!   ! shared memory area tag (from shmget)
          integer(C_INT) :: rank                              !InTf!   ! rank in communicator associated with tag
        end function RPN_Comm_shm_rank                        !InTf!
!InTf!
        function RPN_Comm_shm_mode(tag) result(mode) bind(C,name='RPN_Comm_shm_mode')   !InTf!
        import C_INT                                          !InTf!
        implicit none                                         !InTf!
          integer(C_INT), intent(IN), value :: tag            !InTf!   ! shared memory area tag (from shmget)
          integer(C_INT) :: mode                              !InTf!   ! mode (numa/node) associated with tag
        end function RPN_Comm_shm_mode                        !InTf!
!InTf!
        function RPN_Comm_shm_ptr(tag) result(mem) bind(C,name='RPN_Comm_shm_ptr')   !InTf!
        import C_INT, C_PTR                                   !InTf!
        implicit none                                         !InTf!
          integer(C_INT), intent(IN), value :: tag            !InTf!   ! shared memory area tag (from shmget)
          type(C_PTR) :: mem                                  !InTf!   ! pointer associated with tag
        end function RPN_Comm_shm_ptr                         !InTf!
!InTf!
        function RPN_Comm_shm_malloc(tag,size) result(mem) bind(C,name='RPN_Comm_shm_malloc')   !InTf!
        import C_INT, C_PTR                                   !InTf!
        implicit none                                         !InTf!
          integer(C_INT), intent(IN), value :: tag            !InTf!   ! shared memory area tag (from shmget)
          integer(C_INT), intent(IN), value :: size           !InTf!   ! size to allocate
          type(C_PTR) :: mem                                  !InTf!   ! pointer associated with tag
        end function RPN_Comm_shm_malloc                       !InTf!
#endif

#define MAX_SEGMENT_COMMS 16
static int shared_segments_max_index=-1;
static struct {
  struct {
    MPI_Comm base;     // base communicator
    MPI_Comm node;     // local node communicator
    MPI_Comm node0;    // rank 0 PEs in local node communicator (association of rank 0 PEs)
    MPI_Comm numa;     // local numa communicator
    MPI_Comm numa0;    // rank 0 PEs in local numa communicator (association of rank 0 PEs)
  } com;
  struct {             // rank of this PE in above communicators
    int base;
    int node;
    int node0;
    int numa;
    int numa0;
  } ord;
}com_table[MAX_SEGMENT_COMMS];      // table of communicators used for shared memory tag_table

#define MAX_SHARED_SEGMENTS 128
static int last_shared_segment = -1;
static struct {
  void *mem;       // address of shared memory segment
  int *top;        // address of last allocated area
  int tag;         // shared memory segment tag
  int mode;        // BY_HOST/BY_NUMA/BY_SOCKET
  int slot;        // slot in communicator table 
  int size;        // size of memory segment, in "int" units
  int used;        // allocated space in segment (same units as size)
} tag_table[MAX_SHARED_SEGMENTS];     // shared memory tag table

static long long now = 0;
#if defined(__x86_64__) &&  defined(__linux__)
static unsigned long long RDTSCP(int *socket, int *processor)  // get tsc/socket/processor
{
#if defined(__x86_64__) &&  defined(__linux__)
   unsigned int a, d, c;
   // rdtscp instruction
   // EDX:EAX contain TimeStampCounter
   // ECX contains IA32_TSC_AUX[31:0] (MSR_TSC_AUX value set by OS, lower 32 bits contain socket+processor)
   __asm__ volatile("rdtscp" : "=a" (a), "=d" (d), "=c" (c));
    *socket = (c & 0xFFF000)>>12;
    *processor = c & 0xFFF;

   return ((unsigned long long)a) | (((unsigned long long)d) << 32);;
#else
   *socket = 0;
   *processor = 0;
   return ++now;
#endif
}
#endif

static int ptr_lookup(void *ptr){       // find memory pointer in tag_table
  int i;

  for (i = 0 ; i <= last_shared_segment ; i++) {
    if(tag_table[i].mem == ptr) return(i);   // found in tag_table, return index
  }
  return(-1);                                 // not found in tag_table, return -1
}

static int tag_lookup(int tag){       // find tag in tag_table
  int i;

  for (i = 0 ; i <= last_shared_segment ; i++) {
    if(tag_table[i].tag == tag) return(i);   // found in tag_table, return index
  }
  return(-1);                                 // not found in tag_table, return -1
}

static int comm_lookup(MPI_Comm com){       // find communicator in com_table
  int i;

  for (i = 0 ; i <= shared_segments_max_index ; i++) {
    if(com_table[i].com.base == com) return(i);   // found in com_table, return index
  }
  if (shared_segments_max_index >= MAX_SEGMENT_COMMS-1) return(-1);  // OOPS, com_table is full
  return(shared_segments_max_index+1);                  // not found in com_table, com_table not full, return index of free slot
}

#define BY_HOST    0
#define BY_SOCKET  1
#define BY_NUMA    1

/*=================================== start of user callable functions =============================================*/
MPI_Comm RPN_Comm_shm_comm(int tag){   // get communicator associated with this shared memory segment
  int slot ;
  slot = tag_lookup(tag);
  if(slot <0 ) return (MPI_COMM_NULL);
  return tag_table[slot].mode == BY_NUMA ? 
         com_table[tag_table[slot].slot].com.numa : 
         com_table[tag_table[slot].slot].com.node ;
}
// Fortran interface, translates fortran communicator into c communicator
int F_RPN_Comm_shm_comm(int tag){
  return MPI_Comm_c2f(RPN_Comm_shm_comm(tag));
}

int RPN_Comm_shm_rank(int tag){   // get rank in communicator associated with this shared memory segment
  int slot ;
  slot = tag_lookup(tag);
  if(slot <0 ) return (-1);
  return tag_table[slot].mode == BY_NUMA ? 
         com_table[tag_table[slot].slot].ord.numa : 
         com_table[tag_table[slot].slot].ord.node ;
}

int RPN_Comm_shm_mode(int tag){   // get mode associated with this shared memory segment tag
  int slot ;
  slot = tag_lookup(tag);
  if(slot <0 ) return (-1);
  return tag_table[slot].mode ;
}

void *RPN_Comm_shm_ptr(int tag){   // get pointer associated with this shared memory segment tag
  int slot ;
  slot = tag_lookup(tag);
  if(slot <0 ) return (NULL);
  return tag_table[slot].mem ;
}

int *RPN_Comm_shm_malloc(int tag, int size){   // get pointer to new allocated space associated with this tag
  int slot ;
  int *current ;
  slot = tag_lookup(tag);
  if(slot <0 ) return (NULL);
  current = (int *)tag_table[slot].top ;
  if(tag_table[slot].used + size > tag_table[slot].size) return NULL;  // not enough space, error
  tag_table[slot].used += size ;    // update used
  tag_table[slot].top += size ;     // update top
  return current ;                  // return previous top
}

// broadcast a slice of a shared memory segment
int RPN_Comm_shm_bcast(int tag, int offset, int count, int root, MPI_Comm c_comm){
  int tag_slot = tag_lookup(tag) ;
  int comm_slot ;
  int *ptr;
  MPI_Comm comm ;
  MPI_Status status ;
  int ordinal, ierr ;

  if(tag_slot < 0) {
    fprintf(stderr,"ERROR: (RPN_Comm_shm_bcast) invalid shared segment tag \n");
    return -1 ;
  }
  comm_slot = tag_table[tag_slot].slot ;
  if(c_comm != com_table[comm_slot].com.base){
    fprintf(stderr,"ERROR: (RPN_Comm_shm_bcast) invalid base communicator for  shared segment\n");
    return -1 ;
  }
  ptr = (int *) tag_table[tag_slot].mem ;
  if(count+offset > tag_table[tag_slot].size || offset < 0 || count < 0){
    fprintf(stderr,"ERROR: (RPN_Comm_shm_bcast) invalid bounds for shared segment segment\n");
    return -1;
  }

  if(root != 0) {
    // tentative quick and dirty implementation
    // my pe is com_table[comm_slot].ord.base
    // if my pe is 0, recv from root
    if(com_table[comm_slot].ord.base == 0)    ierr = MPI_Recv(&ptr[offset], count, MPI_INTEGER, root, 0, c_comm, &status) ;
    // if my pe is root, send to 0
    if(com_table[comm_slot].ord.base == root) ierr = MPI_Send(&ptr[offset], count, MPI_INTEGER,    0, 0, c_comm) ;
    // otherwise, nothing to do but participate in later bcast
    fprintf(stderr,"ERROR: (RPN_Comm_shm_bcast) root != 0 not supported yet\n");
    return -1;
  }
  ordinal = tag_table[tag_slot].mode == BY_NUMA ?  // my ordinal in numa or node communicator
	    com_table[comm_slot].ord.numa :        // ordinal in local numa commnicator
	    com_table[comm_slot].ord.node ;        // ordinal in local node commnicator
  if(ordinal == 0) {                               // only for rank 0 PEs in numa/node communicator
    comm = tag_table[tag_slot].mode == BY_NUMA ?   // broadcast communicator
	    com_table[comm_slot].com.numa0 :       // by numa broadcast between rank 0 PEs of numa communicators
	    com_table[comm_slot].com.node0 ;       // by node broadcast between rank 0 PEs of node communicators
    ierr = MPI_Bcast(&ptr[offset], count, MPI_INTEGER, 0, comm) ;   // rank 0 PE of numa/node communicator to all members of rank 0 club
  }
  comm = tag_table[tag_slot].mode == BY_NUMA ?     // barrier communicator
	  com_table[comm_slot].com.numa :
	  com_table[comm_slot].com.node ;
  ierr = MPI_Barrier(comm) ;                       // intra numa/node barrier, wait for data delivery vi rank 0 of local numa/node
  return 0;
}

// Fortran interface, translates fortran communicator into c communicator
int F_RPN_Comm_shm_bcast(int tag, int offset, int count, int root, int f_comm){
  return RPN_Comm_shm_bcast(tag, offset, count, root, MPI_Comm_f2c(f_comm)) ;
}
/*==================================== end of user callable functions ==============================================*/

// mode can be either BY_HOST or BY_NUMA or BY_SOCKET
static int C_RPN_Comm_shmget(MPI_Comm c_comm_in, unsigned int shm_size, int mode)  /* allocate a shared memory segment */
{
 MPI_Comm c_comm, c_comm2, c_node ;
  size_t size;                                   /* size of shared memory segment in bytes */
  int id;
  struct shmid_ds shm_buf;
  void *ptr;
  int ierr, myhost, myhost0, myhost1, all_hosts, myrank, myrank2, myrank0;
  int current_slot, my_numa, my_core;
  unsigned long long dummy;

  if(last_shared_segment >= MAX_SHARED_SEGMENTS) {
    fprintf(stderr,"ERROR: (RPN_Comm_shmget) shared segment address table full \n");
    return -1;    /* error, possibly fatal */
  }
  current_slot = comm_lookup(c_comm_in);
  if(current_slot == -1){
    fprintf(stderr,"ERROR: (RPN_Comm_shmget) shared segment communicator table full \n");
    return -1;    /* error, possibly fatal */
  }
  myhost = gethostid();
  myhost0 = myhost & 0x7FFFFFFF  ; // lower 31 bits
  myhost1 = (myhost >> 31) & 0x1 ; // sign bit
  if(current_slot > shared_segments_max_index) {   // new slot, populate it

    com_table[current_slot].com.base = c_comm_in ;     // base communicator
    ierr = MPI_Comm_rank(c_comm_in,&myrank0) ;         // rank in base communicator
    com_table[current_slot].ord.base = myrank0 ;

    ierr = MPI_Comm_split(c_comm_in,myhost0,myrank0,&c_comm) ;  // split base using lower 31 bits of host id , weight=rank in base
    ierr = MPI_Comm_split(c_comm   ,myhost1,myrank0,&c_node) ;  // re split using upper bit of host id , weight=rank in base
    com_table[current_slot].com.node = c_node ;                 // intra node communicator
    ierr = MPI_Comm_rank(c_node,&myrank) ;                      // rank in intra node communicator
    com_table[current_slot].ord.node = myrank ;

    ierr = MPI_Comm_split(c_comm_in,myrank,myrank0,&c_comm2) ; // split base into noderoots, color=rank in node,  weight=rank in base
    com_table[current_slot].com.node0 = c_comm2 ;              // node roots (rank 0) communicator
    ierr = MPI_Comm_rank(c_comm2,&myrank2) ;                   // rank in above group
    com_table[current_slot].ord.node0 = myrank2 ;
    if(myrank != 0){                                           // not a member of rank 0 club
      com_table[current_slot].com.node0 = MPI_COMM_NULL ;
      com_table[current_slot].ord.node0 = -1 ;
    }

    dummy = RDTSCP(&my_numa, &my_core);  // on non X86 systems, per numaspace becomes per node

    ierr = MPI_Comm_split(c_node,my_numa,myrank0,&c_comm) ;     // re split node communicator into numa spaces, color=numa (my_numa), weight=rank in base
    com_table[current_slot].com.numa = c_comm ;                 // intra numa communicator
    ierr = MPI_Comm_rank(c_comm,&myrank) ;                      // rank in intra numa communicator
    com_table[current_slot].ord.numa = myrank ;

    ierr = MPI_Comm_split(c_comm_in,myrank,myrank0,&c_comm2) ; // split base into numaroots, color=rank in numa,  weight=rank in base
    com_table[current_slot].com.numa0 = c_comm2;               // numa roots (rank 0) communicator
    ierr = MPI_Comm_rank(c_comm2,&myrank2) ;                   // rank in above group
    com_table[current_slot].ord.numa0 = myrank2;
    if(myrank != 0){                                           // not a member of rank 0 club
      com_table[current_slot].com.numa0 = MPI_COMM_NULL ;
      com_table[current_slot].ord.numa0 = -1 ;
    }
  }

  if (mode == BY_HOST) {                      // get appropriate communicator for operatio
    c_comm = com_table[current_slot].com.node ;
    myrank = com_table[current_slot].ord.node ;
  }else{
    c_comm = com_table[current_slot].com.numa ;
    myrank = com_table[current_slot].ord.numa ;
  }
//
// at this point, 
// c_comm is the communicator for shared memory allocation (node or numa)
// myrank is the rank in the c_comm communicator
//
  size = shm_size * 1024;  // shm_size is in KBytes
  if(myrank == 0) {
    id=shmget(IPC_PRIVATE,size,IPC_CREAT|S_IRUSR|S_IWUSR);  /* rank 0 allocates shared memory segment */
    ptr=shmat(id,NULL,0);
#if defined(__linux__)
    shmctl(id,IPC_RMID,&shm_buf);      /* mark segment for deletion preventively (only works on linux) */
#endif
  }
  ierr=MPI_Bcast(&id,1,MPI_INTEGER,0,c_comm);                             /* all processes get segment id */
  if(id == -1) {
    if(myrank == 0) printf("ERROR: (RPN_Comm_shmget) cannot create shared memory segment\n");
    return -1;    /* error, possibly fatal */
  }else{
    if(myrank == 0) printf("INFO: (RPN_Comm_shmget) created shared memory segment of size %d\n",shm_size);
  }

  if(myrank != 0) ptr=shmat(id,NULL,0);             /* all processes attach memory segment, rank 0 has already done it */

  if(ptr == NULL) printf("ERROR: (RPN_Comm_shmget) got a null pointer from shmat, process = %d\n",myrank);
  myhost = (ptr == NULL);   /* this better be zero */
  ierr=MPI_Allreduce(&myhost,&all_hosts,1,MPI_INTEGER,MPI_BOR,c_comm);  /* boolean OR from all members of this comunicator */
//   ierr=MPI_Barrier(c_comm);                                          /* all processes should have attached the segment */

#if ! defined(__linux__)
  if(myrank == 0) shmctl(id,IPC_RMID,&shm_buf);      /* mark segment for deletion to make sure it is released when all processes terminate */
#endif
  if(0 != all_hosts){                                                  /* not zero : OUCH */
    if(myrank == 0) fprintf(stderr,"ERROR: (RPN_Comm_shmget) some processes were not able to attach to segment \n");
    if(ptr == NULL) ierr = shmdt(ptr);               /* detach from segment ia attached */
    return -1;                                     /* error, possibly fatal */
    }
// update tag table, bump last_shared_segment
  last_shared_segment++ ;
  tag_table[last_shared_segment].mem  = ptr;
  tag_table[last_shared_segment].top  = (int *) ptr;
  tag_table[last_shared_segment].tag  = id;
  tag_table[last_shared_segment].mode = mode;
  tag_table[last_shared_segment].slot = current_slot;
  tag_table[last_shared_segment].size = shm_size;
  tag_table[last_shared_segment].used = 0;    // nothing "allocated" in segment yet
// bump shared_segments_max_index if new entry
  shared_segments_max_index = (shared_segments_max_index < current_slot) ? current_slot  : shared_segments_max_index;
  return id;                                        /* return id of shared memory area */
}
/*=================================== start of user callable functions =============================================*/
int RPN_Comm_shmget(MPI_Comm c_comm_in, unsigned int shm_size)  /* allocate a shared memory segment by SMP node */
{
  return C_RPN_Comm_shmget(c_comm_in, shm_size, BY_HOST);
}

int RPN_Comm_shmget_numa(MPI_Comm c_comm_in, unsigned int shm_size)  /* allocate a shared memory segment by numa space */
{
  return C_RPN_Comm_shmget(c_comm_in, shm_size, BY_NUMA);
}

// the following Fortran callable functions expect MPI integer communicators, not RPN_COMM 'character string' communicators
int F_RPN_Comm_shmget(MPI_Fint f_comm, unsigned int shm_size)  /* allocate a shared memory segment on node (all numa spaces) */
{
  return(RPN_Comm_shmget(MPI_Comm_f2c(f_comm), shm_size));  /* translate Fortran communicator into C communicator before call */
}

int F_RPN_Comm_shmget_numa(MPI_Fint f_comm, unsigned int shm_size)  /* allocate a shared memory segment per  numa space*/
{
  return(RPN_Comm_shmget_numa(MPI_Comm_f2c(f_comm), shm_size));  /* translate Fortran communicator into C communicator before call */
}
/*===================================  end of user callable functions  =============================================*/

#if defined(SELF_TEST)
main(int argc, char **argv){
//  int MPI_Init(int *argc, char ***argv)
  int ierr;
  int *sharedmem;
  int localrank;
  int tag;
  int hostid=gethostid();
  MPI_Comm comm = MPI_COMM_WORLD;
  int myrank;
  int my_numa, my_core;
  unsigned long long dummy;

  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&localrank);
//   hostid &= 0x7FFFFFFF;
//   ierr = MPI_Comm_split(MPI_COMM_WORLD,hostid,localrank,&comm);
//   tag = RPN_Comm_shmget_numa(comm,1024*1024);
  tag = RPN_Comm_shmget_numa(comm,1024*1024);
  sharedmem = (int *)RPN_Comm_shm_ptr(tag) ;
  myrank = com_table[0].ord.numa ;
  if(myrank==0) sharedmem[0] = localrank + 110000;
  ierr = MPI_Barrier(MPI_COMM_WORLD);
  dummy = RDTSCP(&my_numa, &my_core);  // on non X86 systems, per numaspace becomes per node
  printf("Rank = %5.5d, Address = %p, contents = %d, core=%2.2d, numa=%d, hostid=%x, ord=%d, mode=%d\n",
	 localrank,sharedmem,sharedmem[0],my_core,my_numa,hostid,RPN_Comm_shm_rank(tag),RPN_Comm_shm_mode(tag));
  ierr = MPI_Finalize();
}
#endif
