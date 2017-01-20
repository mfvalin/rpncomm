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
        function rpn_comm_shmget(comm,size) result(tag) bind(C,name='F_rpn_comm_shmget')    !InTf!
        import C_INT                                          !InTf!
        implicit none                                         !InTf!
          integer(C_INT), intent(IN), value :: comm           !InTf!   ! Fortran communicator
          integer(C_INT), intent(IN), value :: size           !InTf!   ! size in bytes of shared memory area
          integer(C_INT) :: tag                               !InTf!   ! tag of shared memory area
        end function rpn_comm_shmget                          !InTf!
!InTf!
        function rpn_comm_shmget_numa(comm,size) result(tag) bind(C,name='F_rpn_comm_shmget_numa')    !InTf!
        import C_INT                                          !InTf!
        implicit none                                         !InTf!
          integer(C_INT), intent(IN), value :: comm           !InTf!   ! Fortran communicator
          integer(C_INT), intent(IN), value :: size           !InTf!   ! size in bytes of shared memory area
          integer(C_INT) :: tag                               !InTf!   ! tag of shared memory area
        end function rpn_comm_shmget_numa                     !InTf!
!InTf!
        function rpn_comm_shm_rank(tag) result(rank) bind(C,name='rpn_comm_shm_rank')   !InTf!
        import C_INT                                          !InTf!
        implicit none                                         !InTf!
          integer(C_INT), intent(IN), value :: tag            !InTf!   ! shared memory area tag (from shmget)
          integer(C_INT) :: rank                              !InTf!   ! rank in communicator associated with tag
        end function rpn_comm_shm_rank                        !InTf!
!InTf!
        function rpn_comm_shm_mode(tag) result(mode) bind(C,name='rpn_comm_shm_mode')   !InTf!
        import C_INT                                          !InTf!
        implicit none                                         !InTf!
          integer(C_INT), intent(IN), value :: tag            !InTf!   ! shared memory area tag (from shmget)
          integer(C_INT) :: mode                              !InTf!   ! mode (numa/node) associated with tag
        end function rpn_comm_shm_mode                        !InTf!
!InTf!
        function rpn_comm_shm_ptr(tag) result(mem) bind(C,name='rpn_comm_shm_ptr')   !InTf!
        import C_INT, C_PTR                                   !InTf!
        implicit none                                         !InTf!
          integer(C_INT), intent(IN), value :: tag            !InTf!   ! shared memory area tag (from shmget)
          type(C_PTR) :: mem                                  !InTf!   ! pointer associated with tag
        end function rpn_comm_shm_ptr                         !InTf!
#endif

#define MAX_SHARED_SEGMENTS 128
static int shared_segments_max_index=-1;
static struct {
  struct {
    void *mem;
    MPI_Comm base;     // base communicator
    MPI_Comm node;     // node communicator
    MPI_Comm numa;     // numa communicator
    MPI_Comm node0;    // node roots (rank 0) communicator
    MPI_Comm numa0;    // numa roots (rank 0) communicator
  } com;
  struct {
    int base;
    int node;
    int numa;
    int node0;
    int numa0;
    int tag;
    int mode;
  } ord;
//   struct {
//     int base;
//     int node;
//     int numa;
//     int node0;
//     int numa0;
//   } pop;
}table[MAX_SHARED_SEGMENTS];
static int myhostid = -1;

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

static int slot_lookup(int tag){       // find tag in table
  int i;

  for (i = 0 ; i <= shared_segments_max_index ; i++) {
    if(table[i].ord.tag == tag) return(i);   // found in table, return index
  }
  return(-1);                                 // not found in table, return -1
}

static int comm_lookup(MPI_Comm com){       // find communicator in table
  int i;

  for (i = 0 ; i <= shared_segments_max_index ; i++) {
    if(table[i].com.base == com) return(i);   // found in table, return index
  }
  if (shared_segments_max_index >= MAX_SHARED_SEGMENTS-1) return(-1);  // OOPS, table is full
  return(shared_segments_max_index+1);                  // not found in table, table not full, return index of free slot
}

#define BY_HOST    0
#define BY_SOCKET  1
#define BY_NUMA    1

int rpn_comm_shm_rank(int tag){   // get rank in communicator associated with this shared memory segment
  int slot ;
  slot = slot_lookup(tag);
  if(slot <0 ) return (-1);
  return table[slot].ord.mode == BY_NUMA ? table[slot].ord.numa : table[slot].ord.node ;
}

int rpn_comm_shm_mode(int tag){   // get mode associated with this shared memory segment tag
  int slot ;
  slot = slot_lookup(tag);
  if(slot <0 ) return (-1);
  return table[slot].ord.mode ;
}

void *rpn_comm_shm_ptr(int tag){   // get pointer associated with this shared memory segment tag
  int slot ;
  slot = slot_lookup(tag);
  if(slot <0 ) return (NULL);
  return table[slot].com.mem ;
}

// mode can be either BY_HOST or BY_NUMA or BY_SOCKET
static int C_rpn_comm_shmget(MPI_Comm c_comm_in, unsigned int shm_size, int mode)  /* allocate a shared memory segment */
{
 MPI_Comm c_comm, t_comm, c_comm2;
  size_t size;                                   /* size of shared memory segment in bytes */
  int id;
  struct shmid_ds shm_buf;
  void *ptr;
  int ierr, myhost, myhost0, myhost1, all_hosts, myrank, myrank2, myrank0;
  int current_slot, my_numa, my_core;
  unsigned long long dummy;

  current_slot = comm_lookup(c_comm_in);
  if(current_slot == -1){
    fprintf(stderr,"ERROR: (rpn_comm_shmget) shared segments table full \n");
    return -1;    /* error, possibly fatal */
  }
  myhost = gethostid();
  myhost0 = myhost & 0x7FFFFFFF  ; // lower 31 bits
  myhost1 = (myhost >> 31) & 0x1 ; // sign bit
  if(current_slot > shared_segments_max_index) {   // new slot, populate it

    table[current_slot].com.base = c_comm_in ;     // base communicator
    ierr = MPI_Comm_rank(c_comm_in,&myrank0) ;     // rank in base communicator
    table[current_slot].ord.base = myrank0 ;

    ierr = MPI_Comm_split(c_comm_in,myhost0,myrank0,&t_comm) ;  // split base using lower 31 bits of host id , weight=rank in base
    ierr = MPI_Comm_split(t_comm   ,myhost1,myrank0,&c_comm) ;  // re split using upper bit of host id , weight=rank in base
    table[current_slot].com.node = c_comm ;        // intra node communicator
    ierr = MPI_Comm_rank(c_comm,&myrank) ;         // rank in intra node communicator
    table[current_slot].ord.node = myrank ;

    ierr = MPI_Comm_split(c_comm_in,myrank,myrank0,&c_comm2) ; // split base into noderoots, color=rank in node,  weight=rank in base
    table[current_slot].com.node0 = c_comm2 ;      // node roots (rank 0) communicator
    ierr = MPI_Comm_rank(c_comm2,&myrank2) ;       // rank in above group
    table[current_slot].ord.node0 = myrank2 ;

    dummy = RDTSCP(&my_numa, &my_core);  // on non X86 systems, per numaspace becomes per node

    ierr = MPI_Comm_split(c_comm_in,my_numa,myrank0,&c_comm) ;  // split node into numa spaces, color=numa (my_numa), weight=rank in base
    table[current_slot].com.numa = c_comm ;    // intra numa communicator
    ierr = MPI_Comm_rank(c_comm,&myrank) ;       // rank in intra numa communicator
    table[current_slot].ord.numa = myrank ;

    ierr = MPI_Comm_split(c_comm_in,myrank,myrank0,&c_comm2) ; // split base into numaroots, color=rank in numa,  weight=rank in base
    table[current_slot].com.numa0 = c_comm2;   // numa roots (rank 0) communicator
    ierr = MPI_Comm_rank(c_comm2,&myrank2) ;     // rank in above group
    table[current_slot].ord.numa0 = myrank2;
  }

  if (mode == BY_HOST) {                      // get appropriate communicator for operatio
    c_comm = table[current_slot].com.node ;
    myrank = table[current_slot].ord.node ;
  }else{
    c_comm = table[current_slot].com.numa ;
    myrank = table[current_slot].ord.numa ;
  }
//
// at this point, 
// c_comm is the communicator for shared memory allocation (node or numa)
// myrank is the rank in the c_comm communicator
//
  size = shm_size;
  if(myrank == 0) {
    id=shmget(IPC_PRIVATE,size,IPC_CREAT|S_IRUSR|S_IWUSR);  /* rank 0 allocates shared memory segment */
    ptr=shmat(id,NULL,0);
#if defined(__linux__)
    shmctl(id,IPC_RMID,&shm_buf);      /* mark segment for deletion preventively (only works on linux) */
#endif
  }
  ierr=MPI_Bcast(&id,1,MPI_INTEGER,0,c_comm);                             /* all processes get segment id */
  if(id == -1) {
    if(myrank == 0) printf("ERROR: (rpn_comm_shmget) cannot create shared memory segment\n");
    return -1;    /* error, possibly fatal */
  }else{
    if(myrank == 0) printf("INFO: (rpn_comm_shmget) created shared memory segment of size %d\n",shm_size);
  }

  if(myrank != 0) ptr=shmat(id,NULL,0);             /* all processes attach memory segment, rank 0 has already done it */

  if(ptr == NULL) printf("ERROR: (rpn_comm_shmget) got a null pointer from shmat, process = %d\n",myrank);
  myhost = (ptr == NULL);   /* this better be zero */
  ierr=MPI_Allreduce(&myhost,&all_hosts,1,MPI_INTEGER,MPI_BOR,c_comm);  /* boolean OR from all members of this comunicator */
//   ierr=MPI_Barrier(c_comm);                                          /* all processes should have attached the segment */

#if ! defined(__linux__)
  if(myrank == 0) shmctl(id,IPC_RMID,&shm_buf);      /* mark segment for deletion to make sure it is released when all processes terminate */
#endif
  if(0 != all_hosts){                                                  /* not zero : OUCH */
    if(myrank == 0) fprintf(stderr,"ERROR: (rpn_comm_shmget) some processes were not able to attach to segment \n");
    if(ptr == NULL) ierr = shmdt(ptr);               /* detach from segment ia attached */
    return -1;                                     /* error, possibly fatal */
    }
// bump shared_segments_max_index if new entry
  table[current_slot].com.mem = ptr ;
  table[current_slot].ord.tag = id ;
  table[current_slot].ord.mode = mode;
  shared_segments_max_index = (shared_segments_max_index < current_slot) ?current_slot  : shared_segments_max_index;
  return id;                                        /* return pointer to shared memory area */
}
/*=================================== start of user callable functions =============================================*/
int rpn_comm_shmget(MPI_Comm c_comm_in, unsigned int shm_size)  /* allocate a shared memory segment by SMP node */
{
  return C_rpn_comm_shmget(c_comm_in, shm_size, BY_HOST);
}

int rpn_comm_shmget_numa(MPI_Comm c_comm_in, unsigned int shm_size)  /* allocate a shared memory segment by numa space */
{
  return C_rpn_comm_shmget(c_comm_in, shm_size, BY_NUMA);
}

// the following Fortran callable functions expect MPI integer communicators, not RPN_COMM 'character string' communicators
int F_rpn_comm_shmget(MPI_Fint f_comm, unsigned int shm_size)  /* allocate a shared memory segment on node (all numa spaces) */
{
  return(rpn_comm_shmget(MPI_Comm_f2c(f_comm), shm_size));  /* translate Fortran communicator into C communicator before call */
}

int F_rpn_comm_shmget_numa(MPI_Fint f_comm, unsigned int shm_size)  /* allocate a shared memory segment per  numa space*/
{
  return(rpn_comm_shmget_numa(MPI_Comm_f2c(f_comm), shm_size));  /* translate Fortran communicator into C communicator before call */
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
//   tag = rpn_comm_shmget_numa(comm,1024*1024);
  tag = rpn_comm_shmget(comm,1024*1024);
  sharedmem = rpn_comm_shm_ptr(tag) ;
  myrank = table[0].ord.numa ;
  if(myrank==0) sharedmem[0] = localrank + 110000;
  ierr = MPI_Barrier(MPI_COMM_WORLD);
  dummy = RDTSCP(&my_numa, &my_core);  // on non X86 systems, per numaspace becomes per node
  printf("Rank = %5.5d, Address = %p, contents = %d, core=%2.2d, numa=%d, hostid=%x, ord=%d, mode=%d\n",
	 localrank,sharedmem,sharedmem[0],my_core,my_numa,hostid,rpn_comm_shm_rank(tag),rpn_comm_shm_mode(tag));
  ierr = MPI_Finalize();
}
#endif
