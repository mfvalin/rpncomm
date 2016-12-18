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
!InTf!
        function rpn_comm_shmget(comm,size) result(where) bind(C,name='F_rpn_comm_shmget')    !InTf!
        import C_INT, C_PTR                                   !InTf!
        implicit none                                         !InTf!
          integer(C_INT), intent(IN), value :: comm           !InTf!   ! Fortran communicator (all members MUST be on same SMP node if size > 0)
          integer(C_INT), intent(IN), value :: size           !InTf!   ! size in bytes of shared memory area (if <0, split communicator to SMP node)
          type(C_PTR) :: where                                !InTf!   ! pointer to shared memory area
        end function rpn_comm_shmget                          !InTf!
!InTf!
        function rpn_comm_shmget_numa(comm,size) result(where) bind(C,name='F_rpn_comm_shmget_numa')    !InTf!
        import C_INT, C_PTR                                   !InTf!
        implicit none                                         !InTf!
          integer(C_INT), intent(IN), value :: comm           !InTf!   ! Fortran communicator (all members MUST be on same SMP node if size > 0)
          integer(C_INT), intent(IN), value :: size           !InTf!   ! size in bytes of shared memory area (if <0, split communicator to SMP node)
          type(C_PTR) :: where                                !InTf!   ! pointer to shared memory area
        end function rpn_comm_shmget_numa                     !InTf!
#endif

#define MAX_SHARED_SEGMENTS 128
static int shared_segments_max_index=-1;
static struct {
  struct {
    MPI_Comm base;     // base communicator
    MPI_Comm node;     // node communicator
    MPI_Comm socket;   // socket communicator
    MPI_Comm node0;    // node roots (rank 0) communicator
    MPI_Comm socket0;  // socket roots (rank 0) communicator
  } com;
  struct {
    int base;
    int node;
    int socket;
    int node0;
    int socket0;
  } ord;
//   struct {
//     int base;
//     int node;
//     int socket;
//     int node0;
//     int socket0;
//   } pop;
}table[MAX_SHARED_SEGMENTS];
static int myhostid = -1;

// get forced positive host id from gethostid()
static int my_hostid(){
  if(myhostid < 0) {
    myhostid = gethostid();
    myhostid &= 0x7FFFFFFF ;
  }
  return (myhostid);
}

unsigned long long rdtscp(int *socket, int *processor);

static int comm_lookup(MPI_Comm com){       // find communicator in table
  int i;

  for (i = 0 ; i <= shared_segments_max_index ; i++) {
    if(table[i].com.base == com) return(i);   // found in table, return index
  }
  if (shared_segments_max_index >= MAX_SHARED_SEGMENTS-1) return(-1);  // OOPS, table is full
  return(shared_segments_max_index+1);                  // not found in table, table not full, return index of free slot
}

#define BY_HOST    0
#define BY_SOCKET  2
#define BY_NUMA    1

// mode can be either BY_HOST or BY_NUMA or BY_SOCKET
static void *C_rpn_comm_shmget_(MPI_Comm c_comm_in, unsigned int shm_size, int mode)  /* allocate a shared memory segment */
{
 MPI_Comm c_comm, c_comm2;
  size_t size;                                   /* size of shared memory segment in bytes */
  int id;
  struct shmid_ds shm_buf;
  void *ptr;
  int ierr, myhost, myhost2, myrank, myrank2, myrank0;
  int current_slot, my_numa, my_core;
  unsigned long long dummy;

  current_slot = comm_lookup(c_comm_in);
  if(current_slot == -1){
    fprintf(stderr,"ERROR: (rpn_comm_shmget) shared segments table full \n");
    return NULL;    /* error, possibly fatal */
  }
  myhost = my_hostid();
  if(current_slot > shared_segments_max_index) {   // new slot, populate it

    table[current_slot].com.base = c_comm_in ;     // base communicator
    ierr = MPI_Comm_rank(c_comm_in,&myrank0) ;     // rank in base communicator
    table[current_slot].ord.base = myrank0 ;

    ierr = MPI_Comm_split(c_comm_in,myhost,myrank0,&c_comm) ;  // split base into node, color=hostid , weight=rank in base
    table[current_slot].com.node = c_comm ;        // intra node communicator
    ierr = MPI_Comm_rank(c_comm,&myrank) ;         // rank in intra node communicator
    table[current_slot].ord.node = myrank ;

    ierr = MPI_Comm_split(c_comm_in,myrank,myrank0,&c_comm2) ; // split base into noderoots, color=rank in node,  weight=rank in base
    table[current_slot].com.node0 = c_comm2 ;      // node roots (rank 0) communicator
    ierr = MPI_Comm_rank(c_comm2,&myrank2) ;       // rank in above group
    table[current_slot].ord.node0 = myrank2 ;

    if(mode != BY_HOST){    // split by NUMA space
#if defined(__x86_64__) &&  defined(__linux__)
      dummy = rdtscp(&my_numa, &my_core);  // on non X86 systems, per numaspace becomes per node
#else
      my_numa = 0;
      my_core = 0;
#endif
      ierr = MPI_Comm_split(c_comm_in,my_numa,myrank0,&c_comm) ;  // split node into socket, color=socket# , weight=rank in base
      table[current_slot].com.socket = c_comm ;    // intra socket communicator
      ierr = MPI_Comm_rank(c_comm,&myrank) ;       // rank in intra socket communicator
      table[current_slot].ord.socket = myrank ;

      ierr = MPI_Comm_split(c_comm_in,myrank,myrank0,&c_comm2) ; // split base into numaroots, color=rank in numa,  weight=rank in base
      table[current_slot].com.socket0 = c_comm2;   // socket roots (rank 0) communicator
      ierr = MPI_Comm_rank(c_comm2,&myrank2) ;     // rank in above group
      table[current_slot].ord.socket0 = myrank2;
    }else{
      table[current_slot].com.socket = MPI_COMM_NULL ;
      table[current_slot].ord.socket = -1 ;
      table[current_slot].com.socket0 = MPI_COMM_NULL ;
      table[current_slot].ord.socket0 = -1 ;
    }
  }else{                                           // existing slot
    if (mode == BY_HOST) {
      c_comm = table[current_slot].com.node ;
      myrank = table[current_slot].ord.node ;
    }else{
      c_comm = table[current_slot].com.socket ;
      myrank = table[current_slot].ord.socket ;
    }
  }
//
// at this point, 
// c_comm is the communicator for shared memory allocation (node or socket)
// myrank is the rank in the c_comm communicator
//
  size = shm_size;
  if(myrank == 0) {
    id=shmget(IPC_PRIVATE,size,IPC_CREAT|S_IRUSR|S_IWUSR);  /* rank 0 allocates shared memory segment */
    ptr=shmat(id,NULL,0);
//     shmctl(id,IPC_RMID,&shm_buf);      /* mark segment for deletion preventively (only works on linux) */
  }
  ierr=MPI_Bcast(&id,1,MPI_INTEGER,0,c_comm);                             /* all processes get segment id */
  if(id == -1) {
    if(myrank == 0) printf("ERROR: (rpn_comm_shmget) cannot create shared memory segment\n");
    return NULL;    /* error, possibly fatal */
  }else{
    if(myrank == 0) printf("INFO: (rpn_comm_shmget) created shared memory segment of size %d\n",shm_size);
  }

  if(myrank != 0) ptr=shmat(id,NULL,0);             /* all processes attach memory segment, rank 0 has already done it */

  if(ptr == NULL) printf("ERROR: (rpn_comm_shmget) got a null pointer from shmat, process = %d\n",myrank);
  myhost = (ptr == NULL);   /* this better be zero */
//  if(ptr != NULL) printf("DEBUG: (rpn_comm_shmget) process %d attached to segment\n",myrank);
//  printf("DEBUG: addr = %16p\n",ptr);
  ierr=MPI_Allreduce(&myhost,&myhost2,1,MPI_INTEGER,MPI_BOR,c_comm); /* boolean OR from all members of this comunicator */
//   ierr=MPI_Barrier(c_comm);                                          /* all processes have attached the segment */

  if(myrank == 0) shmctl(id,IPC_RMID,&shm_buf);      /* mark segment for deletion to make sure it is released when all processes terminate */
  if(0 != myhost2){                                                  /* not zero : OUCH */
    if(myrank == 0) fprintf(stderr,"ERROR: (rpn_comm_shmget) some processes were not able to attach to segment \n");
    if(ptr == NULL) ierr = shmdt(ptr);               /* detach from segment ia attached */
    return NULL;    /* error, possibly fatal */
    }
// bump shared_segments_max_index if new entry
  shared_segments_max_index = (shared_segments_max_index < current_slot) ?current_slot  : shared_segments_max_index;
  return ptr;                                        /* return pointer to shared memory area */
}

// all members of the c_comm_in communicator MUST be on the same node (same hostid)
void *C_rpn_comm_shmget(MPI_Comm c_comm_in, unsigned int shm_size)  /* allocate a shared memory segment by numa space */
{
  return C_rpn_comm_shmget_(c_comm_in, shm_size, BY_HOST);
}

// all members of the c_comm_in communicator MUST be on the same node (same hostid)
void *C_rpn_comm_shmget_numa(MPI_Comm c_comm_in, unsigned int shm_size)  /* allocate a shared memory segment by numa space */
{
  return C_rpn_comm_shmget_(c_comm_in, shm_size, BY_NUMA);
}

// the following Fortran functions expect MPI integer communicators, not RPN_COMM 'character string' communicators
void *F_rpn_comm_shmget(MPI_Fint f_comm, unsigned int shm_size)  /* allocate a shared memory segment on node (all numa spaces) */
{
//  MPI_Fint f_comm=comm;                                     /* all members of this communicator MUST be on same SMP node */
  return(C_rpn_comm_shmget(MPI_Comm_f2c(f_comm), shm_size));  /* translate Fortran communicator into C communicator before call */
}

void *F_rpn_comm_shmget_numa(MPI_Fint f_comm, unsigned int shm_size)  /* allocate a shared memory segment per  numa space*/
{
//  MPI_Fint f_comm=comm;                                          /* all members of this communicator MUST be on same SMP node */
  return(C_rpn_comm_shmget_numa(MPI_Comm_f2c(f_comm), shm_size));  /* translate Fortran communicator into C communicator before call */
}

#if defined(SELF_TEST)
main(int argc, char **argv){
//  int MPI_Init(int *argc, char ***argv)
  int ierr;
  int *sharedmem;
  int localrank;
  int hostid=gethostid();
  MPI_Comm comm = MPI_COMM_WORLD;
  int myrank;
  int my_numa, my_core;
  unsigned long long dummy;

  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&localrank);
//   hostid &= 0x7FFFFFFF;
//   ierr = MPI_Comm_split(MPI_COMM_WORLD,hostid,localrank,&comm);
  sharedmem = (int *) C_rpn_comm_shmget_numa(comm,1024*1024);
  myrank = table[0].ord.socket ;
  if(myrank==0) sharedmem[0] = localrank + 100000;
  ierr = MPI_Barrier(MPI_COMM_WORLD);
#if defined(__x86_64__) &&  defined(__linux__)
      dummy = rdtscp(&my_numa, &my_core);  // on non X86 systems, per numaspace becomes per node
#else
      my_numa = 0;
      my_core = 0;
#endif
  printf("Rank = %5.5d, Address = %p, contents = %d, core=%2.2d, numa=%d, hostid=%x\n",localrank,sharedmem,sharedmem[0],my_core,my_numa,hostid);
  ierr = MPI_Finalize();
}
#endif
