#include <stdio.h>
#include <stdlib.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <mpi.h>
#include <unistd.h>


/* this function needs an explicit fortran interface using ISO_C_BINDING */
/* because it returns a C pointer and has value type arguments */

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
        end function rpn_comm_shmget                          !InTf!
!InTf!
        function rpn_comm_shmget_all(comm,size) result(where) bind(C,name='F_rpn_comm_shmget_all')    !InTf!
        import C_INT, C_PTR                                   !InTf!
        implicit none                                         !InTf!
          integer(C_INT), intent(IN), value :: comm           !InTf!   ! Fortran communicator (all members MUST be on same SMP node if size > 0)
          integer(C_INT), intent(IN), value :: size           !InTf!   ! size in bytes of shared memory area (if <0, split communicator to SMP node)
          type(C_PTR) :: where                                !InTf!   ! pointer to shared memory area
        end function rpn_comm_shmget                          !InTf!
!InTf!
        function rpn_comm_shmget_all_numa(comm,size) result(where) bind(C,name='F_rpn_comm_shmget_all_numa')    !InTf!
        import C_INT, C_PTR                                   !InTf!
        implicit none                                         !InTf!
          integer(C_INT), intent(IN), value :: comm           !InTf!   ! Fortran communicator (all members MUST be on same SMP node if size > 0)
          integer(C_INT), intent(IN), value :: size           !InTf!   ! size in bytes of shared memory area (if <0, split communicator to SMP node)
          type(C_PTR) :: where                                !InTf!   ! pointer to shared memory area
        end function rpn_comm_shmget                          !InTf!
#endif

static int myrank;

void *C_rpn_comm_shmget(MPI_Comm c_comm_in, unsigned int shm_size)  /* allocate a shared memory segment */
{
 MPI_Comm c_comm;
  size_t size;                                   /* size of shared memory segment in bytes */
  int id;
  struct shmid_ds shm_buf;
  void *ptr;
  int ierr, myhost, myhost2;
  
  myhost=gethostid();
  if(shm_size > 0){
    size = shm_size;
    c_comm = c_comm_in;                   /* all members of this communicator MUST be on same SMP node */
    ierr=MPI_Allreduce(&myhost,&myhost2,1,MPI_INTEGER,MPI_BOR,c_comm); /* boolean OR from all members of this comunicator */
    myhost ^= myhost2;                                                 /* should be zero everywhere */
    ierr=MPI_Allreduce(&myhost,&myhost2,1,MPI_INTEGER,MPI_BOR,c_comm); /* boolean OR from all members of this comunicator */
    if(0 != myhost2){                                                  /* not zero : ERROR */
      if(myrank == 0) fprintf(stderr,"ERROR: (rpn_comm_shmget) processes are not all on same node \n");
      return NULL;    /* error */
      }
  }else{    // make c_comm from c_comm_in
    size = -shm_size;
    ierr = MPI_Comm_rank(c_comm_in,&myrank);
    if(myhost < 0) myhost = -myhost;                          /* color must be positive for MPI_Comm_split */
    ierr = MPI_Comm_split(c_comm_in,myhost,myrank,&c_comm);   /* split input communicator so that all members are on same SMP node*/
  }

  ierr=MPI_Comm_rank(c_comm,&myrank);
  if(myrank == 0) id=shmget(IPC_PRIVATE,size,IPC_CREAT|S_IRUSR|S_IWUSR);  /* rank 0 allocates shared memory segment */
  ierr=MPI_Bcast(&id,1,MPI_INTEGER,0,c_comm);                             /* all processes get id */
  if(id == -1) {
    if(myrank == 0) printf("ERROR: (rpn_comm_shmget) cannot create shared memory segment\n");
    return NULL;    /* error */
  }else{
    if(myrank == 0) printf("INFO: (rpn_comm_shmget) created shared memory segment of size %d\n",shm_size);
  }

  ptr=shmat(id,NULL,0);                                                   /* all processes attach memory segment */
  if(ptr == NULL) printf("ERROR: (rpn_comm_shmget) got a null pointer from shmat, process = %d\n",myrank);
  myhost = (ptr == NULL);   /* this better be zero */
//  if(ptr != NULL) printf("DEBUG: (rpn_comm_shmget) process %d attached to segment\n",myrank);
//  printf("DEBUG: addr = %16p\n",ptr);
  ierr=MPI_Allreduce(&myhost,&myhost2,1,MPI_INTEGER,MPI_BOR,c_comm); /* boolean OR from all members of this comunicator */
//  ierr=MPI_Barrier(c_comm);                                         /* all processes have attached the segment */

  if(myrank == 0) shmctl(id,IPC_RMID,&shm_buf);      /* mark segment for deletion to make sure it is released when all processes terminate */
  if(0 != myhost2){                                                  /* not zero : OUCH */
    if(myrank == 0) fprintf(stderr,"ERROR: (rpn_comm_shmget) some processes were not able to attach to segment \n");
    return NULL;    /* error */
    }
  return ptr;                                        /* return pointer tio shared memory area */
}

unsigned long long rdtscp(int *socket, int *processor);
static  int my_numa = 0;
static  int my_core = 0;

// all members of the c_comm_in communicator MUST be on the same node (same hostid)
void *C_rpn_comm_shmget_numa(MPI_Comm c_comm_in, unsigned int shm_size)  /* allocate a shared memory segment by numa space */
{
  MPI_Comm New_Comm = c_comm_in;
  int status;
  unsigned long long dummy;

#if defined(__x86_64__) &&  defined(__linux__)
  dummy = rdtscp(&my_numa, &my_core);  // on non X86 systems, per numaspace becomes per node
#endif
  // step 1, split communicator c_comm_in by numa space to produce New_Comm
  status = MPI_Comm_split(c_comm_in,my_numa,my_core,&New_Comm);
  return C_rpn_comm_shmget(New_Comm, shm_size);
}

// all members of the c_comm_in communicator need not be be on the same node (same hostid)
void *C_rpn_comm_shmget_all(MPI_Comm c_comm_in, unsigned int shm_size)  /* allocate a shared memory segment by host */
{
  MPI_Comm New_Comm = c_comm_in;
  int status;
  unsigned long long dummy;
  int hostid=gethostid();
  int localrank, ierr;

  hostid &= 0x7FFFFFFF;
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&localrank);
  status = MPI_Comm_split(c_comm_in,hostid,localrank,&New_Comm);  // split communicator c_comm_in by host to produce New_Comm
  return C_rpn_comm_shmget(New_Comm, shm_size);
}

// all members of the c_comm_in communicator need not be be on the same node (same hostid)
void *C_rpn_comm_shmget_all_numa(MPI_Comm c_comm_in, unsigned int shm_size)  /* allocate a shared memory segment by numa space */
{
  MPI_Comm New_Comm = c_comm_in;
  int status;
  unsigned long long dummy;
  int hostid=gethostid();
  int localrank, ierr;

  hostid &= 0x7FFFFFFF;
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&localrank);
  status = MPI_Comm_split(c_comm_in,hostid,localrank,&New_Comm);  // split communicator c_comm_in by host to produce New_Comm
  return C_rpn_comm_shmget_numa(New_Comm, shm_size);
}

void *F_rpn_comm_shmget(MPI_Fint f_comm, unsigned int shm_size)  /* allocate a shared memory segment */
{
//  MPI_Fint f_comm=comm;                                   /* all members of this communicator MUST be on same SMP node */
  return(C_rpn_comm_shmget(MPI_Comm_f2c(f_comm), shm_size));  /* translate Fortran communicator into C communicator before call */
}

void *F_rpn_comm_shmget_all(MPI_Fint f_comm, unsigned int shm_size)  /* allocate a shared memory segment */
{
//  MPI_Fint f_comm=comm;                                   /* all members of this communicator need not be on same SMP node */
  return(C_rpn_comm_shmget_all(MPI_Comm_f2c(f_comm), shm_size));  /* translate Fortran communicator into C communicator before call */
}

void *F_rpn_comm_shmget_numa(MPI_Fint f_comm, unsigned int shm_size)  /* allocate a shared memory segment per  numa space*/
{
//  MPI_Fint f_comm=comm;                                   /* all members of this communicator MUST be on same SMP node */
  return(C_rpn_comm_shmget_numa(MPI_Comm_f2c(f_comm), shm_size));  /* translate Fortran communicator into C communicator before call */
}

void *F_rpn_comm_shmget_all_numa(MPI_Fint f_comm, unsigned int shm_size)  /* allocate a shared memory segment per  numa space*/
{
//  MPI_Fint f_comm=comm;                                   /* all members of this communicator need not be on same SMP node */
  return(C_rpn_comm_shmget_all_numa(MPI_Comm_f2c(f_comm), shm_size));  /* translate Fortran communicator into C communicator before call */
}

#if defined(SELF_TEST)
main(int argc, char **argv){
//  int MPI_Init(int *argc, char ***argv)
  int ierr;
  int *sharedmem;
  int localrank;
  int hostid=gethostid();
  MPI_Comm comm = MPI_COMM_WORLD;

  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&localrank);
//   hostid &= 0x7FFFFFFF;
//   ierr = MPI_Comm_split(MPI_COMM_WORLD,hostid,localrank,&comm);
  sharedmem = (int *) C_rpn_comm_shmget_all_numa(comm,1024*1024);
  if(myrank==0) sharedmem[0] = localrank + 100000;
  ierr = MPI_Barrier(MPI_COMM_WORLD);
  printf("Rank = %5.5d, Address = %p, contents = %d, core=%2.2d, numa=%d, hostid=%x\n",localrank,sharedmem,sharedmem[0],my_core,my_numa,hostid);
  ierr = MPI_Finalize();
}
#endif
