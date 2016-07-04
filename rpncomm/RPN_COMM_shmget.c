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
#endif

void *C_rpn_comm_shmget(MPI_Comm c_comm_in, unsigned int shm_size)  /* allocate a shared memory segment */
{
 MPI_Comm c_comm;
  size_t size;                                   /* size of shared memory segment in bytes */
  int id;
  struct shmid_ds shm_buf;
  void *ptr;
  int ierr, myrank, myhost, myhost2;
  
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

void *F_rpn_comm_shmget(MPI_Fint f_comm, unsigned int shm_size)  /* allocate a shared memory segment */
{
//  MPI_Fint f_comm=comm;                                   /* all members of this communicator MUST be on same SMP node */
  return(C_rpn_comm_shmget(MPI_Comm_f2c(f_comm), shm_size));  /* translate Fortran communicator into C communicator before call */
}

