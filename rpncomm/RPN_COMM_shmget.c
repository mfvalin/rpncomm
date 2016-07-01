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
          integer(C_INT), intent(IN), value :: comm           !InTf!   ! RPN_COMM communicator
          integer(C_INT), intent(IN), value :: size           !InTf!   ! size in bytes of shared memory area
          type(C_PTR) :: where                                !InTf!   ! pointer to shared memory area
        end function rpn_comm_shmget                          !InTf!
#endif

void *F_rpn_comm_shmget(int comm, unsigned int shm_size)  /* allocate a shared memory segment ( < 2 GBytes ) */
{
  size_t size=shm_size;                                   /* size of shared memory segment in bytes */
  int id;
  struct shmid_ds shm_buf;
  void *ptr;
  int ierr, myrank, myhost, myhost2;
  MPI_Fint f_comm=comm;                                   /* all members of this communicator MUST be on same node */
  MPI_Comm c_comm;
  

  c_comm = MPI_Comm_f2c(f_comm);  /* translate Fortran communicator into C communicator */

  myhost=gethostid();
  ierr=MPI_Allreduce(&myhost,&myhost2,1,MPI_INTEGER,MPI_BOR,c_comm); /* boolean OR from all members of this comunicator */
  myhost ^= myhost2;                                                 /* should be zero everywhere */
  ierr=MPI_Allreduce(&myhost,&myhost2,1,MPI_INTEGER,MPI_BOR,c_comm); /* boolean OR from all members of this comunicator */
  if(0 != myhost2){                                                  /* not zero : ERROR */
    if(myrank == 0) fprintf(stderr,"ERROR: (rpn_comm_shmget) processes are not all on same node \n");
    return NULL;    /* error */
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

//  printf("id = %d %x \n",id,id);
  ptr=shmat(id,NULL,0);                                                   /* all processes attach memory segment */
  if(ptr == NULL) printf("ERROR: (rpn_comm_shmget) got a null pointer from shmat, process = %d\n",myrank);
  myhost = (ptr == NULL);   /* this better be zero */
//  if(ptr != NULL) printf("INFO: (rpn_comm_shmget) process %d attached to segment\n",myrank);
  ierr=MPI_Allreduce(&myhost,&myhost2,1,MPI_INTEGER,MPI_BOR,c_comm); /* boolean OR from all members of this comunicator */
//  ierr=MPI_Barrier(c_comm);                                         /* all processes have attached the segment */
//  printf("addr = %16p\n",ptr);

  if(myrank == 0) shmctl(id,IPC_RMID,&shm_buf);      /* mark segment for deletion to make sure it is released when all processes terminate */
  if(0 != myhost2){                                                  /* not zero : OUCH */
    if(myrank == 0) fprintf(stderr,"ERROR: (rpn_comm_shmget) some processes were not able to attach to segment \n");
    return NULL;    /* error */
    }
  return ptr;                                        /* return pointer tio shared memory area */
}

