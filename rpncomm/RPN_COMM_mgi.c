#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <string.h>

static int MPI_mgi_Unpublish_name(const char *service_name, MPI_Info info, const char *port_name);
static int MPI_mgi_Publish_name(const char *service_name, MPI_Info info, const char *port_name);
static int MPI_mgi_Lookup_name(const char *service_name, MPI_Info info, char *port_name);

typedef void (*fptr)(void);
#define MAX_AT_TABLE 10
static fptr table[MAX_AT_TABLE];
static int nf=-1;

int MPI_Finalize(){
//  printf("INFO: custom MPI_Finalize\n");
  while (nf >= 0) (*table[nf--])() ; // call registered routines to be called at mpi_finalize in reverse order
  return PMPI_Finalize();
}
int at_MPI_Finalize(fptr callback){
  if(nf >= MAX_AT_TABLE-1) return -1;
  table[++nf] = callback;
  return(0);
}

#define CHANNEL_ACTIVE 1
typedef struct{
  int control;
  int first;
  int in;
  int out;
  int limit;
  int data_start;
} arena;

typedef struct{
  MPI_Win window;         // one sided MPI window communicator
  MPI_Comm global, local; // inter-communicator, intra-communicator, both have same size (2 members) 
  char *port_name;        // MPI port name (from MPI_Open_port)
  char *channel_name;     // character string used to refer to channel
  void *winbuf;           // one sided window memory address
  int is_server;          // 0 for client, 1 for server
  arena remote;           // copy of memory arena control parameters from other end
} MPI_mgi_channel ;

#define MAX_CHANNELS 8
static MPI_mgi_channel mpi_channel_table[MAX_CHANNELS];
static int last_mpi_channel=-1;

void *MPI_mgi_channel_mem(int channel){
  if(channel >= MAX_CHANNELS) return NULL;
  return mpi_channel_table[last_mpi_channel].winbuf;
}

int MPI_close_mgi_channel(int channel) // close one channel (it is assumed that there is no pending activity)
{
  if(channel >= MAX_CHANNELS || channel < 0) return -1;             // invalid channel number
  if(mpi_channel_table[channel].channel_name == NULL) return -1;    // channel not/no longer active

//   MPI_Barrier(mpi_channel_table[channel].local);

  MPI_Win_free(&mpi_channel_table[channel].window);                 // free window communicator
  free(mpi_channel_table[channel].winbuf);                          // free memory associated with window
  mpi_channel_table[channel].winbuf = NULL;
//   MPI_Barrier(mpi_channel_table[channel].local);
  MPI_Comm_disconnect( &mpi_channel_table[channel].local );         // disconnect intra-communicator
  MPI_Comm_disconnect( &mpi_channel_table[channel].global );        // disconnect iner-communicator
  if(mpi_channel_table[channel].is_server) {                        // on server, unpublish and close port
    MPI_mgi_Unpublish_name(mpi_channel_table[channel].channel_name, MPI_INFO_NULL, "no_port_name");
    MPI_Close_port(mpi_channel_table[channel].port_name);
  }
  free(mpi_channel_table[channel].channel_name);                    // free channel name array
  mpi_channel_table[channel].channel_name = NULL;
  free(mpi_channel_table[channel].port_name);                       // free port name array
  mpi_channel_table[channel].port_name = NULL;
  mpi_channel_table[channel].is_server = -1;                        // invalidate server flag

  return 0;
}

void MPI_mgi_closeall(void)  // close all channels
{
  int i;
  for(i=0 ; i<=last_mpi_channel ; i++) MPI_close_mgi_channel(i);
}

//
// this function MUST BE CALLED ONLY BY A RANK 0 PROCESS
//
// channel_name character string used to refer to the communication channel
// server       0 for client, 1 for server (must not be the same at both ends)
// window_size  is in integer word units (32 bits)
int MPI_open_mgi_channel(const char *channel_name,int server, int window_size)
{
  char *port_name;
  MPI_Win window;
  int rank;
  MPI_Comm global, local;
  MPI_Aint winsize = window_size;
  int dispunit = sizeof(int);
  void *memptr;
  MPI_Aint TargetDisp = 0;
  int data_start ;
  arena *arenaptr;

  if(last_mpi_channel >= MAX_CHANNELS-1) return -1;  // control table is full
  if(window_size < 1024) return -1;                  // window obviously too small

  last_mpi_channel++;           // next channel
  if(last_mpi_channel == 0) {   // first time around
    at_MPI_Finalize(MPI_mgi_closeall) ;// setup for at_MPI_Finalize
  }

  mpi_channel_table[last_mpi_channel].channel_name=malloc(strlen(channel_name)+1);
  strncpy(mpi_channel_table[last_mpi_channel].channel_name,channel_name,strlen(channel_name)+1);
  
  mpi_channel_table[last_mpi_channel].port_name=malloc(MPI_MAX_PORT_NAME+1);
  port_name = mpi_channel_table[last_mpi_channel].port_name;

  mpi_channel_table[last_mpi_channel].is_server = server;                    // keep server flag
  if(server){
    MPI_Open_port(MPI_INFO_NULL, port_name);                                 // create port, get port name
    MPI_mgi_Publish_name(channel_name, MPI_INFO_NULL, port_name);            // publish it under name channel_name
    MPI_Comm_accept(port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF,  &global );  // wait for client to connect
    MPI_Intercomm_merge(global, 0, &local);                                  // create communicator for one-sided window
  }else{
    MPI_mgi_Lookup_name(channel_name, MPI_INFO_NULL, port_name);                 // get port name published under name channel_name
    MPI_Comm_connect(port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF,  &global ); // wait until connected to server
    MPI_Intercomm_merge(global, 1, &local);                                  // create communicator for one-sided window
  }
  mpi_channel_table[last_mpi_channel].global = global;
  mpi_channel_table[last_mpi_channel].local = local;

  winsize = dispunit * winsize;
  MPI_Alloc_mem(winsize,MPI_INFO_NULL,&memptr);                              // allocate memory for one-sided window
  mpi_channel_table[last_mpi_channel].winbuf = memptr;
  MPI_Win_create(memptr, winsize, dispunit, MPI_INFO_NULL, local, &window);  // create one-sided window (local buffer to receive remote writes)
  mpi_channel_table[last_mpi_channel].window = window;
  arenaptr = (arena *)memptr ;                                               // local one-sided memory area
// setup of local arena
  data_start = &(arenaptr->data_start) - &(arenaptr->control);
  arenaptr->control = CHANNEL_ACTIVE;
  arenaptr->first = 0;
  arenaptr->in = arenaptr->first;
  arenaptr->out = arenaptr->in;
  arenaptr->limit = window_size - data_start - 1;
  MPI_Barrier(local);    // sync
// get remote arena parameters (later we will "put" remote in and "get" remote out)
  MPI_Get((int *)memptr,data_start,MPI_INTEGER,1-server,TargetDisp,data_start,MPI_INTEGER,window);
  MPI_Barrier(local);    // sync

  return last_mpi_channel ;                                                  // return mpi channel number;
}

static int MPI_mgi_Unpublish_name(const char *service_name, MPI_Info info, const char *port_name)
{
  char filename[4096];

  snprintf(filename,4096,"%s/%s/%s.channel",getenv("HOME"),".gossip/MPI",service_name);
  unlink(filename);
  return MPI_SUCCESS;
}

static int MPI_mgi_Publish_name(const char *service_name, MPI_Info info, const char *port_name)
{
  FILE *gossip;
  char filename[4096];
  char filenew[4096];

  snprintf(filename,4096,"%s/%s",getenv("HOME"),".gossip");
  mkdir(filename,0755);
  snprintf(filename,4096,"%s/%s",getenv("HOME"),".gossip/MPI");
  mkdir(filename,0755);
  snprintf(filename,4096,"%s/%s/%s.new",getenv("HOME"),".gossip/MPI",service_name);
  unlink(filename);
  snprintf(filenew,4096,"%s/%s/%s.channel",getenv("HOME"),".gossip/MPI",service_name);
  unlink(filenew);

  gossip = fopen(filename,"w");
  fprintf(gossip,"%s",port_name);
  fclose(gossip);

  link(filename,filenew);
  unlink(filename);
  return MPI_SUCCESS;
}

static int MPI_mgi_Lookup_name(const char *service_name, MPI_Info info, char *port_name)
{
  char filename[4096];
  int fd, nc;
  int wait=0;

  snprintf(filename,4096,"%s/%s/%s.channel",getenv("HOME"),".gossip/MPI",service_name);

  while( (fd=open(filename,0)) < 0) { wait++ ; usleep(1000); }
  nc=read(fd,port_name,MPI_MAX_PORT_NAME);
  close(fd);
  port_name[nc]='\0';
  printf("MPI_Lookup_name: wait time = %d msec\n",wait);

  return MPI_SUCCESS;
}

#if defined(SELF_TEST)
int main(int argc, char **argv){
  arena test_arena;
  arena *arenaptr = &test_arena;
  int size, rank;
  int data_offset =  &(arenaptr->data_start) - &(arenaptr->control);
  printf("data offset = %d\n",data_offset);
  MPI_Init( &argc, &argv );
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#if defined(SERVER)
  printf("server rank %d of %d\n",rank+1,size);
#else
  printf("server rank %d of %d\n",rank+1,size);
#endif
  sleep(2);
  MPI_Finalize();
  return 0;
}
#endif
