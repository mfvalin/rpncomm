#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <string.h>

// the following 3 routines are temporarily static and will have to be moved to their own source file
static int MPI_mgi_Unpublish_name(const char *service_name, MPI_Info info, const char *port_name);
static int MPI_mgi_Publish_name(const char *service_name, MPI_Info info, const char *port_name);
static int MPI_mgi_Lookup_name(const char *service_name, MPI_Info info, char *port_name);

// the following 2 routines and associated declarations will have to be moved to their own source file
typedef void (*fptr)(void);
#define MAX_AT_TABLE 10
static fptr table[MAX_AT_TABLE];
static int nf=-1;
int MPI_Finalize();
int at_MPI_Finalize(fptr callback);

/*
 *     data buffer layout
 * 
 *     first : points to first element of buffer            (never updated)
 *     limit : points one past the last element in buffer   (never updated)
 *     in    : points to insertion position                 (updated by process that inserts data)
 *     out   : points to extraction position                (updated by process that extracts data)
 * 
 *     empty buffer, in = out
 * 
 *     +-----------------------------------------------------------------------------------------+
 *     | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | +
 *     +-----------------------------------------------------------------------------------------+
 *      ^                   ^                                                                     ^
 *      | first          in | out                                                                 | limit
 *     ( 0 data elements in buffer, room for (limit - first - 1)
 * 
 *     buffer neither empty nor full
 * 
 *     +-----------------------------------------------------------------------------------------+
 *     | | | | | | | | | | |d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d| | | | | | | | | | | | | | | | +
 *     +-----------------------------------------------------------------------------------------+
 *      ^                   ^                                     ^                               ^
 *      | first             | out                                 | in                            | limit
 *     (in - out) data elements in buffer, room for (limit - in) + (out - in -1) more
 * 
 * 
 *     +-----------------------------------------------------------------------------------------+
 *     |d|d|d|d|d|d|d|d|d|d| | | | | | | | | | | | | | | | | | | |d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d+
 *     +-----------------------------------------------------------------------------------------+
 *      ^                   ^                                     ^                               ^
 *      | first             | in                                  | out                           | limit
 *     (in - first) + limit - out) data elements in buffer, room for (out - in + 1)
 * 
 * 
 *     full buffer, in = (out -1) (modulo limit)
 *     +-----------------------------------------------------------------------------------------+
 *     |d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d| |d|d|d|d|d|d|d|d|d|d|d|d|d|d|d|d+
 *     +-----------------------------------------------------------------------------------------+
 *      ^                                                       ^ ^                               ^
 *      | first                                              in | | out                           | limit
 */
int MPI_mgi_close(int channel);
void MPI_mgi_closeall(void);
int MPI_mgi_open(const char *channel_name,int server, int window_size);
int MPI_mgi_read(int channel, unsigned char *data, unsigned char *dtyp, int nelm);
int MPI_mgi_write(int channel, unsigned char *data, unsigned char *dtyp, int nelm);

static int in_closeall = 0;

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
  int control;    // flags
  int first;      // start of buffer
  int in;         // insertion index
  int out;        // extraction index
  int limit;      // end of buffer + 1
  int data_start[1];
} arena;

typedef struct{
  char *channel_name;     // character string used to refer to channel by name
  char *port_name;        // MPI port name (from MPI_Open_port)
  arena *winbuf;          // local one sided window memory address   (local memory arena)
  MPI_Win window;         // one sided MPI window communicator
  MPI_Comm global, local; // inter-communicator, intra-communicator, both expected to have same size membership (2 members) 
  int is_server;          // 0 for client, 1 for server
  int is_active;          // 1 if active, 0 if inactive
  int thispe;             // my rank in intra-communicator for this channel
  int otherpe;            // rank of remote process in intra-communicator for this channel
  arena remote;           // copy of memory arena control parameters from remote end
} MPI_mgi_channel ;

#define MAX_CHANNELS 8
static MPI_mgi_channel mpi_channel_table[MAX_CHANNELS];
static int last_mpi_channel=-1;

void *MPI_mgi_channel_mem(int channel){
  if(channel >= MAX_CHANNELS) return NULL;
  return mpi_channel_table[last_mpi_channel].winbuf;
}

static int timeout = 1000000;  // 1000 seconds by default
int MPI_mgi_timeout(int new_timeout){
  int old = timeout;
  timeout = new_timeout;
  return old;
}

int MPI_mgi_close(int channel) // close one channel (it is assumed that there is no pending activity)
{
  if(channel >= MAX_CHANNELS || channel < 0) return -1;             // invalid channel number
  if(mpi_channel_table[channel].channel_name == NULL) return -1;    // channel not/no longer active
  if(mpi_channel_table[channel].is_active == 0) return -1;          // channel not/no longer active

  mpi_channel_table[channel].is_active = 0;
  if(in_closeall == 0) return 0;
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
  mpi_channel_table[channel].thispe = -1;
  mpi_channel_table[channel].otherpe = -1;

  return 0;
}

void MPI_mgi_closeall(void)  // close all channels
{
  int i;
  in_closeall = 1;
  for(i=0 ; i<=last_mpi_channel ; i++) MPI_mgi_close(i);
}

//
// this function MUST BE CALLED ONLY BY A RANK 0 PROCESS
//
// channel_name character string used to refer to the communication channel
// server       0 for client, 1 for server (must not be the same at both ends)
// window_size  is in integer word units (32 bits)
//
// result       -1 in case of failure, otherwise channel number (>-0) for read/write/close calls
int MPI_mgi_open(const char *channel_name,int server, int window_size)
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
  int otherpe;

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
  MPI_Comm_rank(local,&mpi_channel_table[last_mpi_channel].thispe);
  otherpe = 1 - mpi_channel_table[last_mpi_channel].thispe;
  mpi_channel_table[last_mpi_channel].otherpe = otherpe;
  mpi_channel_table[last_mpi_channel].global = global;
  mpi_channel_table[last_mpi_channel].local = local;
  mpi_channel_table[last_mpi_channel].is_active = 1;

  winsize = dispunit * winsize;
  MPI_Alloc_mem(winsize,MPI_INFO_NULL,&memptr);                              // allocate memory for one-sided window
  mpi_channel_table[last_mpi_channel].winbuf = memptr;
  MPI_Win_create(memptr, winsize, dispunit, MPI_INFO_NULL, local, &window);  // create one-sided window (local buffer to receive remote writes)
  mpi_channel_table[last_mpi_channel].window = window;
  arenaptr = (arena *)memptr ;                                               // local one-sided memory area
// setup of local memory arena
  data_start = &(arenaptr->data_start[0]) - &(arenaptr->control);   // offset of first data element
  arenaptr->control = CHANNEL_ACTIVE;
  arenaptr->first = data_start;                 // will use memptr[first] to memptr[limit-1]
  arenaptr->in = arenaptr->first;               // in = out = first, buffer starts empty
  arenaptr->out = arenaptr->in;
  arenaptr->limit = window_size - data_start - 1;
  MPI_Barrier(local);    // sync with remote process
// get remote memory arena configuration from "otherpe" (later we will "put" remote in and "get" remote out)
  MPI_Win_lock(MPI_LOCK_SHARED,otherpe,0,window);
  TargetDisp = 0;        //  what we get is right at the beginning of memory area, length = data_start integers
  MPI_Get(memptr, data_start, MPI_INTEGER, otherpe, TargetDisp, data_start, MPI_INTEGER, window);
  MPI_Win_unlock(otherpe,window);
  MPI_Barrier(local);    // sync with remote process

  return last_mpi_channel ;                                                  // return mpi channel number;
}

// read nelm data elements of type *dtyp (I/R/D/C) into data from channel
int MPI_mgi_read(int channel, unsigned char *data, unsigned char *dtyp, int nelm){
  MPI_Datatype mpitype;
  int nitems = nelm;
  int ntok;
  int avail = 0;
  arena *remote;
  arena *memory;
  int *buffer;
  useconds_t wait = 1000;   // 1 millisecond
  int sleepcount = 0;
  int meta;
  MPI_mgi_channel *local;
  int first, in, out, limit, navail;
  size_t nbytes, nbytes1, nbytes2;
  unsigned char rtyp = '?';
  int nmeta;

  if(channel >= MAX_CHANNELS || channel < 0) return -1;             // invalid channel number
  if(mpi_channel_table[channel].channel_name == NULL) return -1;    // channel not/no longer active
  if(mpi_channel_table[channel].is_active == 0) return -1;          // channel not/no longer active

  switch(*dtyp){    // tokens are 32 bit integers
    case 'R':       // 32 bit floating point data
      ntok = nitems;
      nbytes = nitems << 2;   // nitems * 4
      mpitype = MPI_FLOAT;
      break;
    case 'I':      // 32 bit integer data
      ntok = nitems;
      nbytes = nitems << 2;   // nitems * 4
      mpitype = MPI_INT;
      break;
    case 'C':      // character data
      ntok = (nitems + 3) >> 2;  //  (nitems rounded up to a multiple of 4) / 4
      nbytes = nitems;
      mpitype = MPI_CHAR;
      break;
    case 'D':     // 64 bit floating point data
      ntok = nitems << 1;     // nitems * 2
      nbytes = nitems << 3;   // nitems * 8
      mpitype = MPI_DOUBLE;
      break;
    default :     // anything else is an error
      fprintf(stderr,"ERROR: MPI_mgi_read, bad type '%c', valid types are I/R/D/C\n",*dtyp);
      return -1;
      break;
  }
  memory = (arena *)mpi_channel_table[channel].winbuf;  // local memory arena
  buffer = (int *) memory;
  remote = &(mpi_channel_table[channel].remote);

// we need 1 metadata element + ntok data elements
  while(1){                      // wait until there is enough data in local buffer to satisfy request
    first = memory->first;       // index of start of data buffer
    in    = memory->in;          // will be updated by remote writer process
    out   = memory->out;         // will be updated by this process
    limit = memory->limit;       // bufer[limit] is one slot past end of buffer
    navail = (in >= out) ? (in - out) : ((limit - out) + (in - first));  // number of available data elements
    if(navail >= ntok + 1) break; // enough data elements to satisfy request
    usleep(wait);                 // not enough, waiting for remote writer to update "in"
    sleepcount++;
    // timeout ?
    if(sleepcount > timeout) {
      return 0;
    }
  }   // local data buffer : buffer[first] to buffer[limit-1];

  if(out >= limit) out = first;  // wrap around
  meta = buffer[out];            // get metadata
  // validate metadata ( type + 'length' << 8 )   'length' = real length & 0xFFFFFF
  nmeta = meta >> 8;
  rtyp  = meta && 0xFF;
  if(rtyp != *dtyp || nmeta != (nelm & 0xFFFFFF)) {
    return -1;
    fprintf(stderr,"ERROR: MPI_mgi_read, type/length mismatch, expected %d'%c', got %d,'%c'\n",nelm,*dtyp,nmeta,rtyp);
  }
  out++;
  if(out >= limit) out = first;

  if(in > out){                                            // data is in 1 piece
    memcpy(data,&(buffer[out]),nbytes);                    // buffer[out] -> buffer[out + ntok -1] : ntok tokens
    out += ntok;                                           // got data, update out
  }else{                                                   // data is split in 2 pieces
    nbytes1 = (limit - out) << 2;                          // (limit - out) * 4 bytes for 1st piece
    memcpy(data,&(buffer[out]),nbytes1);                   // buffer[out] -> buffer[limit-1] : limit - out tokens
    nbytes2 = nbytes - nbytes1;                            // nbytes = 1st piece
    memcpy(&(data[nbytes1]),&(buffer[first]),nbytes2);     // buffer[first] -> buffer[out-1] : ntok - (limit - out) tokens
    out = first + (ntok - (limit - out));                  // data has been copied, it is safe to update out
  }
  return nelm;     // return number of elements read
}

// write nelm data elements of type *dtyp (I/R/D/C) from data into channel 
int MPI_mgi_write(int channel, unsigned char *data, unsigned char *dtyp, int nelm){
  MPI_Datatype mpitype;
  int nitems = nelm;
  int ntok;
  int avail = 0;
  int first, in, out, limit, navail;
  size_t nbytes, nbytes1, nbytes2;
  int sleepcount = 0;
  int wait = 1000;
  arena *remote;
  MPI_Win window;
  MPI_Aint TargetDisp = 0;
  int meta;
  int *idata = (int *)data;
  int otherpe;

  if(channel >= MAX_CHANNELS || channel < 0) return -1;             // invalid channel number
  if(mpi_channel_table[channel].channel_name == NULL) return -1;    // channel not/no longer active
  if(mpi_channel_table[channel].is_active == 0) return -1;          // channel not/no longer active

  switch(*dtyp){    // tokens are 32 bit integers
    case 'R':       // 32 bit floating point data
      ntok = nitems;
      nbytes = nitems << 2;   // nitems * 4
      mpitype = MPI_FLOAT;
      break;
    case 'I':      // 32 bit integer data
      ntok = nitems;
      nbytes = nitems << 2;   // nitems * 4
      mpitype = MPI_INT;
      break;
    case 'C':      // character data
      ntok = (nitems + 3) >> 2;  //  (nitems rounded up to a multiple of 4) / 4
      nbytes = nitems;
      mpitype = MPI_CHAR;
      break;
    case 'D':     // 64 bit floating point data
      ntok = nitems << 1;     // nitems * 2
      nbytes = nitems << 3;   // nitems * 8
      mpitype = MPI_DOUBLE;
      break;
    default :     // anything else is an error
      fprintf(stderr,"ERROR: MPI_mgi_read, bad type '%c', valid types are I/R/D/C\n",*dtyp);
      return -1;
      break;
  }
// may have to get remote out to update out if not enough space
  remote = &(mpi_channel_table[channel].remote);
  otherpe = mpi_channel_table[channel].otherpe;

  while(1){      // wait until there is enough data in local buffer to satisfy request (ntok + 1) data tokens
    first = remote->first;
    in = remote->in;
    out = remote->out;
    limit = remote->limit;
    navail = (in <= out) ? (out - in) : ((limit - in) + (out - first - 1));  // number of available tokens
    if(navail >= ntok + 1) break;  // enough space is available

    TargetDisp = &(remote->out) - &(remote->control);  // not enough space, get remote 'out' to see if there is more
    MPI_Win_lock(MPI_LOCK_SHARED,otherpe,0,window);
    MPI_Get(&out, 1, MPI_INTEGER, otherpe, TargetDisp, 1, MPI_INTEGER, window);   // get 'out' from remote server
    MPI_Win_unlock(otherpe,window);
    navail = (in <= out) ? (out - in) : ((limit - in) + (out - first - 1));  // number of available tokens
    if(navail >= ntok + 1) break;  // enough space is available

    usleep(wait);                  // not enough space, remote reader needs to further increment increment 'out'
    sleepcount++;
    // timeout ?
    if(sleepcount > timeout) {
      return ntok;
    }
  }
  if(in >= limit) in = first;    // wrap around

  meta = (nelm && 0xFFFFFF) + *dtyp;   // build and send 32 bit token of metadata header
  TargetDisp = in;
  MPI_Win_lock(MPI_LOCK_SHARED,otherpe,0,window);
  MPI_Put(&meta,1,MPI_INTEGER,otherpe,TargetDisp,1,MPI_INTEGER,window);      // remote write header
  MPI_Win_unlock(otherpe,window);
  in++;

  MPI_Win_lock(MPI_LOCK_SHARED,otherpe,0,window);
  TargetDisp = in;
  if(out > in){
    MPI_Put(data, ntok, MPI_INTEGER, otherpe, TargetDisp, ntok, MPI_INTEGER, window);     // remote write data
    in += ntok;
  }else{
    MPI_Put(idata, limit-in, MPI_INTEGER, otherpe, TargetDisp, limit-in, MPI_INTEGER, window);   // remote write data (part 1)
    TargetDisp = first;
    MPI_Put(&(idata[limit-in]), ntok-(limit-in), MPI_INTEGER, otherpe, TargetDisp, ntok-(limit-in), MPI_INTEGER, window);  // remote write data (part 2)
    in = first + (ntok-(limit-in));
  }
  MPI_Win_unlock(otherpe,window);

  TargetDisp = &(remote->in) - &(remote->control);
  MPI_Win_lock(MPI_LOCK_SHARED,otherpe,0,window);
  MPI_Put(&in,1,MPI_INTEGER,otherpe,TargetDisp,1,MPI_INTEGER,window);   // update in on remote server
  MPI_Win_unlock(otherpe,window);
  return 0;
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
