#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <string.h>
#include <malloc.h>

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
int MPI_mgi_init(void);
void MPI_mgi_closeall(void);
int MPI_mgi_open(const char *channel_name);
int MPI_mgi_read(int channel, unsigned char *data, unsigned char *dtyp, int nelm);
int MPI_mgi_write(int channel, unsigned char *data, unsigned char *dtyp, int nelm);
void *MPI_mgi_memptr(int channel);
static int MPI_Can_Publish_name(const char *service_name, int test);

#if defined(WITH_MGI)
ftnword f77_name (mgi_init) (char *channel_name, F2Cl lname);
ftnword f77_name (mgi_open) (ftnword *f_chan, char *mode, F2Cl lmode);
ftnword f77_name (mgi_read) (ftnword *f_chan, void *data, ftnword *f_nelm, char *dtype, F2Cl ltype);
ftnword f77_name (mgi_write) (ftnword *f_chan, void *data, ftnword *f_nelm, char *dtype, F2Cl ltype);
ftnword f77_name (mgi_clos) (ftnword *f_chan);
ftnword f77_name (mgi_term) ();
#endif

static int in_closeall = 0 ;
static int closeall_done = 0 ;

int MPI_Finalize(){                  // interceptor, there are things we may want to do at MPI_Finalize time
//  printf("INFO: custom MPI_Finalize\n");
  while (nf >= 0) (*table[nf--])() ; // call registered routines to be called at mpi_finalize in reverse order
  return PMPI_Finalize();            // now call the real MPI finalize routine
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
  char *channel_name;     // full MPI channel name
  char *port_name;        // MPI port name (from MPI_Open_port)
  char *alias1;           // logical name 1 associated to MPI channel channel_name
  char *alias2;           // logical name 2 associated to MPI channel channel_name
  arena *winbuf;          // local one sided window memory address   (local memory arena)
  MPI_Win window;         // one sided MPI window communicator
  MPI_Comm global, local; // inter-communicator, intra-communicator, both expected to have same size membership (2 members) 
  int is_server;          // 0 for client, 1 for server
  int is_active;          // 1 if active, 0 if inactive
  int winsize;            // window size in KBytes
  int thispe;             // my rank in intra-communicator for this channel
  int otherpe;            // rank of remote process in intra-communicator for this channel
  arena remote;           // copy of memory arena control parameters from remote end
} MPI_mgi_channel ;

// typedef struct{
//   char *channel_name;     // bidirectional channel name
//   char *port_name;        // MPI port name (from MPI_Open_port)
// } MPI_port;

#define MAX_CHANNELS 8
static MPI_mgi_channel mpi_channel_table[MAX_CHANNELS];
static int last_mpi_channel=-1;

// static MPI_port mpi_port_table[MAX_CHANNELS];
// static int last_mpi_port=-1;

void *MPI_mgi_memptr(int channel){
  if(channel >= MAX_CHANNELS) return NULL;
  return mpi_channel_table[last_mpi_channel].winbuf;
}

static int timeout = 1000000;  // 1000 seconds by default
int MPI_mgi_timeout(int new_timeout){
  int old = timeout;
  timeout = new_timeout;
  return old;
}

int MPI_mgi_term() // close the books (we may prefer to wait for MPI_Finalize in order to avoid undue wait for close)
{
//   MPI_mgi_closeall() ;
}

// look into port/channel name table to see if channel name is known
// static char *MPI_mgi_get_port_name(const char *channel_name)
// {
//   char *found = NULL;
//   int i;
//   for(i=0 ; i < MAX_CHANNELS ; i++){
//     if( 0 == strncmp(mpi_port_table[i].channel_name,channel_name,MPI_MAX_PORT_NAME) ) {
//       found = mpi_port_table[i].port_name ;
//       break ;
//     }
//   }
//   return(found);
// }

// create one channel (it is assumed that no mgi activity has started)
// this channel has order number  ordinal (origin 1)
//
// MGI_MPI_CFG=" mastername n : size_1 aname_1 bname_1 : ... : size_n aname_n bname_n "
int MPI_mgi_create(const char *alias) 
{
  int status, status2 ;
  char *cfg ;
  char *service_name ;
  char mastername[128] ;
  int i, ordinal ;

  ordinal = -1 ;

  if(last_mpi_channel != -1){                 // initialize global channel table (do only once)
    cfg = getenv("MGI_MPI_CFG") ;
    if(cfg == NULL) return(-1) ;
    sscanf(cfg,"%64s%d",mastername,&last_mpi_channel) ;   // get mastername and number of mpi channels
    last_mpi_channel--;
    for(i=0 ; i<=last_mpi_channel ; i++) {    // get size, alias1, alias2 for all channels, build channel name
      while(*cfg != ':') {                    // skip to : delimiter
	if(*cfg == '\0') return(-1);          // premature termination of configuration string
	cfg++;
      }
      cfg++;
      mpi_channel_table[i].channel_name = (char *) malloc(257);
      memset(mpi_channel_table[i].channel_name,0,257);
      snprintf(mpi_channel_table[i].channel_name,257,"%s_%d",mastername,i);
      mpi_channel_table[ordinal].port_name = NULL ;

      mpi_channel_table[i].alias1 = (char *) malloc(33);          // alias1 for MPI channel
      memset(mpi_channel_table[i].alias1,0,33);
      mpi_channel_table[i].alias2 = (char *) malloc(33);          // alias2 for MPI channel
      memset(mpi_channel_table[i].alias2,0,33);
      sscanf(cfg,"%d%32s%32s",&mpi_channel_table[i].winsize,mpi_channel_table[i].alias1,mpi_channel_table[i].alias2);

      mpi_channel_table[i].winbuf = NULL ;                        // window not created yet
      mpi_channel_table[i].window = MPI_WIN_NULL ;
      mpi_channel_table[i].global = MPI_COMM_NULL ;
      mpi_channel_table[i].local  = MPI_COMM_NULL ;
      mpi_channel_table[i].is_server = 0 ;
      mpi_channel_table[i].is_active = 0 ;
      mpi_channel_table[i].winsize   = 0 ;                          // window size associated with MPI channel
      mpi_channel_table[i].thispe    = -1 ;
      mpi_channel_table[i].otherpe   = -1 ;
      if( (strcmp(alias,mpi_channel_table[i].alias1) == 0) || (strcmp(alias,mpi_channel_table[i].alias2) == 0) ) {   // alias found
	ordinal = i ;
	service_name = mpi_channel_table[i].channel_name ; // MPI channel name is 'prefix_ordinal'
      }
    }
  }

  if(ordinal < 0) return(-1) ;       // alias not found , OOPS 

  status2 = MPI_Can_Publish_name(service_name,0) ;    // is service_name already published ?
  if(status2 == -1) return(0);                        // yes
  
  mpi_channel_table[ordinal].port_name = malloc(MPI_MAX_PORT_NAME+1);
  mpi_channel_table[ordinal].is_server = 1 ;
  
  status = MPI_Open_port(MPI_INFO_NULL, mpi_channel_table[ordinal].port_name);    // create port, get port name
  if (status == MPI_SUCCESS) {                    // publish port under name channel_name
    status = MPI_mgi_Publish_name(service_name, MPI_INFO_NULL, mpi_channel_table[ordinal].port_name);
  }
  if(status != MPI_SUCCESS) {              // port creation failed
    free(mpi_channel_table[ordinal].channel_name) ;  mpi_channel_table[ordinal].channel_name = NULL ;
    free(mpi_channel_table[ordinal].alias1)    ; mpi_channel_table[ordinal].alias1 = NULL ;
    free(mpi_channel_table[ordinal].alias2)    ; mpi_channel_table[ordinal].alias2 = NULL ;
    free(mpi_channel_table[ordinal].port_name) ;  mpi_channel_table[ordinal].port_name = NULL ;
    return(-1) ; 
  }

  return(0);
}

int MPI_mgi_close(int channel) // close one channel (it is assumed that there is no pending activity)
{
  if(channel >= MAX_CHANNELS || channel < 0) return -1;             // invalid channel number
  if(mpi_channel_table[channel].channel_name == NULL) return -1;    // channel not/no longer active
  if(mpi_channel_table[channel].is_active == 0) return -1;          // channel not/no longer active

  mpi_channel_table[channel].is_active = 0;                         // mark channel as inactive
  if(in_closeall == 0) return 0;                                    // wait for closeall to really close

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
  MPI_Can_Publish_name(mpi_channel_table[channel].channel_name,1) ; // remove lock file if not already done

  free(mpi_channel_table[channel].channel_name);                    // free channel name array
  mpi_channel_table[channel].channel_name = NULL;
  free(mpi_channel_table[channel].port_name);                       // free port name array
  mpi_channel_table[channel].port_name = NULL;
  mpi_channel_table[channel].is_server = -1;                        // invalidate server flag
  mpi_channel_table[channel].thispe = -1;                           // invalidate ranks
  mpi_channel_table[channel].otherpe = -1;

  return 0;
}

void MPI_mgi_closeall(void)  // close all channels if not already done
{
  int i;
  if(closeall_done) return ;  // already done
  in_closeall = 1 ;
  for(i=0 ; i<=last_mpi_channel ; i++) MPI_mgi_close(i);
  in_closeall = 0 ;
  closeall_done = 1 ;
}

//
// secondary init, MUST BE CALLED ONLY BY A RANK 0 PROCESS
//
static int init_2_called = 0;
void MPI_mgi_init_2(void)
{
  if(init_2_called) return ; // active only once, called by MPI_mgi_open
  init_2_called = 1;
  return ;
}

//
// this function MUST BE CALLED ONLY BY A RANK 0 PROCESS
//
// MGI_MPI_CFG=" prefix n : size1 aname_1 bname_1 : ... : aname_n bname_n "
int MPI_mgi_init(void)
{
  char *port_name;
  char *cfg;
  char prefix[257];
  int i, status;

  cfg = getenv("MGI_MPI_CFG");
  if(cfg == NULL) return -1;    // no config found
  sscanf(cfg,"%256s%d",prefix,&last_mpi_channel);  // get name prefix and number of MPI channels
  last_mpi_channel--;      // origin 0
  for(i=0 ; i<=last_mpi_channel ; i++) {    // get size, alias1, alias2 for all channels, build channel name
    while(*cfg != ':') {                    // skip to : delimiter
      if(*cfg == '\0') return(-1);          // premature termination of configuration string
      cfg++;
    }
    cfg++;
    mpi_channel_table[i].channel_name = (char *) malloc(257);
    memset(mpi_channel_table[i].channel_name,0,257);
    mpi_channel_table[i].alias1 = (char *) malloc(33);
    memset(mpi_channel_table[i].alias1,0,33);
    mpi_channel_table[i].alias2 = (char *) malloc(33);
    memset(mpi_channel_table[i].alias2,0,33);
    mpi_channel_table[i].winsize = 0;
    sscanf(cfg,"%d%32s%32s",&mpi_channel_table[i].winsize,mpi_channel_table[i].alias1,mpi_channel_table[i].alias2);
    snprintf(mpi_channel_table[i].channel_name,257,"%s_%d",prefix,i);
  }
//   for(i=0 ; i<=last_mpi_channel ; i++) {
//     status = MPI_Open_port(MPI_INFO_NULL, port_name);    // create port, get port name
//     if(status != MPI_SUCCESS) return (-1);
//     MPI_mgi_Publish_name(mpi_channel_table[i].channel_name, MPI_INFO_NULL, port_name);   // publish it under name channel_name
//   }

  at_MPI_Finalize(MPI_mgi_closeall) ;// setup for at_MPI_Finalize (close everything at finalize time)

  return (0);
}

//
// this function MUST BE CALLED ONLY BY A RANK 0 PROCESS
//
// channel_name character string used to refer to the communication channel
// server       0 for client, 1 for server (must not be the same at both ends)
// window_size  is in integer word units (32 bits)
//
// result       -1 in case of failure, otherwise channel number (>-0) for read/write/close calls
int MPI_mgi_open(const char *alias)
{
  char *port_name;
  char *channel_name ;
  MPI_Win window;
  int rank;
  MPI_Comm global, local;
  MPI_Aint winsize;
  int dispunit = sizeof(int);
  void *memptr;
  arena *remote;
  MPI_Aint TargetDisp;
  int data_start ;
  arena *arenaptr ;
  int otherpe ;
  int server ;
  int mpi_channel, i ;
  int window_size ;

  mpi_channel = -1 ;
  for(i=0 ; i<last_mpi_channel ; i++){
    if( (strcmp(mpi_channel_table[i].alias1,alias) == 0) || (strcmp(mpi_channel_table[i].alias2,alias) == 0) ){
      mpi_channel = i ;
      break ;
    }
  }
  if(mpi_channel == -1) return(-1) ;     // unknown alias

  channel_name = mpi_channel_table[mpi_channel].channel_name ;   // true_name of MPI channel addresses as alias1 or alias2
  window_size  = mpi_channel_table[mpi_channel].winsize ;        // size in "int"
  winsize      = dispunit * window_size ;                        // size in bytes

  server = mpi_channel_table[mpi_channel].is_server;             // get server flag

  if(! server) mpi_channel_table[mpi_channel].port_name=malloc(MPI_MAX_PORT_NAME+1);  // no , allocate space for name from lookup
  port_name = mpi_channel_table[mpi_channel].port_name;

  if(server){
    MPI_Comm_accept(port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF,  &global );  // wait for client to connect (collective over communicator)
    MPI_Intercomm_merge(global, 0, &local);                                  // create communicator for one-sided window
  }else{
    MPI_mgi_Lookup_name(channel_name, MPI_INFO_NULL, port_name);             // get port name published under name channel_name
    MPI_Comm_connect(port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF,  &global ); // wait until connected to server (collective over communicator)
    MPI_Intercomm_merge(global, 1, &local);                                  // create communicator for one-sided window
  }
  MPI_Comm_rank(local,&mpi_channel_table[mpi_channel].thispe);               // size of communicator local should be 2
  otherpe = 1 - mpi_channel_table[mpi_channel].thispe;
  mpi_channel_table[mpi_channel].otherpe = otherpe;                          // otherpe is 0 if i am 1, 1 if i am 0
  mpi_channel_table[mpi_channel].global = global;
  mpi_channel_table[mpi_channel].local  = local;
  mpi_channel_table[mpi_channel].is_active = 1;

  MPI_Alloc_mem(winsize,MPI_INFO_NULL,&memptr);                              // allocate local memory for one-sided window
  mpi_channel_table[mpi_channel].winbuf = memptr;                            // local memory address of 1 sided window
  MPI_Win_create(memptr, winsize, dispunit, MPI_INFO_NULL, local, &window);  // create one-sided window (local buffer to receive remote writes)
  mpi_channel_table[mpi_channel].window = window;

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
  remote = &(mpi_channel_table[mpi_channel].remote) ;
  MPI_Get(remote, data_start, MPI_INTEGER, otherpe, TargetDisp, data_start, MPI_INTEGER, window);
  MPI_Win_unlock(otherpe,window);
  MPI_Barrier(local);    // sync with remote process

  return mpi_channel ;                                                  // return mpi channel number;
}

// read nelm data elements of type *dtyp (I/R/D/C) into data from channel
int MPI_mgi_read(int channel,unsigned char *data, unsigned char *dtyp, int nelm){
  MPI_Datatype mpitype;
  int nitems = nelm;
  int ntok;
  int avail = 0;
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
    case 'D':     // 64 bit floating point data (we might want to downgrade on write, upgrade on read)
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
  memory->out = out ;                                      // update out
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
    case 'D':     // 64 bit floating point data (we might want to downgrade on write, upgrade on read)
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
    remote->out = out;              // update local copy of remote out
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
    MPI_Put(idata, ntok, MPI_INTEGER, otherpe, TargetDisp, ntok, MPI_INTEGER, window);     // remote write data
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
  unlink(filename);                                  // remove channel file
  MPI_Can_Publish_name(service_name,1) ;             // remove lock file

  return MPI_SUCCESS;
}

static int MPI_Can_Publish_name(const char *service_name, int test)
{
  FILE *gossip;
  char filename[4096];
  char filenew[4096];
  int status;

  snprintf(filename,4096,"%s/%s",getenv("HOME"),".gossip");
  mkdir(filename,0755);
  snprintf(filename,4096,"%s/%s",getenv("HOME"),".gossip/MPI");   // make sure that ~/.gossip/MPI directory exists
  mkdir(filename,0755);

  snprintf(filename,4096,"%s/%s/%s.lock",getenv("HOME"),".gossip/MPI",service_name);
  if(test == 0) {
    status = creat(filename,00700);   // test mode, try to create lock file
  }else{
    status = unlink(filename);        // cleanup mode, remove lock file
  }

  return status;
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
  printf("client rank %d of %d\n",rank+1,size);
#endif
  sleep(2);
  MPI_Finalize();
  return 0;
}
#endif
