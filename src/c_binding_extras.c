#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <limits.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <dlfcn.h>
#include <sched.h>

#define SIZEOF(a) s_sizeof[sizeof(a)]

int main(){
  char *s_sizeof[9] = {
    "INVALID",
    "C_INT8_T",
    "C_INT16_T",
    "INVALID",
    "C_INT32_T",
    "INVALID",
    "INVALID",
    "INVALID",
    "C_INT64_T" };

  printf("integer, parameter :: C_SSIZE_T = %s\n",SIZEOF(ssize_t));
  printf("integer, parameter :: C_TIME_T = %s\n",SIZEOF(time_t));
  printf("integer, parameter :: C_SUSECONDS_T = %s\n",SIZEOF(suseconds_t));
  printf("integer, parameter :: C_USECONDS_T = %s\n",SIZEOF(useconds_t));
  printf("integer, parameter :: C_CLOCK_T = %s\n",SIZEOF(clock_t));
  printf("integer, parameter :: C_KEY_T = %s\n",SIZEOF(key_t));
  printf("real,    parameter :: C_CLOCKS_PER_SEC = %ld.0\n",CLOCKS_PER_SEC);
  printf("integer, parameter :: C_MODE_T = %s\n",SIZEOF(mode_t));
  printf("integer, parameter :: C_OFF_T  = %s\n",SIZEOF(off_t));
  printf("integer, parameter :: O_APPEND = %d\n",O_APPEND);
  printf("integer, parameter :: O_CREAT  = %d\n",O_CREAT);
  printf("integer, parameter :: O_EXCL   = %d\n",O_EXCL);
  printf("integer, parameter :: O_TRUNC  = %d\n",O_TRUNC);
  printf("integer, parameter :: O_RDWR   = %d\n",O_RDWR);
  printf("integer, parameter :: O_RDONLY = %d\n",O_RDONLY);
  printf("integer, parameter :: O_WRONLY = %d\n",O_WRONLY);
  printf("integer, parameter :: SEEK_SET = %d\n",SEEK_SET);
  printf("integer, parameter :: SEEK_CUR = %d\n",SEEK_CUR);
  printf("integer, parameter :: SEEK_END = %d\n",SEEK_END);
  printf("integer, parameter :: PATH_MAX = %d\n",PATH_MAX);
  printf("integer, parameter :: IPC_PRIVATE = %d\n",IPC_PRIVATE);
  printf("integer, parameter :: IPC_CREAT = %d\n",IPC_CREAT);
  printf("integer, parameter :: IPC_EXCL = %d\n",IPC_EXCL);
  printf("integer, parameter :: IPC_RMID = %d\n",IPC_RMID);
  printf("integer, parameter :: SIZEOF_SHMID_DS = %ld\n",sizeof(struct shmid_ds));
  printf("integer, parameter :: RTLD_LAZY = %d\n",RTLD_LAZY);
  printf("integer, parameter :: RTLD_NOW = %d\n",RTLD_NOW);
  printf("integer, parameter :: RTLD_GLOBAL = %d\n",RTLD_GLOBAL);
  printf("integer, parameter :: RTLD_LOCAL = %d\n",RTLD_LOCAL);
  printf("type, bind(C) :: CPU_SET_T\n  integer(C_INT32_T), dimension(0:%ld) :: set\nend type\n",sizeof(cpu_set_t)/sizeof(int32_t)-1);
  return 0;
}
