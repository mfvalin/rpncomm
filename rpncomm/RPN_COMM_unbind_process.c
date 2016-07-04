/* RPN_COMM - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2012  Division de Recherche en Prevision Numerique
 *                          Environnement Canada
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */
#if defined(linux)
#define _GNU_SOURCE
#include <stdint.h>
#endif
#include <sched.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>

/*
   unset any processor affinity for the current process
   if the FULL_UNBIND environment variable is defined
   ( FULL_UNBIND would normally be set by r.mpirun / ord_soumet )
   in the case of an MPI launch, only process 0 will print a message
   this code is only applicable to Linux for the time being
   and has been written to counter mpi implementation behaviour in some cases.
   (written specifically for use in the RPN_COMM library)

   this code is Linux only

   Michel Valin , 2011 / 11 / 02 UQAM (openmpi)
                  2011 / 11 / 20 UQAM (mpich / slurm)
*/

/* use pragma weak to create alternate FORTRAN entry points */
/* WARNING: some old versions of gcc do not generate the weak entry points correctly */

#if defined(linux)
static cpu_set_t set;
static int lo_core=-1;
static int hi_core=-1;
#endif

#if defined(linux)
static void X86_cpuid(uint32_t eax, uint32_t ecx, uint32_t* regs)  /* interface to x86 cpuid instruction */
{
#if defined(__x86_64__) || defined( __i386__ )
    uint32_t ebx, edx;
# if defined( __i386__ ) && defined ( __PIC__ )
     /* PIC under 32-bit EBX must not be clobbered */
    __asm__ ( "movl %%ebx, %%edi \n\t cpuid \n\t xchgl %%ebx, %%edi" : "=D" (ebx),
# else
   ebx = 0;
    __asm__ ( "cpuid" : "+b" (ebx),
# endif
              "+a" (eax), "+c" (ecx), "=d" (edx) );
    regs[0] = eax; regs[1] = ebx; regs[2] = ecx; regs[3] = edx;
#else
    regs[0] = 0; regs[1] = 0; regs[2] = 0; regs[3] = 0;
#endif
}     
#endif

#pragma weak rpn_comm_rebind_process__=rpn_comm_rebind_process
#pragma weak rpn_comm_rebind_process_=rpn_comm_rebind_process
void rpn_comm_rebind_process__(void);
void rpn_comm_rebind_process_(void);
void rpn_comm_rebind_process(void)  /* rebind thread(s) of MPI rank */
{
#if defined(linux)
  char *omp;
  int nthreads;
  int ncores, numa, size, rank, i, color, threadpercore, cpucores;
  MPI_Comm comm;
  int t[512];
  uint32_t regs[4];

  if(lo_core == -1) {   /* first call must be by thread 0 in single thread mode */
    nthreads = 1;
    omp = getenv("OMP_NUM_THREADS");
    if(omp != NULL) nthreads=atoi(omp);

    X86_cpuid( 0x0B, 0, regs );
    threadpercore = regs[1] & 0xFFFF;        /* virtual threads per core */
    if(threadpercore < 1) threadpercore = 1;
    X86_cpuid( 0x0B, 1, regs );
    cpucores = regs[1] & 0xFFFF;
    cpucores /= threadpercore;               /* "real" cores */
    if(cpucores < 1) cpucores = 1;
    ncores=sysconf(_SC_NPROCESSORS_CONF);    /* number of cores in node */
    ncores /= threadpercore;                 /* hyperthread factor */
    numa = cpucores;                         /* number of cores per socket */

    color=gethostid();                       /* host id */
    color = (color > 0) ? color : (-color);  /* make color positive */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_split(MPI_COMM_WORLD,color,rank,&comm);
    MPI_Comm_size(comm, &size);              /* size ranks on SMP node */
    MPI_Comm_rank(comm, &rank);              /* my rank on SMP node */

    t[0] = 0;
    for(i=1;i<=size;i++) {
      t[i] = t[i-1] + t[i];
      if(t[i] < numa && t[i] + t[i+1] > numa) t[i]=numa;   /* do no straddle sockets */
      if(t[i] > ncores) t[i] = ncores;                     /* min(t[i],ncores) */
    }
    CPU_ZERO(&set);
    for(i=t[rank] ; i < t[rank+1] ; i++) { CPU_SET(i,&set) ;}
    lo_core = t[rank];
    hi_core = t[rank+1]-1;
  }
  sched_setaffinity(0,sizeof(set),&set);   /* set thread affinity range */
  
#endif
return;
}

#pragma weak rpn_comm_unbind_process__=rpn_comm_unbind_process
#pragma weak rpn_comm_unbind_process_=rpn_comm_unbind_process
void rpn_comm_unbind_process__(void);
void rpn_comm_unbind_process_(void);
void rpn_comm_unbind_process(void)    /* unbind all processes/threads */
{
#if defined(linux)

cpu_set_t set;
int i;
int will_print=1;
char *ompi;
int ncores=sysconf(_SC_NPROCESSORS_CONF);
int nthreads=1;
int nbound=0;
char *omp=getenv("OMP_NUM_THREADS");
int rank;

if(omp != NULL) nthreads=atoi(omp);  /* OMP_NUM_THREADS value */

//                 ompi = getenv("OMPI_COMM_WORLD_RANK");  /* openmpi */
// if(ompi == NULL) ompi = getenv("PMI_RANK");              /* mpich family */
// if(ompi == NULL) ompi = getenv("ALPS_APP_PE");           /* ALPS/mpich family */
// 
// if(ompi != NULL) if(0 != atoi(ompi)) will_print=0;  /* not MPI process 0, no message */

MPI_Comm_rank(MPI_COMM_WORLD, &rank);   /* get MPI rank */
will_print = (rank == 0);

CPU_ZERO(&set);
sched_getaffinity(0,sizeof(set),&set);
i=ncores;
while(--i >=0) { if (CPU_ISSET(i,&set)){ nbound++; } }  /* how many cores are we allowed to run on ? */
  
if(getenv("FULL_UNBIND") != NULL) nbound = 0;  /* FULL_UNBIND variable defined, unbind no matter what */
  
if(nthreads > nbound) {  /* need more threads than cores we can run on , unbind everything */
  if(will_print) printf("INFO: full unbinding will be done, cores=%d, threads needed=%d, usable cores=%d\n",ncores,nthreads,nbound);
  CPU_ZERO(&set);
  i=ncores;
  while(--i >=0) { CPU_SET(i,&set);}
  sched_setaffinity(0,sizeof(set),&set);  /* set affinity to all cores */
}else{
  if(will_print) printf("INFO: no unbinding will be done\n");  /* enough resources available and no forced unbind */
}

#endif

return;
}
