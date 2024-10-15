/************************************************************************/
/* ANN_threads.h	        		                                    */
/*									                                    */
/* 2021/11/26                                                           */
/* 2024/10/11 for from "entropy_ann_threads.h"                          */
/************************************************************************/
/* Nicolas Garnier	nicolas.garnier@ens-lyon.fr			                */
/************************************************************************/
#ifndef ANN_THREADS_H
#define ANN_THREADS_H

#include <pthread.h>
#include <unistd.h>     // to probe the nb of processors
#include <stdio.h>

#ifdef _SC_NPROCESSORS_ONLN
// #define NCORES s_SC_NPROCESSORS_ONLN           // bad on Macos
   #define NCORES (int)sysconf(_SC_NPROCESSORS_ONLN)   // good, but not a compile-time constant
#else
   #define NCORES 1
#endif

#define GET_CORES_AVAILABLE  0x0010
#define GET_CORES_SELECTED   0x0001

#define ANN_PTS_PER_THREAD   (117*16)    // see function "adapt_cores_number"

// global variables (defined for real in "entropy_ann_threads.c") 
extern int USE_PTHREAD;
extern int _n_cores; 

int  get_multithreading_state(int verbosity);
void set_multithreading_state(int do_mp);
int  get_cores_number        (int get_what);
void set_cores_number        (int n);
int  best_cores_number       (int npts_eff);

#endif

