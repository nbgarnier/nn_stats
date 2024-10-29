/*
 *  nn_stats_fixed_k_threads.c
 *
 *  Created by Nicolas Garnier on 2024/10/07.
 *  Copyright 2024 ENS-Lyon - CNRS. All rights reserved.
 *
 *
 */
 
#include <stdlib.h>
#include <stdio.h>      // for printf
#include <string.h>     // for malloc
#include <pthread.h>

#include "ANN_threads.h"
#include "ANN_wrapper.h"
#include "ANN_stats.h"
#include "data.h"       // for data structs

// #include "library_commons.h"   // for is_zero()

// thread arguments and outputs types:
struct thread_args
    {   int core;    // keep track of the current thread number/id
        int i_start; // begining of subset of points to work on
        int i_end;   // end      of subset of points
        int nx;      // dimensionality of location data
        int nA;      // dimensionality of observables
        int k;       // nb of neighbors to search for
    };

struct thread_output
    {   int n_eff;
        int n_errors;
        double mean;
        double var;
    };


// dangerous: global variables:
array pos_out;      // for the locations (examination locations)
array obs_in;       // for the observables at the initial locations
array rad_out;      // for the radius of the k-nn (output)
array obs_mean;     // for the local average, on the outpout locations (output)
array obs_var;      // for the lcoal variance, on the outpout locations (output)

/****************************************************************************************/
/* function to be used by "compute_stats_fixed_k_threads"                               */
/*                                                                                      */
/* 2021-11-26  first multi-threads version                                              */
/****************************************************************************************/
void *threaded_stats_fixed_k_func(void *ptr)
{   struct thread_args  *args = (struct thread_args *)ptr; // cast arguments to the usable struct
    struct thread_output *out = calloc(sizeof(struct thread_output),1); // allocate heap memory for this thread's results
    register int i,d;
    int core     = args->core,
        i_start  = args->i_start,
        i_end    = args->i_end,
        nx       = args->nx,
        nA       = args->nA,
        k        = args->k,
        n_eff    = i_end-i_start; // how many points in this thread
    double rad   = 0;

    double queryPt[nx]; // to be optimized

    for (i=i_start; i<i_end; i++)
    {   for (d=0; d<nx; d++) queryPt[d] = pos_out.A[i + d*pos_out.Npts];
        rad = ANN_compute_stats_single_k(queryPt, obs_in.A, k, rad_out.A + i, obs_mean.A + i, obs_var.A + i, obs_mean.Npts, obs_mean.dim, core);
//        printf("%1.0f %1.0f -> %1.2f\n", pos_out.A[i], pos_out.A[i+pos_out.Npts], rad);
    }
    
    out->n_eff    = n_eff;
//    out->n_errors = l_errors;
        
//    printf("\t\t%d-%d=%d\n", i_start, i_end, i_end-i_start);
//        free(args);
    pthread_exit(out);
}




/****************************************************************************************/
/* computes some local statistics, using nearest neighbors                              */
/*                                                                                      */
/* x        contains all the "position"/"location data, which is of size npts*nx        */
/* A        contains all the observable data, which is of size npts*nA                  */
/* npts     is the number of points (nb of positions/locations)                         */
/* nx       is the dimensionality of the positions/locations                            */
/* nA       is the dimensionality of the observables (usually, nA=1, but may be larger) */
/* k        is the number of neighbors to consider                                      */
/* nb_cores is the nb of threads to use                                                 */
/*          if ==-1, then it will be set to max nb of threads (auto-detect)             */
/*                                                                                      */
/* data is ordered like this :                                                          */
/* x1(t=0)...x1(t=npts-1) x2(t=0) ... x2(t=npts-1) ... xnx(t=0) ... xnx(t=npts-1)       */
/* A1(t=0)...A1(t=npts-1) A2(t=0) ... A2(t=npts-1) ... AnA(t=0) ... AnA(t=npts-1)       */
/*                                                                                      */
/* 2024-10-07  fork from "compute_entropy_ann_threads", no output yet                   */
/****************************************************************************************/
int compute_stats_fixed_k_threads(double *x, double *A, int npts_in, int nx, int nA, double *y, int npts_out, int k, double *A_mean, double *A_std, double *dists)
{	register int core, nb_cores=get_cores_number(GET_CORES_SELECTED), npts_eff_min;
	int n_total=0; // just for sanity check
	int ret;
    
    if (nb_cores<1) set_cores_number(nb_cores); // auto-detect, if asked for
    nb_cores=get_cores_number(GET_CORES_SELECTED);        
    // 2022-12-13: all other threads number manipulation should be done outside of this engine function!
    npts_eff_min   = (npts_out - (npts_out%nb_cores))/nb_cores;  // nb pts mini dans chaque thread

	init_ANN(npts_in, nx, k, nb_cores); 	
    create_kd_tree(x, npts_in, nx);
    
    pthread_t    thread[nb_cores];
    struct thread_args   my_arguments[nb_cores];
    struct thread_output *my_outputs[nb_cores];
    
    // dangerous: we use global variables:
    pos_out.Npts =npts_out; pos_out.dim =nx;    pos_out.A =y;
    obs_in.Npts  =npts_in;  obs_in.dim  =nA;    obs_in.A  =A;
    rad_out.Npts =npts_out; rad_out.dim =k;     rad_out.A =dists;       // note 2024-10-15: only one k (largest) is returned
    obs_mean.Npts=npts_out; obs_mean.dim=nA;    obs_mean.A=A_mean;
    obs_var.Npts =npts_out; obs_var.dim =nA;    obs_var.A =A_std;
    
    for (core=0; core<nb_cores; core++)
    {   my_arguments[core].core    = core;
        my_arguments[core].i_start = core*npts_eff_min;
        my_arguments[core].i_end   = (core+1)*npts_eff_min;
        if (core==(nb_cores-1)) my_arguments[core].i_end = npts_out;    // last thread will work longer!
        my_arguments[core].nx = nx;
        my_arguments[core].nA = nA;
        my_arguments[core].k = k;
        ret=pthread_create(&thread[core], NULL, threaded_stats_fixed_k_func, (void *)&my_arguments[core]);
        if (ret!=0)
        {   printf("[compute_stats_fixed_k_threads] TROUBLE! couldn't create thread!\n");
            return(-1); 
        }
    }
    for (core=0; core<nb_cores; core++)
    {   pthread_join(thread[core], (void**)&my_outputs[core]);

        n_total += my_outputs[core]->n_eff; // just for sanity
        free(my_outputs[core]);
    }
    
    // sanity check:
    if (n_total!=npts_out) printf("[compute_stats_fixed_k_threads] TROUBLE! npts altered!\n");
       
	free_ANN(nb_cores);
    return(0);
} /* end of function "compute_stats_fixed_k_threads" *************************************/



/****************************************************************************************/
/* computes some local statistics, using nearest neighbors                              */
/* same as above function, but for multiple values of k                                 */
/*                                                                                      */
/* x        contains all the "position"/"location data, which is of size npts*nx        */
/* A        contains all the observable data, which is of size npts*nA                  */
/* npts     is the number of points (nb of positions/locations)                         */
/* nx       is the dimensionality of the positions/locations                            */
/* nA       is the dimensionality of the observables (usually, nA=1, but may be larger) */
/* k        contains all the numbers of neighbors to consider                           */
/* nk       is the number of k to consider                                              */
/* nb_cores is the nb of threads to use                                                 */
/*          if ==-1, then it will be set to max nb of threads (auto-detect)             */
/*                                                                                      */
/* data is ordered like this :                                                          */
/* x1(t=0)...x1(t=npts-1) x2(t=0) ... x2(t=npts-1) ... xnx(t=0) ... xnx(t=npts-1)       */
/* A1(t=0)...A1(t=npts-1) A2(t=0) ... A2(t=npts-1) ... AnA(t=0) ... AnA(t=npts-1)       */
/*                                                                                      */
/* 2024-10-07  fork from "compute_entropy_ann_threads", no output yet                   */
/****************************************************************************************/
int compute_stats_multi_k_threads(double *x, double *A, int npts_in, int nx, int nA, double *y, int npts_out, int *k, int nk, double *A_mean, double *A_std, double *dists)
{	register int core, nb_cores=get_cores_number(GET_CORES_SELECTED), npts_eff_min;
	int n_total=0; // just for sanity check
	int ret;
    
    if (nb_cores<1) set_cores_number(nb_cores); // auto-detect, if asked for
    nb_cores=get_cores_number(GET_CORES_SELECTED);        
    // 2022-12-13: all other threads number manipulation should be done outside of this engine function!
    npts_eff_min   = (npts_out - (npts_out%nb_cores))/nb_cores;  // nb pts mini dans chaque thread

	init_ANN(npts_in, nx, k, nb_cores); 	
    create_kd_tree(x, npts_in, nx);
    
    pthread_t    thread[nb_cores];
    struct thread_args   my_arguments[nb_cores];
    struct thread_output *my_outputs[nb_cores];
    
    // dangerous: we use global variables:
    pos_out.Npts =npts_out; pos_out.dim =nx;    pos_out.A =y;
    obs_in.Npts  =npts_in;  obs_in.dim  =nA;    obs_in.A  =A;
    rad_out.Npts =npts_out; rad_out.dim =k;     rad_out.A =dists;       // note 2024-10-15: only one k (largest) is returned
    obs_mean.Npts=npts_out; obs_mean.dim=nA;    obs_mean.A=A_mean;
    obs_var.Npts =npts_out; obs_var.dim =nA;    obs_var.A =A_std;
    
    for (core=0; core<nb_cores; core++)
    {   my_arguments[core].core    = core;
        my_arguments[core].i_start = core*npts_eff_min;
        my_arguments[core].i_end   = (core+1)*npts_eff_min;
        if (core==(nb_cores-1)) my_arguments[core].i_end = npts_out;    // last thread will work longer!
        my_arguments[core].nx = nx;
        my_arguments[core].nA = nA;
        my_arguments[core].k = k;
        ret=pthread_create(&thread[core], NULL, threaded_stats_fixed_k_func, (void *)&my_arguments[core]);
        if (ret!=0)
        {   printf("[compute_stats_multi_k_threads] TROUBLE! couldn't create thread!\n");
            return(-1); 
        }
    }
    for (core=0; core<nb_cores; core++)
    {   pthread_join(thread[core], (void**)&my_outputs[core]);

        n_total += my_outputs[core]->n_eff; // just for sanity
        free(my_outputs[core]);
    }
    
    // sanity check:
    if (n_total!=npts_out) printf("[compute_stats_multi_k_threads] TROUBLE! npts altered!\n");
       
	free_ANN(nb_cores);
    return(0);
} /* end of function "compute_stats_multi_k_threads" *************************************/
