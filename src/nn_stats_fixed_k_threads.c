/*
 *  nn_stats_fixed_k_threads.c
 *
 *  Created by Nicolas Garnier on 2024/10/07.
 *  Copyright 2024 ENS-Lyon - CNRS. All rights reserved.
 *
 *
 *  2024-10-29 - to-do: check if there is a speed improvement when using single_k for multi_k with nk=1
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
    {   int core;       // keep track of the current thread number/id
        int i_start;    // begining of subset of points to work on
        int i_end;      // end      of subset of points
        int nx;         // dimensionality of location data
        int nA;         // dimensionality of observables
        int order_max;  // maximum order of moments to be computed
        int do_center;  // for central moments or moments from the origin
        int k;          // nb of neighbors to search for
    };

struct thread_output
    {   int n_eff;
    };


// dangerous: global variables:
array pos_out;      // for the locations (examination locations)
array obs_in;       // for the observables at the initial locations
k_vector k_in_vec;  // 2024-12-16 to be simplified with declaratino above
array rad_out;      // for the radius of the k-nn (output)
array obs_moments;  // 2025-01-13: for local moments of order 1, 2, ..., order_max, on the outpout locations (output)

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
        order_max= args->order_max,
        do_center= args->do_center,
        n_eff    = i_end-i_start; // how many points in this thread

    double queryPt[nx]; // to be optimized

    for (i=i_start; i<i_end; i++)
    {   for (d=0; d<nx; d++) queryPt[d] = pos_out.A[i + d*pos_out.Npts];
        ANN_compute_stats_single_k(queryPt, obs_in.A, k, rad_out.A + i, obs_moments.A + i, order_max, obs_moments.Npts, nA, do_center, core);
    }
    
    out->n_eff    = n_eff;
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
/* y        contains all the "position"/"location where statistics will be computed     */
/* npts_out is the number of points in y (nb of positions/locations)                    */
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
int compute_stats_fixed_k_threads(double *x, double *A, int npts_in, int nx, int nA, double *y, int npts_out, int k, double *A_moments, int order_max, int do_center, double *dists)
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
    obs_moments.Npts =npts_out; obs_moments.dim =nA;    obs_moments.A =A_moments;
    
    for (core=0; core<nb_cores; core++)
    {   my_arguments[core].core     = core;
        my_arguments[core].i_start  = core*npts_eff_min;
        my_arguments[core].i_end    = (core+1)*npts_eff_min;
        if (core==(nb_cores-1)) my_arguments[core].i_end = npts_out;    // last thread will work longer!
        my_arguments[core].nx       = nx;
        my_arguments[core].nA       = nA;
        my_arguments[core].order_max= order_max;
        my_arguments[core].do_center= do_center;
        my_arguments[core].k        = k;
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
/* function to be used by "compute_stats_multi_k_threads"                               */
/*                                                                                      */
/* 2021-11-26  first multi-threads version                                              */
/****************************************************************************************/
void *threaded_stats_multi_k_func(void *ptr)
{   struct thread_args  *args = (struct thread_args *)ptr; // cast arguments to the usable struct
    struct thread_output *out = calloc(sizeof(struct thread_output),1); // allocate heap memory for this thread's results
    register int i,d;
    int core     = args->core,
        i_start  = args->i_start,
        i_end    = args->i_end,
        nx       = args->nx,
        nA       = args->nA,
        order_max= args->order_max,
        do_center= args->do_center,
        n_eff    = i_end-i_start; // how many points in this thread

    double queryPt[nx]; // to be optimized

    for (i=i_start; i<i_end; i++)
    {   for (d=0; d<nx; d++) queryPt[d] = pos_out.A[i + d*pos_out.Npts];
        ANN_compute_stats_multi_k(queryPt, obs_in.A, k_in_vec, rad_out.A + i, obs_moments.A + i, order_max, obs_moments.Npts, nA, do_center, core);
//        printf("%1.0f %1.0f -> %1.2f\n", pos_out.A[i], pos_out.A[i+pos_out.Npts], rad);
    }
    
    out->n_eff    = n_eff;
    pthread_exit(out);
} /* en of "threaded_stats_multi_k_func" function */



/****************************************************************************************/
/* computes some local statistics, using nearest neighbors                              */
/* same as above function, but for multiple values of k                                 */
/*                                                                                      */
/* x        contains all the "position"/"location data, which is of size npts*nx        */
/* A        contains all the observable data, which is of size npts*nA                  */
/* npts     is the number of points (nb of positions/locations)                         */
/* nx       is the dimensionality of the positions/locations                            */
/* nA       is the dimensionality of the observables (usually, nA=1, but may be larger) */
/* y        contains all the "position"/"location where statistics will be computed     */
/* npts_out is the number of points in y (nb of positions/locations)                    */
/* k        contains all the numbers of neighbors to consider                           */
/* nk       is the number of k to consider                                              */
/*                                                                                      */
/* data is ordered like this :                                                          */
/* x1(t=0)...x1(t=npts-1) x2(t=0) ... x2(t=npts-1) ... xnx(t=0) ... xnx(t=npts-1)       */
/* A1(t=0)...A1(t=npts-1) A2(t=0) ... A2(t=npts-1) ... AnA(t=0) ... AnA(t=npts-1)       */
/*                                                                                      */
/* 2024-10-07  fork from "compute_entropy_ann_threads", no output yet                   */
/* 2024-10-28  fork from "compute_stats_fixed_k_threads", tested OK                     */
/* 2025-01-13 - replaced "mean" and "var" by "moments" and "order_max"                  */
/*              var and mean (and larger order moments) will be saved in moments        */
/*              order_max is the largest order of the moments to be computed            */
/****************************************************************************************/
int compute_stats_multi_k_threads(double *x, double *A, int npts_in, int nx, int nA, double *y, int npts_out, int *k, int nk, double *A_moments, int order_max, int do_center, double *dists)
{	register int core, nb_cores=get_cores_number(GET_CORES_SELECTED), npts_eff_min;
	int n_total=0; // just for sanity check
    int k_max=k[nk-1];
	int ret;
    
    if (nb_cores<1) set_cores_number(nb_cores); // auto-detect, if asked for
    nb_cores=get_cores_number(GET_CORES_SELECTED);        
    // 2022-12-13: all other threads number manipulation should be done outside of this engine function!
    npts_eff_min   = (npts_out - (npts_out%nb_cores))/nb_cores;  // nb pts mini dans chaque thread

	init_ANN(npts_in, nx, k_max, nb_cores);
    create_kd_tree(x, npts_in, nx);
    
    pthread_t    thread[nb_cores];
    struct thread_args   my_arguments[nb_cores];
    struct thread_output *my_outputs[nb_cores];
    
    // dangerous: we use global variables:
    pos_out.Npts =npts_out; pos_out.dim =nx;    pos_out.A =y;
    obs_in.Npts  =npts_in;  obs_in.dim  =nA;    obs_in.A  =A;
    k_in_vec.N   =nk;       k_in_vec.ind_min=0;     k_in_vec.ind_max=nk;        k_in_vec.A = k;
    rad_out.Npts =npts_out; rad_out.dim =nk;    rad_out.A =dists;       // note 2024-10-29: all requested k are returned
    obs_moments.Npts =npts_out; obs_moments.dim =nA;    obs_moments.A =A_moments;

    for (core=0; core<nb_cores; core++)
    {   my_arguments[core].core    = core;
        my_arguments[core].i_start = core*npts_eff_min;
        my_arguments[core].i_end   = (core+1)*npts_eff_min;
        if (core==(nb_cores-1)) my_arguments[core].i_end = npts_out;    // last thread will work longer!
        my_arguments[core].nx = nx;
        my_arguments[core].nA = nA;
        my_arguments[core].order_max = order_max;
        my_arguments[core].do_center = do_center;
        ret=pthread_create(&thread[core], NULL, threaded_stats_multi_k_func, (void *)&my_arguments[core]);
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


#include "nn_stats_kernel_fixed_k_threads.c"
