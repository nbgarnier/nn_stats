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
// #include <string.h>     // for malloc
#include <pthread.h>

#include "ANN_threads.h"
#include "ANN_wrapper.h"
#include "ANN_stats.h"
#include "ANN_stats_kernel.h"
#include "library_commons.h"    // for tree_k_max
#include "kernels.h"            // for kernels definitions
#include "data.h"               // for data structs

#define noDEBUG // set to "DEBUG". can also set to "DEBUG_N"
#define i_look 388 // for debug with "DEBUG_N"

void print_vec_int(int *v, int N)
{   int i;
    printf("(");
    for (i=0; i<N; i++) printf("%d ", v[i]);
    printf(")");
}
void print_vec(double *v, int N)
{   int i;
    printf("(");
    for (i=0; i<N; i++) printf("%1.2f ", v[i]);
    printf(")");
}

// thread arguments and outputs types:
struct thread_args
    {   int core;       // keep track of the current thread number/id
        int i_start;    // begining of subset of points to work on
        int i_end;      // end      of subset of points
        int nx;         // dimensionality of location data
        int nA;         // dimensionality of observables
        int order_max;  // maximum order of moments to be computed
        int do_center;  // for central moments or moments from the origin
        double R;       // radius of ball
    };

struct thread_output
    {   int n_eff;
    };


// dangerous: global variables:
// first ones are already defined in "nn_stats_fixed_k_threads.c":
extern array   pos_out;     // for the locations (examination locations)
extern array   obs_in;      // for the observables at the initial locations
extern array obs_moments;   // 2025-01-13: for local moments of order 1, 2, ..., order_max, on the outpout locations (output)
// last ones are specific to this file:
array   R_in;               // for the radii (input)
arr_int nnn_out;            // for the nb of nn (output)


/****************************************************************************************/
/* function to be used by "compute_stats_fixed_R_threads"                               */
/*                                                                                      */
/* 2021-11-26  first multi-threads version                                              */
/****************************************************************************************/
void *threaded_stats_fixed_R_func(void *ptr)
{   struct thread_args  *args = (struct thread_args *)ptr; // cast arguments to the usable struct
    struct thread_output *out = calloc(sizeof(struct thread_output),1); // allocate heap memory for this thread's results
    register int i,d,j_moment,k;
    int core     = args->core,
        i_start  = args->i_start,
        i_end    = args->i_end,
        nx       = args->nx,
        nA       = args->nA,
        order_max= args->order_max,
        do_center= args->do_center;
    double R     = args->R;
    int    n_eff = i_end-i_start;   // how many points in this thread
    double ratou[2];
    double queryPt[nx];             // to be optimized ?

    for (i=i_start; i<i_end; i++)
    {   for (d=0; d<nx; d++) queryPt[d] = pos_out.A[i + d*pos_out.Npts];
        k = ANN_count_nearest_neighbors(queryPt, R, core);
        nnn_out.A[i] = k;

        if (k>=tree_k_max)
        {   for (d=0; d<nA; d++)
            for (j_moment=0; j_moment<order_max; j_moment++)
                (obs_moments.A+i)[obs_moments.Npts*(nA*j_moment + d)] = my_NAN;
                
        }
        else 
        {   if (current_kernel_type>0)
                ANN_compute_stats_kernel_single_k(queryPt, obs_in.A, k, ratou, obs_moments.A + i, order_max, obs_moments.Npts, nA, do_center, core);
            else
                ANN_compute_stats_single_k       (queryPt, obs_in.A, k, ratou, obs_moments.A + i, order_max, obs_moments.Npts, nA, do_center, core);
        }
    }
    
    out->n_eff    = n_eff;
    pthread_exit(out);
}




/****************************************************************************************/
/* computes some local statistics, using nearest neighbors                              */
/*                                                                                      */
/* x        contains all the "position"/"location data, which is of size npts*nx        */
/* A        contains all the observable data, which is of size npts*nA                  */
/* npts_in  is the number of points (nb of positions/locations)                         */
/* nx       is the dimensionality of the positions/locations                            */
/* nA       is the dimensionality of the observables (usually, nA=1, but may be larger) */
/* y        contains all the "position"/"location where statistics will be computed     */
/* npts_out is the number of points in y (nb of positions/locations)                    */
/* R        the radius to consider                                                      */
/*                                                                                      */
/* data is ordered like this :                                                          */
/* x1(t=0)...x1(t=npts-1) x2(t=0) ... x2(t=npts-1) ... xnx(t=0) ... xnx(t=npts-1)       */
/* A1(t=0)...A1(t=npts-1) A2(t=0) ... A2(t=npts-1) ... AnA(t=0) ... AnA(t=npts-1)       */
/*                                                                                      */
/* 2012-02-27  fork from "compute_entropy_ann"                                          */
/* 2021-11-26  first multi-threads version                                              */
/****************************************************************************************/
int compute_stats_fixed_R_threads(double *x, double *A, int npts_in, int nx, int nA, double *y, int npts_out, double R, double *A_moments, int order_max, int do_center, int *k)
{	register int core, nb_cores=get_cores_number(GET_CORES_SELECTED);
    int npts_eff_min;
	int n_total=0; // just for sanity check
	int ret;
    
    if (nb_cores<1) set_cores_number(nb_cores); // auto-detect, if asked for
    nb_cores=get_cores_number(GET_CORES_SELECTED);          
    // 2022-12-13: all other threads number manipulation should be done outside of this engine function!
    npts_eff_min   = (npts_out - (npts_out%nb_cores))/nb_cores;  // nb pts mini dans chaque thread
    
    if ((tree_k_max<1) || (tree_k_max>=npts_in)) tree_k_max = npts_in/10; // 2024-10-08, quick security, not optimized (memory)
	init_ANN(npts_in, nx, tree_k_max, nb_cores); 	// pb with unknow k...
    create_kd_tree(x, npts_in, nx);

    pthread_t    thread[nb_cores];
    struct thread_args   my_arguments[nb_cores];
    struct thread_output *my_outputs[nb_cores];
    
    // dangerous: we use global variables:
    pos_out.Npts =npts_out; pos_out.dim =nx;    pos_out.A =y;
    obs_in.Npts  =npts_in;  obs_in.dim  =nA;    obs_in.A  =A;
    nnn_out.Npts =npts_out; nnn_out.dim =1;     nnn_out.A =k;       // note 2024-10-15: only one nb (largest R) is returned
    obs_moments.Npts =npts_out; obs_moments.dim =nA;    obs_moments.A =A_moments;
    
    
    for (core=0; core<nb_cores; core++)
    {   my_arguments[core].core     = core;
        my_arguments[core].i_start  = core*npts_eff_min;
        my_arguments[core].i_end    = (core+1)*npts_eff_min;
        if (core==(nb_cores-1)) my_arguments[core].i_end = npts_out;  // last thread will work longer!
        my_arguments[core].nx       = nx;
        my_arguments[core].nA       = nA;
        my_arguments[core].order_max= order_max;
        my_arguments[core].do_center= do_center;
        my_arguments[core].R        = R;
        ret=pthread_create(&thread[core], NULL, threaded_stats_fixed_R_func, (void *)&my_arguments[core]);
        if (ret!=0)
        {   printf("[compute_stats_fixed_R_threads] TROUBLE! couldn't create thread!\n");
            return(-1); 
        }
    }
    for (core=0; core<nb_cores; core++)
    {   pthread_join(thread[core], (void**)&my_outputs[core]);

        n_total += my_outputs[core]->n_eff; // just for sanity
        free(my_outputs[core]);
    }
    
    // sanity check:
    if (n_total!=npts_out) printf("[compute_stats_fixed_R_threads] TROUBLE! npts altered!\n");
       
	free_ANN(nb_cores);
    return(0);
} /* end of function "compute_stats_fixed_R_threads" *************************************/



/****************************************************************************************/
/* function to be used by "compute_stats_multi_R_threads"                               */
/*                                                                                      */
/* 2021-11-26 - first multi-threads version                                             */
/* 2025-01-20 - now with specific kernel, if resuired (globaly)                         */
/****************************************************************************************/
void *threaded_stats_multi_R_func(void *ptr)
{   struct thread_args  *args = (struct thread_args *)ptr; // cast arguments to the usable struct
    struct thread_output *out = calloc(sizeof(struct thread_output),1); // allocate heap memory for this thread's results
    register int i,j,j_moment,d;
    int ind_k_min, ind_k_max;
    int core     = args->core,
        i_start  = args->i_start,
        i_end    = args->i_end,
        nx       = args->nx,
        nA       = args->nA,
        order_max= args->order_max,
        do_center= args->do_center;
    int n_eff    = i_end-i_start;   // how many points in this thread
    k_vector    k_vec;

    double queryPt[nx];             // to be optimized
    k_vec.A = calloc(R_in.dim, sizeof(int));
    k_vec.N = R_in.dim;

    for (i=i_start; i<i_end; i++)
    {   for (d=0; d<nx; d++) queryPt[d] = pos_out.A[i + d*pos_out.Npts];

        j=0;
        do
        {   k_vec.A[j] = ANN_count_nearest_neighbors(queryPt, R_in.A[j], core);
            nnn_out.A[i + obs_moments.Npts*j]= k_vec.A[j];
            j++;
        }
        while ( (j<R_in.dim) && (k_vec.A[j-1]<tree_k_max) );

        // search for the range of valid k values and keep their index:
        ind_k_min=0;
        while ((k_vec.A[ind_k_min]<1) && (ind_k_min<R_in.dim)) ind_k_min++;
        ind_k_max=ind_k_min;
        while ((k_vec.A[ind_k_max]<tree_k_max) && (ind_k_max<R_in.dim)) ind_k_max++;
        
#ifdef DEBUG        
        printf("point %d ", i); print_vec(queryPt, nx); printf(" -> j=%d (k=%d)\t k = ", j, k_vec.A[j-1]);
        print_vec_int(k_vec.A, k_vec.N);
        printf(" k in [ %d - %d [\t", ind_k_min, ind_k_max);
#endif

#ifdef DEBUG_N
        if ((i%i_look)==0)
        {   printf("point %d ( ", i); 
            for (d=0; d<nx; d++) printf("%1.2f ,", queryPt[d]);
            printf(")\t");
            for (j=0; j<R_in.dim; j++) printf("  R=%1.2f -> k=%d", R_in.A[j], k_vec.A[j]);
            printf("\n");
            printf("\t ind_k_min = %d => k_min = %d, \tind_k_max = %d => k_max = %d\n", ind_k_min, k_vec.A[ind_k_min], ind_k_max, k_vec.A[ind_k_max-1]);
        }
#endif

        // work in the range of valid k:
        if ( (ind_k_min<R_in.dim) && (ind_k_min<ind_k_max) )
        {   k_vec.ind_min = ind_k_min;
            k_vec.ind_max = ind_k_max;
            if (current_kernel_type>0)      // then we use a special kernel
                ANN_compute_stats_kernel_multi_k(queryPt, obs_in.A, k_vec, NULL, obs_moments.A+i, order_max, obs_moments.Npts, nA, do_center, core);
            else
                ANN_compute_stats_multi_k       (queryPt, obs_in.A, k_vec, NULL, obs_moments.A+i, order_max, obs_moments.Npts, nA, do_center, core);
                
        }

        int my_min = (ind_k_min<R_in.dim) ? ind_k_min : R_in.dim;
        for (j=0; j<my_min; j++)                          
        {   for (d=0; d<nA; d++)
            for (j_moment=0; j_moment<order_max; j_moment++)
                (obs_moments.A+i)[obs_moments.Npts* (nA*(R_in.dim*j_moment + j) + d)] = my_NAN;
        }
        for (j=ind_k_max; j<R_in.dim; j++)                     
        {   for (d=0; d<nA; d++)
            for (j_moment=0; j_moment<order_max; j_moment++)
                (obs_moments.A+i)[obs_moments.Npts* (nA*(R_in.dim*j_moment + j) + d)] = my_NAN;
        }
    }
    
    out->n_eff    = n_eff;
    free(k_vec.A);
    pthread_exit(out);
}


/****************************************************************************************/
/* computes some local statistics, using nearest neighbors                              */
/* same as above function, but for multiple values of R                                 */
/*                                                                                      */
/* x        contains all the "position"/"location data, which is of size npts*nx        */
/* A        contains all the observable data, which is of size npts*nA                  */
/* npts_in  is the number of points (nb of positions/locations)                         */
/* nx       is the dimensionality of the positions/locations                            */
/* nA       is the dimensionality of the observables (usually, nA=1, but may be larger) */
/* y        contains all the "position"/"location where statistics will be computed     */
/* npts_out is the number of points in y (nb of positions/locations)                    */
/* R        contains all the radii to consider                                          */
/* nR       is the number of radii to consider                                          */
/*                                                                                      */
/* data is ordered like this :                                                          */
/* x1(t=0)...x1(t=nx-1) x2(t=0) ... x2(t=nx-1) ... xn(t=0) ... xn(t=nx-1)               */
/*                                                                                      */
/* 2012-02-27  fork from "compute_entropy_ann"                                          */
/* 2021-11-26  first multi-threads version                                              */
/****************************************************************************************/
int compute_stats_multi_R_threads(double *x, double *A, int npts_in, int nx, int nA, double *y, int npts_out, double *R, int nR, double *A_moments, int order_max, int do_center, int *k)
{	register int core, nb_cores=get_cores_number(GET_CORES_SELECTED);
    int npts_eff_min;
	int n_total=0; // just for sanity check
	int ret;
    
    if (nb_cores<1) set_cores_number(nb_cores); // auto-detect, if asked for
    nb_cores=get_cores_number(GET_CORES_SELECTED);            
    // 2022-12-13: all other threads number manipulation should be done outside of this engine function!
    npts_eff_min   = (npts_out - (npts_out%nb_cores))/nb_cores;  // nb pts mini dans chaque thread
    
    if ((tree_k_max<1) || (tree_k_max>=npts_in)) tree_k_max = npts_in/10; // 2024-10-08, quick security, not optimized (memory)
	init_ANN(npts_in, nx, tree_k_max, nb_cores); 	// pb with unknow k...
    create_kd_tree(x, npts_in, nx);

    pthread_t    thread[nb_cores];
    struct thread_args   my_arguments[nb_cores];
    struct thread_output *my_outputs[nb_cores];
    
    // dangerous: we use global variables:
    pos_out.Npts =npts_out; pos_out.dim =nx;    pos_out.A =y;
    obs_in.Npts  =npts_in;  obs_in.dim  =nA;    obs_in.A  =A;
    R_in.Npts    =1;        R_in.dim    =nR;    R_in.A    =R;
    nnn_out.Npts =npts_out; nnn_out.dim =nR;    nnn_out.A =k;       // note 2024-12-16: all R results are returned
    obs_moments.Npts =npts_out; obs_moments.dim =nA;    obs_moments.A =A_moments;
    
    
    for (core=0; core<nb_cores; core++)
    {   my_arguments[core].core      = core;
        my_arguments[core].i_start   = core*npts_eff_min;
        my_arguments[core].i_end     = (core+1)*npts_eff_min;
        if (core==(nb_cores-1)) my_arguments[core].i_end = npts_out;  // last thread will work longer!
        my_arguments[core].nx        = nx;
        my_arguments[core].nA        = nA;
        my_arguments[core].order_max = order_max;
        my_arguments[core].do_center = do_center;
        ret=pthread_create(&thread[core], NULL, threaded_stats_multi_R_func, (void *)&my_arguments[core]);
        if (ret!=0)
        {   printf("[compute_stats_multi_R_threads] TROUBLE! couldn't create thread!\n");
            return(-1); 
        }
    }
    for (core=0; core<nb_cores; core++)
    {   pthread_join(thread[core], (void**)&my_outputs[core]);
        n_total += my_outputs[core]->n_eff; // just for sanity
        free(my_outputs[core]);
    }
    
    // sanity check:
    if (n_total!=npts_out) printf("[compute_stats_multi_R_threads] TROUBLE! npts altered!\n");
       
	free_ANN(nb_cores);
    return(0);
} /* end of function "compute_stats_multi_R_threads" *************************************/
