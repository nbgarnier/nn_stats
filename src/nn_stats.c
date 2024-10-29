/*
 *  nn_stats.c
 *  
 *  to compute local statistics with k-nn algorithms
 *  using ANN library : http://www.cs.umd.edu/~mount/ANN/
 *
 *  Created by Nicolas Garnier on 2024/10/07 (fork from the entropy library).
 *  Copyright 2012-2024 ENS-Lyon - CNRS. All rights reserved.
 *
 *  2017-11-29 : renamed "search_ANN" as "search_ANN_internal" for future extensions
 *  2019-01-21 : added test for distance==0 inside ANN_wrapper.c (returned value!=0 if problem)
 *               and rewritten all tests for (nb_errors>=npts)
 *  2019-01-22 : rewritten some ANN functions, and cleaned up a little (tests moved back here)
 *  2021-11-25 : first multithreading tests
 *  2024-10-07 : fork from "entropy_ann.c"
 */

#include <math.h>                   // for fabs and others
#include <string.h>

#include "library_commons.h"        // for definitions of nb_errors, and stds
//#include "library_matlab.h"         // compilation for Matlab
#include "ANN_wrapper.h"            // for ANN library functions (in C++)
// #include "nns_count.h"              // NBG counting functions (2019-01-23)
//#include "math_tools.h"
#include "nn_stats.h"
#include "nn_stats_fixed_k_threads.h"
#include "nn_stats_fixed_R_threads.h"

#define noDEBUG	    // for debug information, replace "noDEBUG" by "DEBUG"
#define noDEBUG_EXPORT
#define LOOK 17 	// for debug also (of which point(s) will we save the data ?)



#ifdef DEBUG
/* 2018-04-11/3: test of memory alignment: */
void debug_trace(char *text, double *x, int npts, int m, int p, int stride, int k)
{    register int i, j;
#define NPTS_T 9
     printf("%s - npts=%d, m=%d, p=%d, stride=%d, k=%d\n", text, npts, m, p, stride, k);
     for (j=0; j<m; j++)
     {    printf("  vector %d: [ ", j);
          for (i=0; i<(NPTS_T>npts?npts:NPTS_T); i++) printf("%f ",x[i+j*npts]);
          printf("] along time\n");
     }
     return;
}
#endif



/****************************************************************************************/
/* computes Shannon entropy, using nearest neighbor statistics (Grassberger 2004)	    */
/*                                                                                      */
/* this version is for m-dimentional systems, with eventually some stride/embedding	    */
/* here, N_eff is imposed, as is tau_Theiler                                            */
/*																			            */
/* x      contains all the data, which is of size nx in time					    	*/
/* nx     is the number of points in time											    */
/* m	  is the (initial) dimensionality of x								            */
/* p	  indicates how many points to take in time (in the past) (embedding)           */
/* tau    indicates the time lag between 2 consecutive points to be considered in time	*/
/* tau_Theiler    : mimimal stride between two sets of points                           */
/* N_eff          : number of points to use (ie, to consider for statistics)            */
/* N_realizations : number of realizations to consider (expected independant)           */
/* k      nb of neighbors to be considered										        */
/*																			            */
/* data is ordered like this :													        */
/* x1(t=0)...x1(t=nx-1) x2(t=0) ... x2(t=nx-1) ... xn(t=0) ... xn(t=nx-1)				*/
/*																			            */
/* this function is a  wrapper to the functions :									    */
/*	- compute_entropy_nd_ann                                                            */
/*																			            */
/* 2011-11-15 : first version													        */
/* 2012-05-05 : added Theiler correction                                                */
/* 2012-06-01 : Theiler correction should also work with p=1 (stride>1) : new test      */
/* 2020-02-26 : added output of standard deviation                                      */
/* 2021-12-08 : using pthread                                                           */
/* 2021-12-17 : new function for embedding                                              */
/* 2022-04-14 : function "compute_entropy_ann_N", forked from "compute_entropy_ann"     */
/* 2023-11-28 : renamed "compute_entropy_ann" (old one is now "_legacy")                */
/****************************************************************************************/
int compute_stats_fixed_k(double *x, double *A, int npts, int nx, int nA, int k, double *A_mean, double *A_var)
{	register int j;
	
#ifdef DEBUG
    debug_trace("[compute_stats_fixed_k] signal x", x, npts, m, p, tau, k);
#endif

    j = compute_stats_fixed_k_threads(x, A, npts, nx, nA, k, A_mean, A_var, get_cores_number(GET_CORES_SELECTED));

	return(j);
} /* end of function "compute_stats_fixed_k" **********************************************/

