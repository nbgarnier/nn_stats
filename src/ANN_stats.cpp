//----------------------------------------------------------------------
//	ANN_stats.cpp
//
// do to a bug in some compilers (especially on macos), symbols from 
// multiple objects cannot be imported properly in python
// => as a quick fix, I import this file into ANN_wrapper.cpp
//
// 2024-10-07: new functions for local averaging ("nn_stats")
//----------------------------------------------------------------------
#define noDEBUG
#define DEBUG_N 37

#include "ANN/ANN.h"
#include "ANN_wrapper.h"        // definitions of functions only
#include "ANN_stats.h"          // definitions of functions only
#include <stdio.h>              // for printf, to be removed

// #define UNUSED(expr) do { (void)(expr); } while (0)
// https://stackoverflow.com/questions/1486904/how-do-i-best-silence-a-warning-about-unused-variables


// #include "ann_1.1.2/src/pr_queue_k.h"	// 2021-12-01, k-element priority queue, for definition of ANNmin_k

// global variables for (internal, but lower level so "private") operations on the main tree:
// note that for these variables, moved from other .cpp files, the names have been kept for
// consistency accross the library
                                    
// global variables for the main tree, internal to the C++ code, but exposed only in ANN_Wrapper.cpp:
extern double	        ANN_eps;    // error bound, exact search if 0
// extern ANNpointArray	dataPts;        // data points, type (*(*double))
extern ANNidxArray	   *nnIdx;      // k-nn indices    // 2021, adapted for pthread 
extern ANNdistArray   *dists;       // k-nn distances  // 2021, adapted for pthread
extern ANNkd_tree*	    kdTree;     // search structure



/***************************************************************************************/
/* below : piece of code to average observables over nearest neighbors of a point      */
/*         using a previously computed kd-tree (with ANN library)                      */
/* faster and memory efficient coding                                                  */
/*                                                                                     */
/* input parameters:                                                                   */
/* x    : is a d-dimensional query point (the one to search neighbors of)              */
/* A    : is a vector of observable of size Npts,                                      */
/*          where Npts is the same as the location data the tree is built upon         */
/* k    : is the rank of the neighbor to search for                                    */
/* npts_out: nb of points in the output                                                */
/* nA      : nb of observables                                                         */
/* output parameters:                                                                  */
/* R    : distance of the k-nn from x                                                  */
/* mean : expected values of observables                                               */
/* var  : variance of observables                                                      */
/*                                                                                     */
/* 2024-10-07 - first version                                                          */
/* 2024-10-14 - draft for returning distance : to do: check max index                  */
/* 2024-10-15 - factorized loop on observables                                         */
/***************************************************************************************/
double ANN_compute_stats_single_k(double *x, double *A, int k, double *R, double *mean, double *var, int npts_out, int nA, int core)
{   int i, d, ind, N=k-1+ANN_ALLOW_SELF_MATCH;
    int npts=kdTree->nPoints();
    double tmp, m=0., v=0.;

    kdTree->annkSearch(x,                       // query point
                       k+ANN_ALLOW_SELF_MATCH,  // number of near neighbors (including or excluding central point)
                       nnIdx[core],             // nearest neighbors (returned)
                       dists[core],             // distance (returned)
                       ANN_eps);

    for (d=0; d<nA; d++)
    {   m=0.; v=0.;
        for (i=0; i<N; i++)
        {   ind = nnIdx[core][i];
            tmp = (A+npts*d)[ind]; 
            m  += tmp;  v += tmp*tmp;
        }
        v -= m*m/N;
        m /= N;
        v /= N-1; // unbiased estimator

        mean[npts_out*d] = m;
        var [npts_out*d] = v;

    }
    R[0] = (double)dists[core][N-1];
//    printf("k=%d, N=%d, allow=%d   ", k, N, ANN_ALLOW_SELF_MATCH);
//    return(0);
    return((double)dists[core][N-1]);
} /* end of function "ANN_compute_stats_single_k" ***********************************************/




/****************************************************************************************/
/* same as above, but for a set of prescribed values of numbers k                       */
/*                                                                                      */
/* 2024-10-15 - factorized loop on observables                                          */
/****************************************************************************************/
double ANN_compute_stats_multi_k(double *x, double *A, int *k, int Nk, double *R, double *mean, double *var, int npts_out, int nA, int core)
{   int i, d, ind, N;
    int npts=kdTree->nPoints();
    double tmp, m=0., v=0.;
    int ind_k, k_max=k[Nk-1];

    kdTree->annkSearch(x,                       // query point
                       k_max+ANN_ALLOW_SELF_MATCH,  // number of near neighbors (including or excluding central point)
                       nnIdx[core],             // nearest neighbors (returned)
                       dists[core],             // distance (returned)
                       ANN_eps);

    for (ind_k=0; ind_k<Nk; ind_k++)
    {   N=k[ind_k]-1+ANN_ALLOW_SELF_MATCH;
        for (d=0; d<nA; d++)
        {   m=0.; v=0.;
            for (i=0; i<N; i++)
            {   ind = nnIdx[core][i];
                tmp = (A+npts*d)[ind]; 
                m  += tmp;  v += tmp*tmp;
            }
            v -= m*m/N;
            m /= N;
            v /= N-1; // unbiased estimator

            mean[npts_out*d] = m;
            var [npts_out*d] = v;

        }
        R[0] = (double)dists[core][N-1];
//    printf("k=%d, N=%d, allow=%d   ", k, N, ANN_ALLOW_SELF_MATCH);
//    return(0);
    }

    return((double)dists[core][N-1]);
} /* end of function "ANN_compute_stats_multi_k" ***********************************************/
