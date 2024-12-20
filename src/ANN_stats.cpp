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
#include <iostream>
#include <vector>               // for Debian compile

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
{   int i, d, N=k+ANN_ALLOW_SELF_MATCH-1;
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
        {   tmp = (A+npts*d)[nnIdx[core][i]]; 
            m  += tmp;  v += tmp*tmp;
        }
        v -= m*m/N;
        m /= N;
        v /= N-1; // unbiased estimator

        mean[npts_out*d] = m;
        var [npts_out*d] = v;

    }
    R[0] = (double)dists[core][N-1];
#ifdef DEBUG    
    std::printf("k=%d, N=%d, allow=%d    R[0]=%1.2f  R[1]=%1.2f  R[2]=%1.2f  R[N-1]=%1.2f  R[N]=%1.2f", 
                k, N, ANN_ALLOW_SELF_MATCH, dists[core][0], dists[core][1], dists[core][2], dists[core][N-1], dists[core][N]);       //  fflush(stdout);
//    std::cout << "k=" << k << ", N=" << N << ", allow=" << ANN_ALLOW_SELF_MATCH; fflush(stdout);
#endif 
    return((double)dists[core][N-1]);
} /* end of function "ANN_compute_stats_single_k" ***********************************************/




/****************************************************************************************/
/* same as above, but for a set of prescribed values of numbers k                       */
/*                                                                                      */
/* beware memory usage!!!                                                               */
/*                                                                                      */
/* 2024-10-17 - initial fork                                                            */
/* 2024-10-29 - full rewritting, for optimization                                       */
/* 2024-12-16 - replaced int *k, int Nk by k_vector k_vec                               */
/****************************************************************************************/
double ANN_compute_stats_multi_k(double *x, double *A, k_vector k_vec, double *R, double *mean, double *var, int npts_out, int nA, int core)
{   int i, d, N=1, N_old;
    int npts=kdTree->nPoints();
    double tmp; 
    int ind_k, k_max=k_vec.A[k_vec.ind_max-1];  // !!! k must be sorted, we take the largest

    std::vector<double> m, v;                   // 2024/10/28: to optimize computations, C++ allocation 
    m.resize(nA); v.resize(nA);

    kdTree->annkSearch(x,                       // query point
                       k_max+ANN_ALLOW_SELF_MATCH,  // number of near neighbors (including or excluding central point)
                       nnIdx[core],             // nearest neighbors (returned)
                       dists[core],             // distance (returned)
                       ANN_eps);

    N_old=0;
    for (d=0; d<nA; d++)
    {   m[d]=0.; v[d]=0.;
    }

    for (ind_k=k_vec.ind_min; ind_k<k_vec.ind_max; ind_k++)
    {   N=k_vec.A[ind_k]-1+ANN_ALLOW_SELF_MATCH;
        if (R!=NULL) R[npts_out*ind_k] = (double)dists[core][N-1];

        for (d=0; d<nA; d++)
        {   for (i=N_old; i<N; i++) 
            {   tmp = (A+npts*d)[nnIdx[core][i]]; 
                m[d]  += tmp;  v[d] += tmp*tmp;
            }
            mean[npts_out*(ind_k+d*k_vec.N)] = m[d]/N;
            var [npts_out*(ind_k+d*k_vec.N)] = (v[d]-m[d]*m[d]/N)/(N-1);  // unbiased estimator
        }

        N_old=N;
    }

    return((double)dists[core][N-1]);
} /* end of function "ANN_compute_stats_multi_k" ***********************************************/
