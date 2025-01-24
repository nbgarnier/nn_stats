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

#include <stdio.h>              // for printf, to be removed
#include <iostream>
#include <vector>               // for Debian compile
#include "ANN/ANN.h"
#include "ANN_wrapper.h"        // definitions of functions only
#include "ANN_stats.h"          // definitions of functions only


// #define UNUSED(expr) do { (void)(expr); } while (0)
// https://stackoverflow.com/questions/1486904/how-do-i-best-silence-a-warning-about-unused-variables


// #include "ann_1.1.2/src/pr_queue_k.h"	// 2021-12-01, k-element priority queue, for definition of ANNmin_k

// global variables for (internal, but lower level so "private") operations on the main tree:
// note that for these variables, moved from other .cpp files, the names have been kept for
// consistency accross the library
                                    
// global variables for the main tree, internal to the C++ code, but exposed only in ANN_Wrapper.cpp:
extern double	        ANN_eps;    // error bound, exact search if 0
extern ANNidxArray	   *nnIdx;      // k-nn indices    // 2021, adapted for pthread 
extern ANNdistArray    *dists;      // k-nn distances  // 2021, adapted for pthread
extern ANNkd_tree*	    kdTree;     // search structure





int binomial_1[2] = {1, 1};
int binomial_2[3] = {1, 2,  1};
int binomial_3[4] = {1, 3,  3,  1};
int binomial_4[5] = {1, 4,  6,  4,  1};
int binomial_5[6] = {1, 5, 10, 10,  5,  1};
int binomial_6[7] = {1, 6, 15, 20, 15,  6, 1};  
int binomial_7[8] = {1, 7, 21, 35, 35, 21, 7, 1};  

int *get_binomial(int order)
{   if      (order==1) return binomial_1; // a trick
    else if (order==2) return binomial_2;
    else if (order==3) return binomial_3;
    else if (order==4) return binomial_4;
    else if (order==5) return binomial_5;
    else if (order==6) return binomial_6;
    else if (order==7) return binomial_7;
    else return NULL;
}



/***************************************************************************************/
/* below : piece of code to average observables over nearest neighbors of a point      */
/*         using a previously computed kd-tree (with ANN library)                      */
/* faster and memory efficient coding                                                  */
/*                                                                                     */
/* input parameters:                                                                   */
/* x         : is a d-dimensional query point (the one to search neighbors of)         */
/* A         : is a vector of observable of size Npts,                                 */
/*              where Npts is the same as the location data the tree is built upon     */
/* k         : is the rank of the neighbor to search for                               */
/* npts_out  : nb of points in the output                                              */
/* nA        : nb of observables                                                       */
/* order_max : maximal order of moments to be computed                                 */
/* do_center : central moments if ==1, or not-centered moments if ==0                  */
/* core      : which core to run on                                                    */
/* output parameters:                                                                  */
/* R         : distance of the k-nn from x                                             */
/* moments   : expected values of observables                                          */
/* var       : variance of observables                                                 */
/*                                                                                     */
/* 2024-10-07 - first version                                                          */
/* 2024-10-14 - draft for returning distance : to do: check max index                  */
/* 2024-10-15 - factorized loop on observables                                         */
/***************************************************************************************/
double ANN_compute_stats_single_k(double *x, double *A, int k, double *R, double *moments, int order_max, int npts_out, int nA, int do_center, int core)
{   int i, d, N=k+ANN_ALLOW_SELF_MATCH-1, j_moments, l;
    int npts=kdTree->nPoints();
    double tmp, prod, mean;

    std::vector<double> mom;                    // 2025/01/13: to optimize computations, C++ allocation 
    mom.resize(order_max+1);                    // 2025/01/14: we add the moment of order 0 for ease of use

    kdTree->annkSearch(x,                       // query point
                       k+ANN_ALLOW_SELF_MATCH,  // number of near neighbors (including or excluding central point)
                       nnIdx[core],             // nearest neighbors (returned)
                       dists[core],             // distance (returned)
                       ANN_eps);

    R[0] = (double)dists[core][N-1];
    
    if (order_max>0)                            // added 2025-01-16 for robustness
    for (d=0; d<nA; d++)
    {   mom[0] = (double)N;                            // convention for moment of order 0
        for (j_moments=1; j_moments<=order_max; j_moments++)
        {   mom[j_moments]=0.;
        }
        
        for (i=0; i<N; i++)
        {   tmp  = (A+npts*d)[nnIdx[core][i]];

            prod = 1.;
            for (j_moments=1; j_moments<=order_max; j_moments++)
            {   prod *= tmp;
                mom[j_moments] += prod;
            }
        }

        if (do_center==0)       // non-centered moments (moments about the origin):
        {   for (j_moments=0; j_moments<order_max; j_moments++)
            {   moments[npts_out*(nA*j_moments + d)] = mom[j_moments+1]/N;
            }
        }
        else                    // central moments: https://en.wikipedia.org/wiki/Central_moment
        {   if (order_max>0)    
            {   mean = mom[1]/N;                        // mean (expected value) or initial data
                moments[npts_out*(nA*0 + d)] = mean;      // mean of centered data is 0, but we return the real mean (2025-01-24)
            }
            for (j_moments=1; j_moments<order_max; j_moments++)
            {   moments[npts_out*(nA*j_moments + d)] = 0.0;
                for (l=0; l<=j_moments+1; l++)
                   moments[npts_out*(nA*j_moments + d)] += get_binomial(j_moments+1)[l] * mom[l]/N * pow(-mean, j_moments+1-l);
            }
        }

    }

#ifdef DEBUG    
    std::printf("k=%d, N=%d, allow=%d    R[0]=%1.2f  R[1]=%1.2f  R[2]=%1.2f  R[N-1]=%1.2f  R[N]=%1.2f", 
                k, N, ANN_ALLOW_SELF_MATCH, dists[core][0], dists[core][1], dists[core][2], dists[core][N-1], dists[core][N]);       //  fflush(stdout);
//    std::cout << "k=" << k << ", N=" << N << ", allow=" << ANN_ALLOW_SELF_MATCH; fflush(stdout);
#endif 
    return((double)dists[core][N-1]);
} /* end of function "ANN_compute_stats_single_k" ***********************************************/




/****************************************************************************************/
/* same as above, but for a set of prescribed values of numbers k                       */
/* parameter k is replaced by k_vec                                                     */
/*                                                                                      */
/* beware memory usage!!!                                                               */
/*                                                                                      */
/* 2024-10-17 - initial fork                                                            */
/* 2024-10-29 - full rewritting, for optimization                                       */
/* 2024-12-16 - replaced int *k, int Nk by k_vector k_vec                               */
/* 2025-01-13 - replaced "mean" and "var" by "moments" and "order_max"                  */
/*              var and mean (and larger order moments) will be returned in moments     */
/*              order_max is the largest order of the moments to be computed            */
/****************************************************************************************/
double ANN_compute_stats_multi_k(double *x, double *A, k_vector k_vec, double *R, double *moments, int order_max, int npts_out, int nA, int do_center, int core)
{   int i, d, N=1, N_old, j_moments, l;
    int npts=kdTree->nPoints();
    double tmp, prod, mean; 
    int ind_k, k_max=k_vec.A[k_vec.ind_max-1];  // !!! k must be sorted, we take the largest

    std::vector<double> mom;                    // 2024/10/28: to optimize computations, C++ allocation 
    mom.resize(nA*(order_max+1));               // 2025/01/14: to compute central moments, we include moment of order 0
    for (d=0; d<nA*(order_max+1); d++)  mom[d]=0.;

    kdTree->annkSearch(x,                       // query point
                       k_max+ANN_ALLOW_SELF_MATCH,  // number of near neighbors (including or excluding central point)
                       nnIdx[core],             // nearest neighbors (returned)
                       dists[core],             // distance (returned)
                       ANN_eps);

    N_old=0;
    for (ind_k=k_vec.ind_min; ind_k<k_vec.ind_max; ind_k++)
    {   N=k_vec.A[ind_k]-1+ANN_ALLOW_SELF_MATCH;
        if (R!=NULL) R[npts_out*ind_k] = (double)dists[core][N-1];

        for (d=0; d<nA; d++)
        {   
            for (i=N_old; i<N; i++) 
            {   tmp        = (A+npts*d)[nnIdx[core][i]]; 
                prod       = 1.;
                for (j_moments=1; j_moments<=order_max; j_moments++)
                {   prod *= tmp;
                    mom[d + j_moments*nA] += prod;
                }
            }

            if (do_center==0)       // natural moments (moments about the origin):
            {   for (j_moments=0; j_moments<order_max; j_moments++)
                {   moments[npts_out*(nA*(k_vec.N*j_moments + ind_k) + d)] = mom[d + (j_moments+1)*nA]/N;
                }
            }
            else                    // central moments: https://en.wikipedia.org/wiki/Central_moment
            {   
                if (order_max>0)    
                {   mom[d] = (double)N;             // convention for moment of order 0
                    mean   = mom[d + 1*nA] /N;         
                    moments[npts_out*(nA*(k_vec.N*0 + ind_k) + d)] = mom[d];    // mean of centered data is 0, but we return the real mean (2025-01-24)
                }
                for (j_moments=1; j_moments<order_max; j_moments++)
                {   moments[npts_out*(nA*(k_vec.N*j_moments + ind_k) + d)]  = 0.0;
                    for (l=0; l<=j_moments+1; l++)
                    {   moments[npts_out*(nA*(k_vec.N*j_moments + ind_k) + d)] += get_binomial(j_moments+1)[l] * mom[d + l*nA]/N * pow(-mean, j_moments+1-l);
                    }
                }
            }
        }

        N_old=N;
    }

    return((double)dists[core][N-1]);
} /* end of function "ANN_compute_stats_multi_k" ***********************************************/
