//----------------------------------------------------------------------
// File:			ANN_wrapper.cpp
//		
// Description:	header file for ANN_wrapper.c, to use ANN in C
// This is a C++ header
//
// lines 84 and 89 are important depending on the algo choice
//
// 2017-11-29: added function "search_ANN_external" (for relative entropy)
// 2017-11-29: to-do: simplify parameters call of "search_ANN"
// 2017-11-29: to-do: use global variable "dists" instead of re-computing epsilon in "search_ANN*"
//
// 2021-11-29: multithread adaptations started...
// 2021-12-01: multithread entropy OK
// 2024-10-07: new functions for local averaging ("nn_stats")
//----------------------------------------------------------------------
#define noDEBUG
#define DEBUG_N 37
 
#include "ANN/ANN.h"
#include "ANN_wrapper.h"        // definitions of functions only
//#include <stdio.h>              // for printf, to be removed
//#include <iostream>

#define UNUSED(expr) do { (void)(expr); } while (0)
// https://stackoverflow.com/questions/1486904/how-do-i-best-silence-a-warning-about-unused-variables

// #include "ann_1.1.2/src/pr_queue_k.h"	// 2021-12-01, k-element priority queue, for definition of ANNmin_k

// global variables for (internal, but lower level so "private") operations on the main tree:
// note that for these variables, moved from other .cpp files, the names have been kept for
// consistency accross the library

// global variables for the main tree, internal to the C++ code, but exposed only in ANN_Wrapper.cpp:
double	        ANN_eps = 0.f;  // error bound, exact search if 0
ANNpointArray	dataPts;        // data points, type (*(*double))
ANNidxArray	   *nnIdx;          // k-nn indices    // 2021, adapted for pthread 
ANNdistArray   *dists;          // k-nn distances  // 2021, adapted for pthread
ANNkd_tree*	    kdTree;         // search structure
 
// allocate function (housekeeping) :
int init_ANN(int maxPts, int dim, int max_k, int NCORES=1)
{   if (NCORES<1) return(-1);

 	dataPts = annAllocPts(maxPts, dim);			// allocate data points 	
	nnIdx   = new ANNidx*[NCORES];				// allocate nearest neighbors indices 
	for (int i=0; i<NCORES; i++) nnIdx[i] = new ANNidx[max_k+1];
	dists   = new ANNdist*[NCORES];             // allocate nearest neighbors distances
	for (int i=0; i<NCORES; i++) dists[i] = new ANNdist[max_k+1];
				
//    std::cout << "[init_ANN] OK, using " << NCORES <<" core(s)\n";
    return(0);
}


// following two functions are new, to set the ANN_ALLOW_SELF_MATCH variable/singleton
// new 2014/02
// commented out 2014-06-03
/*
int set_ANN_state(char CHOICE)
{   int tmp;

    tmp = set_ANN_ALLOW_SELF_MATCH((char)CHOICE);
    std::cout << "set to " << tmp <<"\n";
    return(set_ANN_ALLOW_SELF_MATCH((char)CHOICE));
}

int get_ANN_state(void)
{   int tmp = get_ANN_ALLOW_SELF_MATCH();
    std::cout << "got " << tmp << "\n";
    return(tmp);
}
*/


// to do in create_kd_tree()
// dataPts must be allocated !!!
int create_kd_tree(double *x, int npts, int n)
{   int i,j;

    for (i=0; i<npts; i++)
    for (j=0; j<n;  j++)
        dataPts[i][j] = x[i + j*npts]; // dataPts is an ANNpointArray, so an array of ANNpoints
        
    // build search structure :
    kdTree = new ANNkd_tree(dataPts, npts, n);
    return(0);
}


/****************************************************************************************/
/* below : piece of code to search for nearest neighbors of an external point           */
/*         using a previously computed kd-tree (with ANN library)                       */
/*         This function returns the marginal distances epsilon_z of the                */
/*         corresponding ball, to be used by the Kraskov et al. "second" algorithm      */
/*                                                                                      */
/* input parameters:                                                                    */
/* x    : is the query point (the one to search neighbors of)                           */
/* n    : is the dimension of vector space (dimension of x)                             */
/* k    : is the rank of the neighbor to search for                                     */
/* core : index of thread or core on which to operate                                   */
/* output parameters:                                                                   */
/* epsilon_z : n-dimensional vector of maximal 1-d distances from x                     */
/*                                                                                      */
/* 2019-01-28: first version                                                            */
/* 2019-12-17: unchanged, but tested against "search_ANN_external": much better!        */
/* 2021-12-02: threaded version (parameter "core" introduced)                           */
/****************************************************************************************/
int ANN_marginal_distances_ex(double *x, int n, int k, double *epsilon_z, int core)
{   double eps_local;
    int d, l;
    
    kdTree->annkSearch(                     // search
                     x,                     // query point
                     k+ANN_ALLOW_SELF_MATCH,// number of near neighbors (including or excluding central point)
                     nnIdx[core],           // nearest neighbors (returned)
                     dists[core],           // distance (returned)
                     ANN_eps,
                     core);
    
    for (d=0; d<n; d++)  // for each dimension
    {   epsilon_z[d] = 0.0;
        for (l=0; l<k+ANN_ALLOW_SELF_MATCH; l++) // loop over the neighbors
        {   // nnIdx[l] is the index of the l-th nearest neighbor
            eps_local = fabs(x[d] - dataPts[nnIdx[core][l]][d]); // distance du l-ieme voisin, dans la direction d
            if (eps_local>epsilon_z[d]) epsilon_z[d] = eps_local; // search the max in this direction
        }
    }
    // 2019-01-29: at this stage, there is a different epsilon in each dimension
    return(0);
} /* end of function "ANN_marginal_distances_ex" ****************************************/



/***************************************************************************************/
/* below : piece of code to search for nearest neighbors of a point                    */
/*         using a previously computed kd-tree (with ANN library)                      */
/* faster and memory efficient coding                                                  */
/*                                                                                     */
/* input parameters:                                                                   */
/* i    : is the index of the query point (the one to search neighbors of)             */
/* n    : is the dimension of vector space (dimension of x)                            */
/* k    : is the rank of the neighbor to search for                                    */
/* output parameters:                                                                  */
/* (returned) : distance of the k-nn neighbor from x                                   */
/*                                                                                     */
/* 2019-01-22 - first version                                                          */
/* 2021-12-01 - thread version                                                         */
/***************************************************************************************/
double ANN_find_distance_in(int i, int n, int k, int core) // 2021-12-01: core=0 if not multithread!
{   UNUSED(n);
    kdTree->annkSearch(            // search
                       dataPts[i], //queryPt,            // query point
                       k+ANN_ALLOW_SELF_MATCH,  // number of near neighbors (including or excluding central point)
                       nnIdx[core],                // nearest neighbors (returned)
                       dists[core],                // distance (returned)
                       ANN_eps,
                       core);
// 2021-12-13, test   
#ifdef DEBUG 
    if (i==DEBUG_N)
    {   std::cout << "\t[ANN_find_distance_in] with k=" << k+ANN_ALLOW_SELF_MATCH << "\n";
        std::cout << "\tnnIdx = [ ";
        for (int q=0; q<k+ANN_ALLOW_SELF_MATCH; q++) std::cout << nnIdx[core][q] << " ";
        std::cout << "]\n\tdists = [ ";
        for (int q=0; q<k+ANN_ALLOW_SELF_MATCH; q++) std::cout << dists[core][q] << " ";
        std::cout << "]\n";
    }
#endif    
    return( ((nnIdx[core][k-1+ANN_ALLOW_SELF_MATCH]<0) ? 0. : (double)dists[core][k-1+ANN_ALLOW_SELF_MATCH]) ); 
    // distance from central point
} /* end of function "ANN_find_distance_in" ****************************************/



/***************************************************************************************/
/* below : piece of code to search for nearest neighbors of a point                    */
/*         using a previously computed kd-tree (with ANN library)                      */
/* faster and memory efficient coding                                                  */
/*                                                                                     */
/* input parameters:                                                                   */
/* x    : is a d-dimensional query point (the one to search neighbors of)              */
/* n    : is the dimension of vector space (dimension of x)                            */
/* k    : is the rank of the neighbor to search for                                    */
/* output parameters:                                                                  */
/* (returned) : distance of the k-nn from x                                            */
/*                                                                                     */
/* 2019-01-22 - first version                                                          */
/***************************************************************************************/
double ANN_find_distance_ex(double *x, int n, int k, int core)
{   UNUSED(n);
    //int d;
    /* composantes du point central : */
    //for (d=0; d<n; d++) queryPt[d] = x[d];
    
    kdTree->annkSearch(x, // queryPt,            // query point
                       k+ANN_ALLOW_SELF_MATCH,  // number of near neighbors (including or excluding central point)
                       nnIdx[core],                // nearest neighbors (returned)
                       dists[core],                // distance (returned)
                       ANN_eps);

    return((double)dists[core][k-1+ANN_ALLOW_SELF_MATCH]); // distance from central point
} /* end of function "ANN_find_distance_ex" ****************************************/



/**************************************************************************************
 * ANN_count_nearest_neighbors
 * to count nearest neighbors of a point x0 in a ball of radius epsilon in n-dimensions
 *
 * x0 is the central point (of dimension d, the dimension of the tree)
 * epsilon is the ball radius
 *
 * 2019-01-29 : first version
 * 2021-12-02 : multithread-safe versions
 * 2024-10-07 : adapted version
 *************************************************************************************/
int ANN_count_nearest_neighbors(double *x, double epsilon, int core)
{   return(kdTree->annkFRSearch(x,      // query point
                        epsilon,        // squared radius (same as radius for L^\infty norm)
                        0,              // (number of near neighbors to return), k=0 to search and count
                        NULL, // nnIdx_1[core],  // nearest neighbors (returned if !=NULL)
                        NULL, // dists_1[core],  // distance (returned if !=NULL)
                        ANN_eps,
                        core)); 
}



// free_function (housekeeping) :
void free_ANN(int NCORES)
{   
//    delete [] ANNkdDim;
//    annDeallocPts(ANNkdQ);
//    delete [] ANNkdPointMK;
    
    for (int i=0; i<NCORES; i++) { delete nnIdx[i]; }   delete [] nnIdx;
    for (int i=0; i<NCORES; i++) { delete dists[i]; }   delete [] dists;
    annDeallocPts(dataPts);
    delete kdTree;             // clean main tree
    
    annClose();                // done with ANN
}


// a quick fix to a compiler/linker bug:
// https://stackoverflow.com/questions/75139508/c-extension-for-python-throws-symbol-not-found-in-flat-namespace-zl10func1
#include "ANN_stats.cpp"
#include "ANN_stats_kernel.cpp"
