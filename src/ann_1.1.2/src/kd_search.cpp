//----------------------------------------------------------------------
// File:			kd_search.cpp
// Programmer:		Sunil Arya and David Mount
// Description:		Standard kd-tree search
// Last modified:	01/12/21 (Version 1.0)
//----------------------------------------------------------------------
// Copyright (c) 1997-2005 University of Maryland and Sunil Arya and
// David Mount.  All Rights Reserved.
// 
// This software and related documentation is part of the Approximate
// Nearest Neighbor Library (ANN).  This software is provided under
// the provisions of the Lesser GNU Public License (LGPL).  See the
// file ../ReadMe.txt for further information.
// 
// The University of Maryland (U.M.) and the authors make no
// representations about the suitability or fitness of this software for
// any purpose.  It is provided "as is" without express or implied
// warranty.
//----------------------------------------------------------------------
// History:
//	Revision 0.1  03/04/98
//		Initial release
//	Revision 1.0  04/01/05
//		Changed names LO, HI to ANN_LO, ANN_HI
//  Modifications 01/12/21 - Nicolas B. Garnier
//      made procedures thread-safe by isolating global variables per thread
//----------------------------------------------------------------------

#include "kd_search.h"					// kd-search declarations

#define UNUSED(expr) do { (void)(expr); } while (0)
// https://stackoverflow.com/questions/1486904/how-do-i-best-silence-a-warning-about-unused-variables

//----------------------------------------------------------------------
//	Approximate nearest neighbor searching by kd-tree search
//		The kd-tree is searched for an approximate nearest neighbor.
//		The point is returned through one of the arguments, and the
//		distance returned is the squared distance to this point.
//
//		The method used for searching the kd-tree is an approximate
//		adaptation of the search algorithm described by Friedman,
//		Bentley, and Finkel, ``An algorithm for finding best matches
//		in logarithmic expected time,'' ACM Transactions on Mathematical
//		Software, 3(3):209-226, 1977).
//
//		The algorithm operates recursively.  When first encountering a
//		node of the kd-tree we first visit the child which is closest to
//		the query point.  On return, we decide whether we want to visit
//		the other child.  If the box containing the other child exceeds
//		1/(1+eps) times the current best distance, then we skip it (since
//		any point found in this child cannot be closer to the query point
//		by more than this factor.)  Otherwise, we visit it recursively.
//		The distance between a box and the query point is computed exactly
//		(not approximated as is often done in kd-tree), using incremental
//		distance updates, as described by Arya and Mount in ``Algorithms
//		for fast vector quantization,'' Proc.  of DCC '93: Data Compression
//		Conference, eds. J. A. Storer and M. Cohn, IEEE Press, 1993,
//		381-390.
//
//		The main entry points is annkSearch() which sets things up and
//		then call the recursive routine ann_search().  This is a recursive
//		routine which performs the processing for one node in the kd-tree.
//		There are two versions of this virtual procedure, one for splitting
//		nodes and one for leaves.  When a splitting node is visited, we
//		determine which child to visit first (the closer one), and visit
//		the other child on return.  When a leaf is visited, we compute
//		the distances to the points in the buckets, and update information
//		on the closest points.
//
//		Some trickery is used to incrementally update the distance from
//		a kd-tree rectangle to the query point.  This comes about from
//		the fact that which each successive split, only one component
//		(along the dimension that is split) of the squared distance to
//		the child rectangle is different from the squared distance to
//		the parent rectangle.
//----------------------------------------------------------------------

//----------------------------------------------------------------------
//		To keep argument lists short, a number of global variables
//		are maintained which are common to all the recursive calls.
//		These are given below. 
// 
// 		2021-12-01: unfortunately, these global declarations make the 
// 		search being not thread-safe! 
//		Tentative corrections are started to be applied... 
//		which led to a modification of global variables
//		(some moved to "ANN_wrapper.cpp")
//----------------------------------------------------------------------

__thread int	    ANNkdDim;	    // dimension of space (untouched, so thread-safe)
__thread ANNpoint   ANNkdQ;	   	    // query point // 2021-12-01: adapted for pthread (made thread-safe)
									// needs to be allocated: this is done by init_ANN in ANN_wrapper.cpp
double			ANNkdMaxErr;		// max tolerable squared error (untouched, so thread-safe)
ANNpointArray	ANNkdPts;			// the points (untouched, so thread-safe)
__thread ANNmin_k   *ANNkdPointMK;	// set of k closest points // 2021-12-01: adapted for pthread (made thread-safe)
									// needs to be allocated: this is done by init_ANN 
									// (now defined for real in ANN_wrapper.cpp)

//----------------------------------------------------------------------
//	annkSearch - search for the k nearest neighbors
//  2021-12-01 : now thread-safe, using specified core index
//----------------------------------------------------------------------

void ANNkd_tree::annkSearch(
	ANNpoint			q,			// the query point
	int					k,			// number of near neighbors to return
	ANNidxArray			nn_idx,		// nearest neighbor indices (returned)
	ANNdistArray		dd,			// the approximate nearest neighbor
	double				eps,		// the error bound
	int					core)		// added 2021-12-01: an index to make call thread-safe
{	
	ANNkdDim = dim;			// copy arguments to static equivs
	ANNkdQ = q;
	ANNkdPts = pts;
	ANNptsVisited = 0;				// initialize count of points visited

	if (k > n_pts) {				// too many near neighbors?
		annError("Requesting more near neighbors than data points", ANNabort);
	}

	// DEBUG CODE:
//	std::cout << "-\tdim of space:" << ANNkdDim[core] << "\t max error " << ANNkdMaxErr << "\n";
//	std::cout << "- data:\t" << ANNkdPts[0][0] << " " << ANNkdPts[1][0] << " " << ANNkdPts[2][0] << "\n" ;
//	std::cout << "- query pt:";
//	for (register int j=0; j<dim; j++) std::cout << "\t" << ANNkdQ[core][j];
//	std::cout << "\n";

	ANNkdMaxErr = ANN_POW(1.0 + eps);
	ANN_FLOP(2)						// increment floating op count
	
// std::cout << "[kdTree::annkSearch] - working in core " << core << "\n";
	ANNkdPointMK = new ANNmin_k(k);		// create set for closest k points
// std::cout << "[kdTree::annkSearch] - allocation of pointer OK\n";	
									// search starting at the root
	root->ann_search(annBoxDistance(q, bnd_box_lo, bnd_box_hi, dim), core);

	for (int i = 0; i < k; i++) {	// extract the k-th closest points
		dd[i] = ANNkdPointMK->ith_smallest_key(i);
		nn_idx[i] = ANNkdPointMK->ith_smallest_info(i);
	}
// debug trace:
//	if (core>0) 
//	{	for (int i = 0; i < k; i++) std::cout << "\t" << dd[i] << " " << nn_idx[i];		
//		std::cout << "\n";	
//	}
	
	delete ANNkdPointMK;		// deallocate closest point set
}

//----------------------------------------------------------------------
//	kd_split::ann_search - search a splitting node
//----------------------------------------------------------------------

void ANNkd_split::ann_search(ANNdist box_dist, int core)
{
										// check dist calc term condition
	if (ANNmaxPtsVisited != 0 && ANNptsVisited > ANNmaxPtsVisited) return;

										// distance to cutting plane
	ANNcoord cut_diff = ANNkdQ[cut_dim] - cut_val;

	if (cut_diff < 0) {					// left of cutting plane
		child[ANN_LO]->ann_search(box_dist, core);// visit closer child first

		ANNcoord box_diff = cd_bnds[ANN_LO] - ANNkdQ[cut_dim];
		if (box_diff < 0)				// within bounds - ignore
			box_diff = 0;
										// distance to further box
		box_dist = (ANNdist) ANN_SUM(box_dist,
				ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

										// visit further child if close enough
		if (box_dist * ANNkdMaxErr < ANNkdPointMK->max_key())
			child[ANN_HI]->ann_search(box_dist, core);

	}
	else {								// right of cutting plane
		child[ANN_HI]->ann_search(box_dist, core);// visit closer child first

		ANNcoord box_diff = ANNkdQ[cut_dim] - cd_bnds[ANN_HI];
		if (box_diff < 0)				// within bounds - ignore
			box_diff = 0;
										// distance to further box
		box_dist = (ANNdist) ANN_SUM(box_dist,
				ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

										// visit further child if close enough
		if (box_dist * ANNkdMaxErr < ANNkdPointMK->max_key())
			child[ANN_LO]->ann_search(box_dist, core);

	}
	ANN_FLOP(10)						// increment floating ops
	ANN_SPL(1)							// one more splitting node visited
}

//----------------------------------------------------------------------
//	kd_leaf::ann_search - search points in a leaf node
//		Note: The unreadability of this code is the result of
//		some fine tuning to replace indexing by pointer operations.
//----------------------------------------------------------------------

void ANNkd_leaf::ann_search(ANNdist box_dist, int core)
{
    ANNdist dist;				// distance to data point
    ANNcoord* pp;				// data coordinate pointer
    ANNcoord* qq;				// query coordinate pointer
    ANNdist min_dist;			// distance to k-th closest point
    ANNcoord t;
    int d;

	UNUSED(box_dist);
	UNUSED(core);

	min_dist = ANNkdPointMK->max_key(); // k-th smallest distance so far

	for (int i = 0; i < n_pts; i++) {	// check points in bucket

		pp = ANNkdPts[bkt[i]];			// first coord of next data point
		qq = ANNkdQ;					// first coord of query point // 2021-12-01, adapted for pthread
		dist = 0;

		for(d = 0; d < ANNkdDim; d++) {
			ANN_COORD(1)				// one more coordinate hit
			ANN_FLOP(4)					// increment floating ops

			t = *(qq++) - *(pp++);		// compute length and adv coordinate
										// exceeds dist to k-th smallest?
			if( (dist = ANN_SUM(dist, ANN_POW(t))) > min_dist) {
				break;
			}
		}

		if (d >= ANNkdDim &&				// among the k best?
		   (ANN_ALLOW_SELF_MATCH || dist!=0)) { // and no self-match problem
												// add it to the list
			ANNkdPointMK->insert(dist, bkt[i]);
			min_dist = ANNkdPointMK->max_key();
		}
	}
	ANN_LEAF(1)							// one more leaf node visited
	ANN_PTS(n_pts)						// increment points visited
	ANNptsVisited += n_pts;				// increment number of points visited
}
