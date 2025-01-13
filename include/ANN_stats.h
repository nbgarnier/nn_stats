//----------------------------------------------------------------------
// File:			ANN_stats.h
//		
// Description:	header file for ANN_stats.cpp
//
//----------------------------------------------------------------------
#ifndef ANN_STATS_H
#define ANN_STATS_H

#ifdef __cplusplus
extern "C" {
#endif

    struct int_range
    {   int ind_min;    // min index
        int ind_max;    // max index
        int N;          // total nb of index
        int *A;         // index set
    };
    typedef struct int_range k_vector;

    double ANN_compute_stats_single_k(double *x, double *A, int k,          double *R, double *moments, int order_max, int npts_out, int nA, int core);

    double ANN_compute_stats_multi_k (double *x, double *A, k_vector k_vec, double *R, double *moments, int order_max, int npts_out, int nA, int core);

#ifdef __cplusplus
}
#endif


#endif
