//----------------------------------------------------------------------
// File:			ANN_stats_kernel.h
//		
// Description:	header file for ANN_stats_kernel.cpp
//
//----------------------------------------------------------------------
#ifndef ANN_STATS_KERNEL_H
#define ANN_STATS_KERNEL_H

#ifdef __cplusplus
extern "C" {
#endif

    double ANN_compute_stats_kernel_single_k(double *x, double *A, int k,          double *R, double *moments, int order_max, int npts_out, int nA, int do_center, int core);

    double ANN_compute_stats_kernel_multi_k (double *x, double *A, k_vector k_vec, double *R, double *moments, int order_max, int npts_out, int nA, int do_center, int core);

#ifdef __cplusplus
}
#endif


#endif
