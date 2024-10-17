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
    double ANN_compute_stats_single_k(double *x, double *A, int k,          double *R, double *mean, double *var, int npts_out, int nA, int core); // new 2024
    double ANN_compute_stats_multi_k (double *x, double *A, int *k, int Nk, double *R, double *mean, double *var, int npts_out, int nA, int core); // new 2024

#ifdef __cplusplus
}
#endif


#endif