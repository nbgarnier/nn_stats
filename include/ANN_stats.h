//----------------------------------------------------------------------
// File:			ANN_stats.h
//		
// Description:	header file for ANN_wrapper.c, to use ANN in C
// This is a C header, and ANN_wrapper.cpp is a cpp file
//
// 2021-12-01 : added parameters "core" and "nb_cores" for multithreading
//----------------------------------------------------------------------
#ifndef ANN_STATS_H
#define ANN_STATS_H

#ifdef __cplusplus
extern "C" {
#endif
    double ANN_compute_stats            (double *x, double *A, int k, double *R, double *mean, double *var, int npts_out, int nA, int core); // new 2024

#ifdef __cplusplus
}
#endif


#endif
