/*
 *  nn_stats.h
 *  
 *  Created by Nicolas Garnier on 2024/10/07 (fork from the entropy library).
 *  Copyright 2012-2024 ENS-Lyon - CNRS. All rights reserved.
 *
 *  The method used in these functions is from Grassberger (2004), 
 *  using nearest neighbor statistics and ANN library 
 *  http://www.cs.umd.edu/~mount/ANN/
 *
 *  This set of functions operates on continuous data !!!
 *
 *  2024-10-07 : created for "nn_stats.c", itself forked from "entropy_ann.c"
 */

double compute_stats_fixed_k(double *x, double *A, int npts, int nx, int nA, int k,    double *A_mean, double *A_var);
double compute_stats_fixed_R(double *x, double *A, int npts, int nx, int nA, double R, double *A_mean, double *A_var);


