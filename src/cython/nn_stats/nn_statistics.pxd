# file nn_statistics.pxd
#
# 2024-10-08

cdef extern from "nn_stats_fixed_k_threads.h":
    int compute_stats_fixed_k_threads(double *x, double *A, int npts_in, int nx, int nA, double *y, int npts_out, int k, double *A_mean, double *A_std, double *dists)

cdef extern from "nn_stats_fixed_R_threads.h":
    int compute_stats_fixed_R_threads(double *x, double *A, int npts_in, int nx, int nA, double *y, int npts_out, double R, double *A_mean, double *A_std, int *k)

	