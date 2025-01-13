# file nn_statistics.pxd
#
# 2024-10-08

cdef extern from "library_commons.h":
    int k_default
    int tree_k_max

cdef extern from "nn_stats_fixed_k_threads.h":
    int compute_stats_fixed_k_threads(double *x, double *A, int npts_in, int nx, int nA, double *y, int npts_out, int k,          double *A_moments, int order_max, double *dists)
    int compute_stats_multi_k_threads(double *x, double *A, int npts_in, int nx, int nA, double *y, int npts_out, int *k, int nk, double *A_moments, int order_max, double *dists)

cdef extern from "nn_stats_fixed_R_threads.h":
    int compute_stats_fixed_R_threads(double *x, double *A, int npts_in, int nx, int nA, double *y, int npts_out, double R,          double *A_moments, int order_max, int *k)
    int compute_stats_multi_R_threads(double *x, double *A, int npts_in, int nx, int nA, double *y, int npts_out, double *R, int nR, double *A_moments, int order_max, int *k)

	