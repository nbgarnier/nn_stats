/************************************************************************/
/* nn_stats_fixed_k_threads.h			                                */
/*									                                    */
/* 26/11/2021 	                                                        */
/************************************************************************/
/* Nicolas Garnier	nicolas.garnier@ens-lyon.fr			                */
/************************************************************************/

int compute_stats_fixed_k_threads(double *x, double *A, int npts_in, int nx, int nA, double *y, int npts_out, int k,          double *A_moments, int order_max, int do_center, double *dists);

int compute_stats_multi_k_threads(double *x, double *A, int npts_in, int nx, int nA, double *y, int npts_out, int *k, int nk, double *A_moments, int order_max, int do_center, double *dists);

