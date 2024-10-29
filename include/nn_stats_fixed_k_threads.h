/************************************************************************/
/* nn_stats_fixed_k_threads.h			                                */
/*									                                    */
/* 26/11/2021 	                                                        */
/************************************************************************/
/* Nicolas Garnier	nicolas.garnier@ens-lyon.fr			                */
/************************************************************************/

int compute_stats_fixed_k_threads(double *x, double *A, int npts, int nx, int nA, double *y, int npts_out, int k,          double *A_mean, double *A_std, double *dists);

int compute_stats_multi_k_threads(double *x, double *A, int npts, int nx, int nA, double *y, int npts_out, int *k, int nk, double *A_mean, double *A_std, double *dists);

