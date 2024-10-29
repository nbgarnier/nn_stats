/************************************************************************/
/* nn_stats_fixed_k_threads.h			                                */
/*									                                    */
/* 26/11/2021 	                                                        */
/************************************************************************/
/* Nicolas Garnier	nicolas.garnier@ens-lyon.fr			                */
/************************************************************************/

int compute_stats_fixed_R_threads(double *x, double *A, int npts_in, int nx, int nA, double *y, int npts_out, double R,          double *A_mean, double *A_std, int *k);

int compute_stats_fixed_R_threads(double *x, double *A, int npts_in, int nx, int nA, double *y, int npts_out, double *R, int nR, double *A_mean, double *A_std, int *k);

