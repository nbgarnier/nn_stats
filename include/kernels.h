// kernels.h
//
// definitions of kernels for local averagings
//
// 2025-01-20, N.B.G.

//#pragma GCC diagnostic push
//#pragma GCC diagnostic ignored "-Wgnu-statement-expression"
//#pragma GCC diagnostic pop

#ifndef KERNELS_H 
#define KERNELS_H

#define KERNEL_BRICKWALL 0x0001
#define KERNEL_GAUSSIAN  0x0002
#define KERNEL_TRIANGLE  0x0003

double kernel_brickwall(double x, double d);
double kernel_triangle (double x, double d);
double kernel_Gaussian (double x, double d);

// global variables : the kernel in use, and the (bservation) scale for the kernel 
// (all are defined in "kernel.c")
extern int    current_kernel_type;
extern double (*current_kernel)(double, double); //= kernel_brickwall;    // default kernel
extern double current_obs_scale;

void select_kernel(int kernel_type, double prescribed_scale);

#endif
