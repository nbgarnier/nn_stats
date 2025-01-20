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

#define KERNEL_BRICKWALL 0x0000
#define KERNEL_TRIANGLE  0x0001
#define KERNEL_GAUSSIAN  0x0002

double kernel_brickwall(double x, double d);
double kernel_triangle (double x, double d);
double kernel_Gaussian (double x, double d);

// global variable : the kernel in use, defined in kernel.c
extern double (*current_kernel)(double, double); //= kernel_brickwall;    // default kernel

void select_kernel(int kernel_type);

#endif
