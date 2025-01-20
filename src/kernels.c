// kernels.c
//
// definitions of kernels for local averagings
//
// 2025-01-20, N.B.G.

#include <math.h>
#include "kernels.h"


double kernel_brickwall(double x, double d) 
{   return 1.0 + 0.*(x - d);
}

double kernel_triangle(double x, double d) 
{   return 1-x/d;
}

double kernel_Gaussian(double x, double d) 
{   return exp(-(x/d)*(x/d)/2);
}

// global variable : the kernel in use, defined in kernel.c
int    current_kernel_type = 0;
double (*current_kernel)(double, double) = kernel_brickwall;    // default kernel
double obs_scale = 1.;                                          // default observation scale for the kernel

void select_kernel(int kernel_type, double prescribed_scale)
{   switch (kernel_type)
    {   case 0 : 
            current_kernel_type = 0;
            current_kernel      = &kernel_brickwall;            // select regular kernel
        break;
        case 1 : 
            current_kernel_type = 1;
            current_kernel      = &kernel_triangle;             // select triangle kernel
        break;
        case 2 : 
            current_kernel_type = 2;
            current_kernel      = &kernel_Gaussian;             // select Gaussian kernel
        break;
        default: 
            current_kernel_type = 0;
            current_kernel      = &kernel_brickwall;            // select regular kernel
        break;  
    }
    if (prescribed_scale>0) obs_scale = prescribed_scale;
}
