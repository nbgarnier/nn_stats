// kernels.c
//
// definitions of kernels for local averagings
//
// 2025-01-20, N.B.G.

#include <math.h>
#include "kernels.h"


double kernel_brickwall(double x, double d) 
{   (void)x; (void)d;
    return 1.;
}

double kernel_Gaussian(double x, double d) 
{   return exp(-(x/d)*(x/d)/2);
}

double kernel_triangle(double x, double d) 
{   if (fabs(x/d) <1.) return 1-x/d;
    else return 0.;
}

double kernel_quartique(double x, double d) 
{   double a=kernel_Epanechnikov(x, d);
    return a*a;
}

double kernel_Epanechnikov(double x, double d) 
{   if (fabs(x/d) <1.) return 1.-(x/d)*(x/d);
    else return 0.;
}

double kernel_exponential(double x, double d) 
{   return exp(-(x/d));
}

// global variables : the kernel in use, defined in kernel.c, and related parameters
int    current_kernel_type = 0;
double (*current_kernel)(double, double) = kernel_brickwall;    // default kernel
double current_obs_scale   = 1.;                                // default observation scale for the kernel

void select_kernel(int kernel_type, double prescribed_scale)
{   switch (kernel_type)
    {   case 0 :                                                // no kernel (faster functions)
        case 1 :                                                // regular kernel
            current_kernel_type = kernel_type;                  // we keep the distinction between 0 (no kernel) and 1 (regular kernel)
            current_kernel      = &kernel_brickwall;            // select regular kernel
        break;
        case 2 : 
            current_kernel_type = 2;
            current_kernel      = &kernel_Gaussian;             // select Gaussian kernel
        break;
        case 3 : 
            current_kernel_type = 3;
            current_kernel      = &kernel_triangle;             // select triangle kernel
        break;
        case 4 : 
            current_kernel_type = 4;
            current_kernel      = &kernel_quartique;            // select quartique kernel
        break;
        case 5 : 
            current_kernel_type = 5;
            current_kernel      = &kernel_Epanechnikov;          // select Epanechikov kernel
        break;
        case 6 : 
            current_kernel_type = 6;
            current_kernel      = &kernel_exponential;          // select exponential kernel
        break;
        default: 
            current_kernel_type = 0;                            // no kernel
            current_kernel      = &kernel_brickwall;            // select regular kernel
        break;  
    }
    if (prescribed_scale>0) current_obs_scale = prescribed_scale;
}
