// kernels.c
//
// definitions of kernels for local averagings
//
// 2025-01-20, N.B.G.

#include <math.h>
#include "kernels.h"
#include "library_commons.h"    // for PI


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

double kernel_triweight(double x, double d) 
{   double a=kernel_Epanechnikov(x, d);
    return a*a*a;
}

double kernel_tricube(double x, double d) 
{   double u=x/d, a=0.;
    if (u<1.) 
    {   a=1.-u*u*u;
        return a*a*a;
    }
    else return 0.;
}

double kernel_cosine(double x, double d) 
{   double u=x/d;
    if (u<1.) return cos(PI*u/2);
    else return 0.;
}

double kernel_exponential(double x, double d) 
{   return exp(-(x/d));
}

double kernel_logistic(double x, double d) 
{   double u=x/d;
    return 1./(exp(u) + 2. + exp(-u));
}

double kernel_sigmoid(double x, double d) 
{   double u=x/d;
    return 1./(exp(u) + exp(-u));
}

double kernel_Silverman(double x, double d) 
{   double u=x/d/sqrt(2.);
    return (exp(-u) * sin(u+PI/4.));
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
            current_kernel      = &kernel_Epanechnikov;         // select Epanechikov kernel
        break;
        case 6 : 
            current_kernel_type = 6;
            current_kernel      = &kernel_triweight;         
        break;
        case 7 : 
            current_kernel_type = 7;
            current_kernel      = &kernel_tricube;        
        break;
        case 8 : 
            current_kernel_type = 8;
            current_kernel      = &kernel_cosine;        
        break;
        case 9 : 
            current_kernel_type = 9;
            current_kernel      = &kernel_exponential;          // select exponential kernel
        break;
        case 10 : 
            current_kernel_type = 10;
            current_kernel      = &kernel_logistic;             // select logistic kernel
        break;
        case 11 : 
            current_kernel_type = 11;
            current_kernel      = &kernel_sigmoid;              // select sigmoid kernel
        break;
        case 12 : 
            current_kernel_type = 12;
            current_kernel      = &kernel_Silverman;            // select Silverman kernel
        break;
        default: 
            current_kernel_type = 0;                            // no kernel
            current_kernel      = &kernel_brickwall;            // select regular kernel
        break;  
    }
    if (prescribed_scale>0) current_obs_scale = prescribed_scale;
}
