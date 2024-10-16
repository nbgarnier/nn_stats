/*
 * library_commons.c
 *
 *  Created by Nicolas Garnier on 2020-03-05.
 *  Copyright 2012-2022 ENS-Lyon - CNRS. All rights reserved.
 *
 *
 */
#include <math.h>       // for fabs
#include <float.h>      // for DBL_MIN (2020-07-17)
#include <stdio.h>      // for printf

#include "library_commons.h" // for definitions of nb_errors, and stds

// global variables that changes the behavior of the library
int lib_verbosity=1;    // <0 : no messages, even if an error is encountered (not recommended!)
                        // 0  : messages only if an error is encountered (default)
                        // 1  : important warnings only
                        // 2 or more: more and more warnings
int lib_warning_level=1; // 0 : no physical checks, code trusts user / user is responsible
                        // 1  : physical checks raises warnings 
                        // 2  : physical checks raises errors 
                        //      (1 or 2 helps the user, but limits advanced use of the library)

int tree_k_max=5;

/****************************************************************************************/
/* selects the verbosity level                                                          */
/****************************************************************************************/
void ANN_set_verbosity(int level)
{   lib_verbosity=level;
    return;
}


int is_equal(double x, double y)
{   const double epsilon=1e-7; // arbitrary
    return(fabs(x - y) <= (epsilon*fabs(x)));
}

int is_zero(double x)
{   const double epsilon=1e-15; // arbitrary
//    const double epsilon=2*DBL_MIN; // DBL_MIN is the min double value on the machine
    return(fabs(x) <= epsilon);
}

