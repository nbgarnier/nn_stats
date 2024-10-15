/*
 * data.h
 *
 * prototypes for some data struct
 * N.B.Garnier 2024-10-15
 * 
 */

struct process 
    {   int Npts;   // nb of points
        int dim;    // nb of observables
        double *A;  // observables data
    };
typedef struct process array;

struct discrete_process 
    {   int Npts;   // nb of points
        int dim;    // nb of observables
        int *A;     // observables data
    };
typedef struct discrete_process arr_int;

