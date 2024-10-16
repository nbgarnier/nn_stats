#cython: language_level=3
# Cython wrappers for C functions
#
# 2024-10-08 fork from entropy library

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

import  cython
import  numpy as PNP # shoud be useless here
cimport numpy as CNP
cimport commons
cimport nn_statistics
CNP.import_array()

include "commons.pyx"   # for basic library manipulation

@cython.boundscheck(False)
@cython.wraparound(False)
def compute_local_stats( double[:, ::1] x, double[:, ::1] A, double [:, ::1] y, int k=-1, double R=-1):
    """     
    compute local averages (and corresponding stds) of observables A (possibly multi-dimensional)
    given at locations x (usually 2-dimensional).
    
    The averages are computed at new locations y.   
    
    :param x: initial locations/positions (NumPy array with ndim=2, all coordinates along 1st dimension)
    :param A: embedding dimension :math:`m` (default=1)
    :param y: locations/positions (NumPy array with ndim=2, all coordinates along 1st dimension) where statistics will be computed.
    :param k: number of neighbors to consider for a fixed-k computation.
    :param R: radius to consider for a fixed-radius computation.
                 
    :returns: the local averages of A over y, using either fixed-k or fixed-R.
    """
    
    cdef int npts_in =x.shape[1], nx=x.shape[0], ratou # 2018-04-13: carefull with ordering of dimensions!
    cdef int npts_A  =A.shape[1], nA=A.shape[0]
    cdef int npts_out=y.shape[1], ny=y.shape[0]
    
    cdef CNP.ndarray[dtype=double, ndim=2, mode='c'] dists  #= void;
    cdef CNP.ndarray[dtype=int,    ndim=2, mode='c'] nnn
    cdef CNP.ndarray[dtype=double, ndim=2, mode='c'] A_mean = PNP.zeros((nA,npts_out), dtype=PNP.float64) 
    cdef CNP.ndarray[dtype=double, ndim=2, mode='c'] A_var  = PNP.zeros((nA,npts_out), dtype=PNP.float64) 
     
    if (npts_in<nx):      raise ValueError("please transpose x")
    if (npts_A<nA):       raise ValueError("please transpose A")
    if (npts_out<ny):     raise ValueError("please transpose y")
    if (npts_A!=npts_in): raise ValueError("A and x do not have the same nb of points!")
    if (ny!=nx):          raise ValueError("y and x do not have the same dimensionality!")
    if (k<=1) and (R<=0): raise ValueError("specify at least k or radius R!")
    if (k>1) and (R>0):   raise ValueError("specify either k or radius R, but not both!")
    
    if (k>1):   
#                print("fixed k computation")
#                dists = PNP.zeros((k,npts_out), dtype=PNP.float64) # final version
        dists = PNP.zeros((1,npts_out), dtype=PNP.float64) # tmp version
        ratou = nn_statistics.compute_stats_fixed_k_threads(&x[0,0], &A[0,0], npts_in, nx, nA, &y[0,0], npts_out, k, &A_mean[0,0], &A_var[0,0], &dists[0,0])
        return  PNP.asarray(A_mean), PNP.asarray(A_var), PNP.asarray(dists)
    if (R>0):   
        nnn   = PNP.zeros((1,npts_out), dtype=PNP.intc) # tmp version
        #ratou = nn_statistics.compute_stats_fixed_R_threads(&x[0,0], &A[0,0], npts_in, nx, nA, &y[0,0], npts_out, R, &A_mean[0,0], &A_var[0,0], &nnn[0,0])
        print(A_mean)
        print(nnn)
        return  PNP.asarray(A_mean), PNP.asarray(A_var), PNP.asarray(nnn)

