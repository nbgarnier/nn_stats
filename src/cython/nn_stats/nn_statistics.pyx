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
def compute_local_stats( double[:, ::1] x, double[:, ::1] A, double [:, ::1] y, 
                            int[::1] k=PNP.zeros(shape=(1),dtype=PNP.intc), 
                            double[::1] R=PNP.zeros(shape=(1),dtype=PNP.float64)):
    """     
    compute local averages (and corresponding stds) of observables A (possibly multi-dimensional)
    given at locations x (usually 2-dimensional).
    
    The averages are computed at new locations y, using either:
    - the k-nn with prescribed k,
    - all points in a ball of radius R, with prescribed R.
    
    :param x: initial locations/positions (NumPy array with ndim=2, all coordinates along 1st dimension)
    :param A: observables (NumPy array with ndim=2) 
    :param y: locations/positions (NumPy array with ndim=2, all coordinates along 1st dimension) where statistics will be computed.
    :param k: 1d-array (int) of number of neighbors to consider for a fixed-k computation.
    :param R: 1d-array of radii to consider for a fixed-radius computation.
                 
    :returns: the local averages of A over y, using either fixed-k or fixed-R.
    """
    
    cdef int npts_in =x.shape[1], nx=x.shape[0], ratou # 2018-04-13: carefull with ordering of dimensions!
    cdef int npts_A  =A.shape[1], nA=A.shape[0]
    cdef int npts_out=y.shape[1], ny=y.shape[0]
    cdef int nb_k    =k.size
    cdef int nb_R    =R.size
    
    cdef CNP.ndarray[dtype=double, ndim=2, mode='c'] dists  #= void;
    cdef CNP.ndarray[dtype=int,    ndim=2, mode='c'] nnn
    cdef CNP.ndarray[dtype=double, ndim=2, mode='c'] A_mean 
    cdef CNP.ndarray[dtype=double, ndim=2, mode='c'] A_var 
     
    if (npts_in<nx):      raise ValueError("please transpose x")
    if (npts_A<nA):       raise ValueError("please transpose A")
    if (npts_out<ny):     raise ValueError("please transpose y")
    if (npts_A!=npts_in): raise ValueError("A and x do not have the same nb of points!")
    if (ny!=nx):          raise ValueError("y and x do not have the same dimensionality!")
    if ( (nb_k>1) or (k[0]>0) ) and ( (nb_R>1) or (R[0]>0) ):
                          raise ValueError("specify either k or radius R, but not both!")
    if (nb_k>1) and (PNP.min(PNP.diff(k))<0):   raise ValueError("k should be sorted (with increasing values)")
    if (nb_R>1) and (PNP.min(PNP.diff(R))<0):   raise ValueError("R should be sorted (with increasing values)")
    if (k[0]==0) and (R[0]==0):                 raise ValueError("specify at least k or radius R!")
    
    if (k[0]>0):
        if (k[nb_k-1]>=npts_in):    raise ValueError("imposed k is larger than the number of input points!")
        print("fixed k computation", end=" ")
        dists  = PNP.zeros((nb_k,npts_out), dtype=PNP.float64) 
        A_mean = PNP.zeros((nb_k*nA,npts_out), dtype=PNP.float64)
        A_var  = PNP.zeros((nb_k*nA,npts_out), dtype=PNP.float64) 
        if (nb_k==1):
            print("1 value of k :", k[0])
            ratou  = nn_statistics.compute_stats_fixed_k_threads(&x[0,0], &A[0,0], npts_in, nx, nA, &y[0,0], npts_out, k[0], &A_mean[0,0], &A_var[0,0], &dists[0,0])
        else:
            print("multiple values of k :", PNP.array(k))
            ratou  = nn_statistics.compute_stats_multi_k_threads(&x[0,0], &A[0,0], npts_in, nx, nA, &y[0,0], npts_out, &k[0], nb_k, &A_mean[0,0], &A_var[0,0], &dists[0,0])
        return  PNP.asarray(A_mean), PNP.asarray(A_var), PNP.sqrt(PNP.asarray(dists))
    if (nb_R>0):   
        print("fixed R computation", end=" ")
        nnn    = PNP.zeros((nb_R,npts_out), dtype=PNP.intc)
        A_mean = PNP.zeros((nb_R*nA,npts_out), dtype=PNP.float64)
        A_var  = PNP.zeros((nb_R*nA,npts_out), dtype=PNP.float64)    
        if (nb_R==1):
            print("1 value of R :", R[0]) 
            ratou = nn_statistics.compute_stats_fixed_R_threads(&x[0,0], &A[0,0], npts_in, nx, nA, &y[0,0], npts_out, R[0], &A_mean[0,0], &A_var[0,0], &nnn[0,0])
        else:
            print("multiple values of R :", PNP.array(R))
            ratou = nn_statistics.compute_stats_multi_R_threads(&x[0,0], &A[0,0], npts_in, nx, nA, &y[0,0], npts_out, &R[0], nb_R, &A_mean[0,0], &A_var[0,0], &nnn[0,0])
            return(PNP.asarray(A_mean).reshape(nb_R, nA, npts_out), PNP.asarray(A_var.reshape((nb_R, nA, npts_out))), PNP.asarray(nnn) )
        return  PNP.asarray(A_mean), PNP.asarray(A_var), PNP.asarray(nnn)



def set_nn_max(int k=nn_statistics.k_default):
    """
    sets the maximal value of allowed number of nearest neighbors
    
    :param k: an integer (default=5)
    :returns: no output
    """
    nn_statistics.tree_k_max=k
    
def get_nn():
    """
    gets the current maximal value of allowed number of nearest neighbors
    
    :param none:
    :returns: an integer
    """
    return(nn_statistics.tree_k_max)


