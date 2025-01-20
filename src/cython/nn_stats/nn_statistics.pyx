#cython: language_level=3
# Cython wrappers for C functions
#
# 2024-10-08 fork from entropy library

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

import  cython
import  numpy as PNP    # shoud be useless here
cimport numpy as CNP
cimport commons
cimport nn_statistics
import commons as PC
CNP.import_array()

include "commons.pyx"   # for basic library manipulation

@cython.boundscheck(False)
@cython.wraparound(False)
def compute_local_stats(double[:, ::1]  x, 
                        double[:, ::1]  A = PNP.zeros(shape=(0,1),dtype=PNP.float64), 
                        double [:, ::1] y = PNP.zeros(shape=(0,1),dtype=PNP.float64), 
                        int[::1]        k = PNP.zeros(shape=(1),dtype=PNP.intc), 
                        double[::1]     R = PNP.zeros(shape=(1),dtype=PNP.float64), 
                        int             order_max = 2,
                        bint            centered = False,
                        int             kernel_type = -1, 
                        double          scale = -1.,
                        int             nn_max = -1, 
                        int             verbosity = -1):
    """     
    compute local moments (mean, variance, etc.) of observables A (possibly multi-dimensional) initially 
    given at locations x (usually 2-dimensional).
    
    The averages are computed at new locations y, using either:
    - the k-nearest neighbors with prescribed k,
    - all points in a ball of radius R, with prescribed R.
    
    :param x: initial locations/positions (NumPy array with ndim=2, all coordinates along 1st dimension).
    :param A: observables (NumPy array with ndim=2) 
                If this parameter is not provided, no local average will be computed, but the distances (for fixed R)
                or the numbers of neighbors (for fixed k) will be returned.
    :param y: locations/positions (NumPy array with ndim=2, all coordinates along 1st dimension) where statistics will be computed. 
                If this parameter is not provided, the initial locations/positions will be used.
    :param k: 1d-array (int) of number of neighbors to consider for a fixed-k computation.
    :param R: 1d-array of radii to consider for a fixed-radius computation.
    :param order_max: maximal order of the moments to be computed (default=2)
    :param centered: Boolean to indicate if moments are centered (True) or not centered (False) (default=False)
    :param kernel_type: integer to indicate the kernel to use (default is set with the function 'set_kernel')
    :param scale: (observation) scale for the kernel (default is set with the fuunction 'set_kernel')
    :param nn_max: maximal nb of neighbors to consider when performing a fixed-R search (default : automatic, 10% of available points)
    :param verbosity: 0 to operate quietly without any message or larger value for more messages
                (default value can be set by function "set_verbosity")
                 
    :returns (1+order_max) variables: 
        - the set of R if k is provided OR the set of k if R is provided
        - the local moments of A (central or not) computed at y, using either fixed-k or fixed-R.
    """
    
    # analysing and formating input parameters
    cdef int npts_in =x.shape[1], nx=x.shape[0], i, ratou  # 2018-04-13: carefull with ordering of dimensions!
    if (A.shape[0]==0): 
        if (verbosity>1): print("A was not provided, shape", A.shape, end=' ')
        A=PNP.zeros( (0, npts_in))      # 2025-01-13: if A is not provided, we set it to appropriate shape
        if (verbosity>1): print(" -> new shape", A.shape)
    cdef int npts_A  =A.shape[1], nA=A.shape[0]
    if (y.shape[0]==0): 
        if (verbosity>1): print("y was not provided, shape", y.shape, end=' ')
        y=x.copy()                      # 2025-01-13: if destination locations are not provided, we use initial locations
        if (verbosity>1): print(" -> new shape", y.shape)
    cdef int npts_out=y.shape[1], ny=y.shape[0]
    cdef int nb_k    =k.size
    cdef int nb_R    =R.size
    cdef CNP.ndarray[dtype=double, ndim=1, mode='c'] R2     # squared radii

    if (verbosity>2):   print("function called with verbosity", verbosity, end=' ')
    if (verbosity==-1): verbosity=get_verbosity(0)
    if (verbosity>2):   print("so using verbosity", verbosity)
    
    if (verbosity>1): 
        print("k has shape", k.shape, "and R has shape", R.shape, "\t-> nb_k", nb_k, "  nb_R", nb_R)
    
    cdef CNP.ndarray[dtype=double, ndim=2, mode='c'] dists  #= void;
    cdef CNP.ndarray[dtype=int,    ndim=2, mode='c'] nnn
    cdef CNP.ndarray[dtype=double, ndim=2, mode='c'] moments 
     
    if (npts_in<nx):      raise ValueError("please transpose x")
    if (npts_A<nA):       raise ValueError("please transpose A")
    if (npts_out<ny):     raise ValueError("please transpose y")
    if (npts_A!=npts_in): raise ValueError("A and x do not have the same nb of points!")
    if (ny!=nx):          raise ValueError("y and x do not have the same dimensionality!")
    if (order_max<0):     raise ValueError("order_max cannot be negative!")
    if ( (nb_k>1) or (k[0]>0) ) and ( (nb_R>1) or (R[0]>0) ):
                          raise ValueError("specify either k or radius R, but not both!")
    if (nb_k>1) and (PNP.min(PNP.diff(k))<0):   raise ValueError("k should be sorted (with increasing values)")
    if (nb_R>1) and (PNP.min(PNP.diff(R))<0):   raise ValueError("R should be sorted (with increasing values)")
    if (k[0]==0) and (R[0]==0):                 raise ValueError("specify at least k or radius R!")
    if centered and (order_max>7):              raise ValueError("maximal order of central moments cannot exceed 7; use non-centered moments instead")

    if (kernel_type<0): kernel_type = commons.current_kernel_type
    if (scale<=0):      scale = commons.current_obs_scale;
    set_kernel(kernel_type, scale)                      # kernel should be set and parameterized before (with "set_kernel")
    
    nn_statistics.tree_k_max=nn_max                     # if (-1) then auto set to 1/10 of available points in x

    if verbosity: 
        print("computing moments of order 1 up to", order_max, "end=, ")
        PC.get_kernel()

    if (k[0]>0):
        if (k[nb_k-1]>=npts_in):    raise ValueError("imposed k is larger than the number of input points!")
        if verbosity: print("fixed k computation", end=" ")
        dists   = PNP.zeros((nb_k,npts_out), dtype=PNP.float64) 
        moments = PNP.zeros((order_max*nb_k*nA,npts_out), dtype=PNP.float64)
        if (nb_k==1):
            if verbosity: print("1 value of k :", k[0], end=" ")
#            if use_kernel:
#                if verbosity: print("using special kernel, with scale %2.2f", obs_scale)
#                nn_statistics.compute_stats_kernel_fixed_k_threads(&x[0,0], &A[0,0], npts_in, nx, nA, &y[0,0], npts_out, k[0], &moments[0,0], order_max, centered, &dists[0,0])
#            else:
#                if verbosity: print("no kernel")
            nn_statistics.compute_stats_fixed_k_threads(&x[0,0], &A[0,0], npts_in, nx, nA, &y[0,0], npts_out, k[0], &moments[0,0], order_max, centered, &dists[0,0])
            mom = PNP.asarray(moments).reshape(order_max, nA, npts_out)
        else:
            if verbosity: print("multiple values of k :", PNP.array(k))
            nn_statistics.compute_stats_multi_k_threads(&x[0,0], &A[0,0], npts_in, nx, nA, &y[0,0], npts_out, &k[0], nb_k, &moments[0,0], order_max, centered, &dists[0,0])
            mom = PNP.asarray(moments).reshape(order_max, nb_k, nA, npts_out)
        ret = [PNP.sqrt(PNP.asarray(dists))]            # we return the distances
        
    elif (R[0]>0):   
        if verbosity: print("fixed R computation", end=" ")
        R2      = PNP.zeros(nb_R, dtype=PNP.float64)
        for ratou in range(nb_R): R2[ratou]=R[ratou]**2 # internal code expects squared distances
        nnn     = PNP.zeros((nb_R,npts_out), dtype=PNP.intc)
        moments = PNP.zeros((order_max*nb_R*nA,npts_out), dtype=PNP.float64)
        if (nb_R==1):
            if verbosity: print("1 value of R :", R[0]) 
            nn_statistics.compute_stats_fixed_R_threads(&x[0,0], &A[0,0], npts_in, nx, nA, &y[0,0], npts_out, R2[0], &moments[0,0], order_max, centered, &nnn[0,0])
            mom = PNP.asarray(moments).reshape(order_max, nA, npts_out)
        else:
            if verbosity: print("multiple values of R :", PNP.array(R))
            nn_statistics.compute_stats_multi_R_threads(&x[0,0], &A[0,0], npts_in, nx, nA, &y[0,0], npts_out, &R2[0], nb_R, &moments[0,0], order_max, centered, &nnn[0,0])
            mom = PNP.asarray(moments).reshape(order_max, nb_R, nA, npts_out)
        ret = [PNP.asarray(nnn)]                        # we return the nb of neighbors
            
    for i in range(mom.shape[0]): ret.append(mom[i])    # we append all required moments
    return ret
        


def set_nn_max(int k=nn_statistics.k_default):
    """
    sets the maximal value of allowed number of nearest neighbors
    
    :param k: an integer (default=5)
    :returns: no output
    """
    nn_statistics.tree_k_max=k
    
def get_nn_max():
    """
    gets the current maximal value of allowed number of nearest neighbors
    
    :param none:
    :returns: an integer
    """
    return(nn_statistics.tree_k_max)


