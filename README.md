# nn_stats: nearest-neighbors statistics
an efficient C/C++ library integrated with Python to estimate local averages (and their corresponding variance) using nearest neighbors estimates.
This library relies on the [ANN library](http://www.cs.umd.edu/~mount/ANN/) by David Mount and Sunil Arya.

Please cite any use of this library with the following DOI:  
[![DOI](https://zenodo.org/badge/873066948.svg)](https://doi.org/10.5281/zenodo.14523934)



# compilation and installation
- run ./configure and eventually solve the issues by installing missing programs and libraries (e.g.: "apt install libtool-bin" on Linux if asked to do so)
- then run "make python" to produce the library. This will both compile the library and install it in your python path, which depends on you current environment. You should select your environment first, then run "./configure" and "make python", in order to have the library and its functions available in your favored environment.
  
# how to use in Python
There is a single function, called **"compute_local_stats"**, which can be invoked in 2 different ways:
  * by imposing a set of values of k (numbers of neighbors to consider)
  * by imposing a set of values of R (radii to consider)

 The function expects (these are mandatory):
  * a set of initial locations (parameter "x") in a n-dimensional space
  * either a set of alues of k (parameterr "k") or a set of values of R (parameter "R")

The following parameters can also be provided:
  * a set of observables values taken on the initial locations (parameter "A"); by default, this set is empty and no moments are computed.
  * a set of "destination" locations (parameter "y") where the statistics of observables will be computed; by default, if you do not provide these "destination" locations, the "initial" locations will be used.
  
There are examples in the examples/ subdirectory: please look at them to learn how to import and use the library, which should be as easy as:
<pre><code>
import numpy as np
import nn_stats as ns

Npts = 100000
locations = np.random.randn(2,Npts)
values    = np.random.randn(1,Npts)

Npts_new  = 100
loc_new   = np.random.randn(2, Npts_new)

k=np.array([5, 10, 15], dtype=np.intc)    # multiple values of k
R=np.array([0.5])                         # single value of R

mean, var, R1 = ns.compute_local_stats(locations, values, loc_new, k=k)   # imposed k (indeed 3 values of k)
mean, var, k1 = ns.compute_local_stats(locations, values, loc_new, R=R)   # imposed R (1 value of R)
</code></pre>

# important remarks

As the library is built for maximum efficiency (in both speed and memory usage), you should respect the following:

- all parameters for the function "compute_local_stats" are expected to be Numpy arrays. If you have a list, convert it first to a nd-array:
<pre><code>
k = [5, 10, 15]                             # a Python list, not efficient and will throw an exception
mean, var, R = ns.compute_local_stats(locations, values, loc_new, k=k)   # throws an exception

k = np.array([5, 10, 15], dtype=np.intc)    # a nd-array, as expected
mean, var, R = ns.compute_local_stats(locations, values, loc_new, k=k)   # works OK
</code></pre>

- parameter k (if used) is expected of type "intc". You can set it to be this way like this:
<pre><code>
k = [5, 10, 15]                 # a Python list, easy to read
k = np.array(k, dtype=np.intc)  # a NumPy array, of type "intc"
</code></pre>

- parameters k and R are expected to be sorted, i.e., their values must be increasing with the index: 
<pre><code>
k[i-1] <= k[i] # True for any valid index 1 <= i < size(k) 
</code></pre>

- parameters "locations" and "observables" are *2d-arrays*, with their .shape[0] being respectively the space dimension (i.e., the number of coordinates) in "positions" and the number of observables. 
Their .shape[1] is simply the number of available points, which should be the same for "positions" and "observables".

- if there is just 1 observable, i.e., if values.shape[0] equals 1 (as in the example above), then the returned values "mean" and "var" have shape (k.size, loc_new.shape[1]) while the third returned value (R) has shape (1, loc_new.shape[1]).

- if there are more than just 1 observable, i.e., if values.shape[0] is larger than 1 in the above example, then the returned values "mean" and "var" have shape (k.size, values.shape[0], loc_new.shape[1]) while the third returned value (R) has shape (1, loc_new.shape[1]) (i.e., the same shape as if just 1 observable was provided).


# other remarks

 if you are just interested in the number of neighbors (given a fixed radius R) or in the radius where the k-th neighbors lies, you can invoke the function "compute_local_stats" without providing any observable. This is for example done with:
<pre><code> 
mean, var, R = ns.compute_local_stats(locations, y=loc_new, k=k)   # imposed k -> returns R at new locations loc_new
mean, var, k = ns.compute_local_stats(locations, y=loc_new, R=R)   # imposed R -> returns k at new locations loc_new
</code></pre>
note that there are still 3 output variables, even in that case (mean and var should be empty).

# notes
this is still under development, but has been tested OK in most common situations (with 1 observable only though).

