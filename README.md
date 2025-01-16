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
  * either a set of values of k (parameter "k") or a set of values of R (parameter "R")

The following parameters can also be provided:
  * a set of observables values taken on the initial locations (parameter "A"). By default, this set is empty and no moments are computed.
  * a set of "destination" locations (parameter "y") where the statistics of observables will be computed. By default, he "initial" locations are used.
  * a maximal order for moments computation (parameter "order_max", new in v0.7.2). Moments of order 1, 2, ..., order_max will be computed. By default, order_max=2. Note that if you provide order_max=0, no moments will be computed, even if some observables (parameter "A") are provided.
  * a boolean (parameter "centered") to indicate central moments (True) or natural moments (False) are requested. By default, observables are not centered, and natural moments are returned. (new in v0.7.3)
  * a maximal number of neighbors to consider when performing a fixed-radius search (parameter "nn_max") (this may slow down the library, do it at your own risk!). By default, the function search for at most x.shape[1]/10 neighbors, i.e., 10% of the available data. You can specify to search for more points with this parameter. (new in v0.8.0)
  
The function returns (1 + order_max) Numpy arrays, in the following order:
  * first, the radii (if k was provided) or the k (if R was provided).
  * then the moments of observable(s) "A" computed locally around the 'destination' locations provided with parameter "y", starting with the first moment (expected value), and up to the moment or order "order_max".

# additional documentation

There are detailed examples in the examples/ subdirectory: please look at them to learn how to import and use the library, which should be as easy as:
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

R1, mean, var = ns.compute_local_stats(locations, values, loc_new, k=k)   # imposed k (indeed 3 values of k)
k1, mean, var = ns.compute_local_stats(locations, values, loc_new, R=R)   # imposed R (1 value of R)
</code></pre>

# important remarks

As the library is built for maximum efficiency (in both speed and memory usage), you should respect the following:

- all parameters for the function "compute_local_stats" are expected to be Numpy arrays. If you have a list, convert it first to a nd-array:
<pre><code>
k = [5, 10, 15]                             # a Python list, not efficient and will throw an exception
R, mean, var = ns.compute_local_stats(locations, values, loc_new, k=k)   # throws an exception

k = np.array([5, 10, 15], dtype=np.intc)    # a nd-array, as expected
R, mean, var = ns.compute_local_stats(locations, values, loc_new, k=k)   # works OK
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

if you are just interested in the number of neighbors (given a fixed radius R) or in the radius where the k-th neighbors lies, you can invoke the function "compute_local_stats" in the following two different ways:
 
* without providing any observable. This is for example done with:
<pre><code> 
R, _, _ = ns.compute_local_stats(locations, y=loc_new, k=k)   # imposed k -> returns R at new locations loc_new
k, _, _ = ns.compute_local_stats(locations, y=loc_new, R=R)   # imposed R -> returns k at new locations loc_new
</code></pre>
Note that in that case, if you want to provide a set of "destination" locations, you have to explicitly prefix them with "y=" as in the example code above. If you do not do so, the function will expect its second parameter to be the observables, while you provided "destination" locations.

* imposing order_max=0. This is for example done with:
<pre><code> 
R = ns.compute_local_stats(locations, A=values, y=loc_new, k=k, order_max=0)   # imposed k -> returns R at new locations loc_new
k = ns.compute_local_stats(locations, A=values, y=loc_new, R=R, order_max=0)   # imposed R -> returns k at new locations loc_new
</code></pre>

Note using the first way (and having a non-zero order_max) there will be (1+order_max) output variables (so 1+2=3 by default), which may be desirable depending on the style of your script, although returned moments variables should be empty. 

# notes
This is still under development, but has been tested OK in most common situations. Open an issue if some trouble arises that puzzles you.

