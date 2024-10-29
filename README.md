# nn_stats: nearest-neighbors statistics
an efficient C/C++ library integrated with Python to estimate local averages (and their corresponding variance) using nearest neighbors estimates.
This library relies on the [ANN library](http://www.cs.umd.edu/~mount/ANN/) by David Mount and Sunil Arya.

# compilation and installation
- run ./configure and eventuallly solve the issues by installing missing programs and libraries (e.g.: "apt install libtool-bin" on Linux if asked to do so)
- then run "make python" to produce the library.
  
# how to use in Python
- "make python" will both compile the library and install it in your python path, which depends on you current environment. You should select your environment first, then run "./configure" and "make python", in order to have the library and its functions available in your favored environment.
- there are (or will be) examples in the bin/python subdirectory: please look at them to learn how to import and use the library, which should be as easy as:
<pre><code>
import numpy as np
import nn_stats as ns

Npts = 100000
locations = np.random.randn(2,Npts)
values    = np.random.randn(1,Npts)

Npts_new  = 100
loc_new   = np.random.randn(2, Npts_new)

k=np.array([5])
R=np.array([0.5])

mean, var, = ns.compute_local_stats(pos, val, y, k=k)   # imposed k
mean, var, = ns.compute_local_stats(pos, val, y, R=R)   # imposed R
</code></pre>

# important remarks

As the library is built for maximum efficiency (in both speed and memory usage), you should respect the following:

- all parameters for the function "compute_local_stats" are expected to be Numpy arrays. If you have a list, convert it first to a nd-array:
<pre><code>
k = [5, 10, 15]         # a Python list, not efficient and will throw an exception
mean, var, = ns.compute_local_stats(pos, val, y, k=k)   # throws an exception

k = np.array([5, 10, 15])     # a nd-array, as expected
mean, var, = ns.compute_local_stats(pos, val, y, k=k)   # works OK
</code></pre>

- parameters k and R are expected to be sorted, i.e., their values are increasing with the index: 
<pre><code>
k[i-1] <= k[i] # True for any valid index 1 <= i < size(k) 
</code></pre>

- parameters "positions" and "observables" are *2d-arrays*, with their .shape[0] being respectively the space dimension (i.e., the number of coordinates) in "positions" and the nb of observables. 
Their .shape[1] is simply the number of availabe points, which should be the same for "positions" and "observables".

# notes
this is still under develpment...

