# nn_stats: nearest-neighbors statistics
an efficient C/C++ library integrated with Python to estimate local averages (and their corresponding variance) using nearest neighbors estimates.
This library relies on the [ANN library](http://www.cs.umd.edu/~mount/ANN/) by David Mount and Sunil Arya.

# compilation and installation
- run ./configure and eventuallly solve the issues by installing missing programs and libraries (e.g.: "apt install libtool-bin" on Linux if asked to do so)
- then run "make python" to produce the library.
  
# how to use in Python
- "make python" will both compile the library and install it in your python path, which depends on you current environment. You should select your environment first, then run "./configure" and "make python", in order to have the library and its functions available in your favored environment.
- there are examples in the bin/python subdirectory: please look at them to learn how to import and use the library, which should be as easy as:
<pre><code>
import numpy as np
import nn_stats as ns

Npts = 100000
locations = np.random.randn(2,Npts)
values    = np.random.randn(1,Npts)

Npts_new  = 100
loc_new   = np.random.randn(2, Npts_new)

k=5
R=0.5

mean, var, = ns.compute_local_stats(pos, val, y, k=k)   # imposed k
mean, var, = ns.compute_local_stats(pos, val, y, R=R)   # imposed R
</code></pre>

# notes
this is still untder develpment...

