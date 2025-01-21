# -*- coding: ascii -*-
import numpy as np
import nn_stats as ns 

from time import time
from numpy.random import default_rng
np.random.seed(1234)
rng = default_rng()

import matplotlib
import matplotlib.pyplot as plt

#from matplotlib.path import Path
plt.rc('text', usetex=True)
font = {'family' : 'serif',  # was 'normal'
        'sans-serif':'Helvetica',
        'serif'  : 'cmr10',    # LaTeX default
        'weight' : 'normal',   # was 'bold'
        'size'   : 16}
matplotlib.rc('font', **font)
matplotlib.rcParams['axes.formatter.use_mathtext'] = True


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# set parameters
Npts    = 500 # 100000
k       = Npts//2       # for fixed k subplot
R       = 0.2         # for fixed R subplot
d       = 0.5
sigma_x = 1

verbosity=0

# here, we can play with nb of cores:
#ns.multithreading(1)       # single core
ns.multithreading("auto")   # max cores
ns.multithreading("info")   # prints informations
ns.set_verbosity(verbosity)         # set default verbosity to 1, to have some little messages
ns.get_verbosity()

radius  = 5.   # for inital points and output grid

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# prepare initial set of points, and output grid (they are different in this example)
def Disk(radius, Npts=100):
    theta=2*np.pi*np.random.uniform(0,1,Npts); #angular coordinates 
    rho  =radius*np.sqrt(np.random.uniform(0,1,Npts)); #radial coordinates 
     
    #Convert from polar to Cartesian coordinates
    xx = rho * np.cos(theta)
    yy = rho * np.sin(theta)
    return xx, yy

# input set, with values for the observables:
pos=np.array(Disk(radius, Npts=Npts), dtype=float)
val=np.ones((1,pos.shape[1]), dtype=float)
print("input positions :", pos.shape, "with observables :", val.shape)

# output grid (will be of size Nx_out^2):
Nx_out = 50
a = np.arange(1,Nx_out+1)*2-(Nx_out+1)
b = np.ones(Nx_out+1)
y = np.array((np.outer(a,b).flatten(),np.outer(b,a).flatten()), dtype=float)
y = y/Nx_out*radius
print("output positions :", y.shape)
print()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# prepare plot
Fig = plt.figure(figsize=(15,10))
P1  = Fig.add_subplot(2,3,1); P1.set_title("initial locations (%d pts)" %Npts); P1.set_xlim(-radius, radius); P1.set_ylim(-radius, radius)
P2  = Fig.add_subplot(2,3,2); P2.set_title("R=%1.1f, no kernel" %R);      P2.set_xlim(-radius, radius); P2.set_ylim(-radius, radius)
P3  = Fig.add_subplot(2,3,3); P3.set_title("k=%d, no kernel" %k);         P3.set_xlim(-radius, radius); P3.set_ylim(-radius, radius)
P4  = Fig.add_subplot(2,3,4); P4.set_title("R=%1.1f, Gaussian kernel (d=%1.1f)" %(R, d));P4.set_xlim(-radius, radius); P4.set_ylim(-radius, radius)
P5  = Fig.add_subplot(2,3,5); P5.set_title("R=%1.1f, regular kernel" %R); P5.set_xlim(-radius, radius); P5.set_ylim(-radius, radius)
P6  = Fig.add_subplot(2,3,6); P6.set_title("k=%d, regular kernel" %k);    P6.set_xlim(-radius, radius); P6.set_ylim(-radius, radius)

P1.scatter(pos[0], pos[1], marker='.', edgecolor='b', facecolor='none', alpha=0.5 )

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# computations, and timing measurements

# single value of k, no kernel:
print("fixed k =", k, "no kernel       ", end=" ")
t1=time()
[dists, A_mean] = ns.compute_local_stats(pos, val, y, k=np.array([k], dtype=np.intc), order_max=1)
print("\telapsed time", time()-t1)
P3.scatter(y[0,:].flatten(), y[1,:].flatten(), c=A_mean/dists, marker='o')

# single value of k, regular kernel:
ns.set_kernel(1, 1); ns.get_kernel()
print("fixed k =", k, "regular kernel  ", end=" ")
t1=time()
[dists, A_mean, A_var] = ns.compute_local_stats(pos, val, y, k=np.array([k], dtype=np.intc))
print("\telapsed time", time()-t1)
P6.scatter(y[0,:].flatten(), y[1,:].flatten(), c=A_mean/dists, marker='o')

# single value of k, Gaussian kernel:
# ns.set_kernel(2, d); ns.get_kernel()
# print("fixed k =", k, "Gaussian kernel ", end=" ")
# t1=time()
# [dists, A_mean, A_var] = ns.compute_local_stats(pos, val, y, k=np.array([k], dtype=np.intc))
# print("\telapsed time", time()-t1)
# P4.scatter(y[0,:].flatten(), y[1,:].flatten(), c=A_mean/dists, marker='o')

# single value of R, no kernel:
ns.set_kernel(); ns.get_kernel()
print("fixed R =", R, "  no kernel       ", end=" ")
t1=time()
[nn, A_mean, A_var] = ns.compute_local_stats(pos, val, y, R=np.array([R]))
print("\telapsed time", time()-t1)
P2.scatter(y[0,:].flatten(), y[1,:].flatten(), c=A_mean*nn, marker='o')

# single value of R, regular kernel:
ns.set_kernel(1, 1); ns.get_kernel()
print("fixed R =", R, "  regular kernel  ", end=" ")
t1=time()
[n, A_mean, A_var] = ns.compute_local_stats(pos, val, y, R=np.array([R]))
print("\telapsed time", time()-t1)
P5.scatter(y[0,:].flatten(), y[1,:].flatten(), c=A_mean*nn, marker='o')

# single value of R, Gaussian kernel:
ns.set_kernel(2, d); ns.get_kernel()
print("fixed R =", R, "  Gaussian kernel ", end=" ")
t1=time()
[nn, A_mean, A_var] = ns.compute_local_stats(pos, val, y, R=np.array([R]))
print("\telapsed time", time()-t1)
P4.scatter(y[0,:].flatten(), y[1,:].flatten(), c=A_mean*nn, marker='o')

Fig.savefig("test_kernel.pdf")

# %%
