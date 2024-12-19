# -*- coding: ascii -*-
import numpy as np
import nn_stats as ns 

from time import time
from numpy.random import default_rng
np.random.seed(1234)
rng = default_rng()

import matplotlib
import matplotlib.pyplot as plt
#import pylab # for color manipulation

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
k=30
R=0.5

Npts    = 2500
ndim    = 1
sigma_x = 1

radius  = 5.   # for inital points and output grid

# here, we can play with nb of cores:
#ns.multithreading(1)       # single core
ns.multithreading("auto")  # max cores
ns.multithreading("info")  # prints informations
ns.set_verbosity(1)  

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# prepare initial set of points, and output grid (they are different in this example)
def Disk(radius, Npts=100):
    theta=2*np.pi*np.random.uniform(0,1,Npts); #angular coordinates 
    rho  =radius*np.sqrt(np.random.uniform(0,1,Npts)); #radial coordinates 
     
    #Convert from polar to Cartesian coordinates
    xx = rho * np.cos(theta)
    yy = rho * np.sin(theta)
    return xx, yy


pos=np.array(Disk(radius, Npts=Npts), dtype=float)
val=np.ones((1,pos.shape[1]), dtype=float)
print("input positions :", pos.shape, "with observables :", val.shape)

# output grid (will be of size Nx_out^2):
Nx_out = 45
a = np.arange(1,Nx_out+1)*2-(Nx_out+1)
b = np.ones(Nx_out+1)
y = np.array((np.outer(a,b).flatten(),np.outer(b,a).flatten()), dtype=float)
y = y/Nx_out*radius
print("output positions :", y.shape)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# prepare plot
Fig = plt.figure(figsize=(15,10))
P1  = Fig.add_subplot(2,3,1); P1.set_title("initial locations");    P1.set_xlim(-radius, radius); P1.set_ylim(-radius, radius)
P2  = Fig.add_subplot(2,3,2); P2.set_title("fixed R=%1.1f" %R);     P2.set_xlim(-radius, radius); P2.set_ylim(-radius, radius)
P3  = Fig.add_subplot(2,3,3); P3.set_title("fixed R=%1.1f" %(2*R)); P3.set_xlim(-radius, radius); P3.set_ylim(-radius, radius)
P4  = Fig.add_subplot(2,3,4); P4.set_title("fixed k=%d" %k);        P4.set_xlim(-radius, radius); P4.set_ylim(-radius, radius)
P5  = Fig.add_subplot(2,3,5); P5.set_title("fixed R=%1.1f" %(3*R)); P5.set_xlim(-radius, radius); P5.set_ylim(-radius, radius)
P6  = Fig.add_subplot(2,3,6); P6.set_title("fixed R=%1.1f" %(4*R)); P6.set_xlim(-radius, radius); P6.set_ylim(-radius, radius)

P1.scatter(pos[0], pos[1], marker='.', edgecolor='b', facecolor='none', alpha=0.5 )

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# computations, and timing measurements

# single value of k (for comparison):
print("\nfixed k =", k, end=" ")
t1=time()
[A_mean, A_var, dists] = ns.compute_local_stats(pos, val, y, k=np.array([k], dtype=np.intc))
print("\telapsed time", time()-t1)
P4.scatter(y[0,:].flatten(), y[1,:].flatten(), c=A_mean/dists, marker='o')

# multiple values of R:
print("\nfixed R =", [R, 2*R, 3*R, 4*R], end=" ")
t1=time()
[A_mean, A_var, nn] = ns.compute_local_stats(pos, val, y, R=np.array([R, 2*R, 3*R, 4*R]))
print("\telapsed time", time()-t1)
P2.scatter(y[0,:].flatten(), y[1,:].flatten(), c=A_mean[0,0,:]*nn[0,:], marker='o')
P3.scatter(y[0,:].flatten(), y[1,:].flatten(), c=A_mean[1,0,:]*nn[1,:], marker='o')
P5.scatter(y[0,:].flatten(), y[1,:].flatten(), c=A_mean[2]*nn[2], marker='o')
P6.scatter(y[0,:].flatten(), y[1,:].flatten(), c=A_mean[3]*nn[3], marker='o')
print("output values of size", A_mean.shape, A_var.shape, "and", nn.shape)

Fig.savefig("test_fixed_R.pdf")

