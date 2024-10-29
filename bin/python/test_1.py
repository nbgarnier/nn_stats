# -*- coding: ascii -*-
import numpy as np
import nn_stats as ns # python

from time import time
from numpy.random import default_rng
np.random.seed(1234)
rng = default_rng()

pi = np.pi
k=150
R=1.8

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

#%load_ext Cython
#%load_ext wurlitzer

#%%
def print_result(message, value_from_function):
#    time1=time()
#    value_from_function=function
#    time2=time()-time1
    print(message, "\t%2.5f" %value_from_function, end="")
    tmp=entropy.get_last_info()
    print(" +/- %2.5f (%d (%d) eff. pts, %d errors)" %(tmp[0], tmp[3], tmp[4], tmp[2]), end="")
    
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# test with fixed k or fixed radius
###############################################################################
Npts    = 1000
ndim    = 1
sigma_x = 1

n_embed = 2       # embedding
stride  = 5      # stride (tau)

ns.multithreading(1)

def Disk(radius, Npts=100, center=(0,0)):
    #Simulate Poisson point process
 #   Npts = np.random.poisson(lambda0*areaTotal);#Poisson number of points
    theta=2*np.pi*np.random.uniform(0,1,Npts); #angular coordinates 
    rho  =radius*np.sqrt(np.random.uniform(0,1,Npts)); #radial coordinates 
     
    #Convert from polar to Cartesian coordinates
    xx = rho * np.cos(theta)
    yy = rho * np.sin(theta)
     
    #Shift centre of disk to (xx0,yy0) 
    return xx+center[0], yy+center[1]
#%%

pos=np.array(Disk(5, Npts=Npts), dtype=float)
val=np.ones((1,pos.shape[1]), dtype=float)
print("input positions :", pos.shape, "with observables :", val.shape)

# output grid:
Nx_out = 5
a = np.arange(1,Nx_out+1)*2-(Nx_out+1)
b = np.ones(Nx_out+1)
y = np.array((np.outer(a,b).flatten(),np.outer(b,a).flatten()), dtype=float)
print("output positions :", y.shape)
#print(y)

# prepare plot
Fig = plt.figure(figsize=(15,10))
P1  = Fig.add_subplot(2,3,1); P1.set_title("initial locations")
P2  = Fig.add_subplot(2,3,2); P2.set_title("fixed k=%d" %k)
P3  = Fig.add_subplot(2,3,3); P3.set_title("fixed k=%d" %(2*k))
P4  = Fig.add_subplot(2,3,4); P4.set_title("fixed R=%1.1f" %R)
P5  = Fig.add_subplot(2,3,5); P5.set_title("fixed k=%d" %k)
P6  = Fig.add_subplot(2,3,6); P6.set_title("fixed k=%d" %(2*k))

P1.scatter(pos[0], pos[1], marker='.', edgecolor='b', facecolor='none', alpha=0.5 )

# single value of k:
print("\nfixed k =", k, end=" ")
[A_mean, A_var, dists] = ns.compute_local_stats(pos, val, y, k=np.array([k], dtype=np.intc))
P2.scatter(y[0,:].flatten(), y[1,:].flatten(), c=A_mean/dists, marker='o')
print(dists)
print(A_mean)
print("\nfixed k =", k, end=" ")
[A_mean, A_var, dists] = ns.compute_local_stats(pos, val, y, k=np.array([2*k], dtype=np.intc))
P3.scatter(y[0,:].flatten(), y[1,:].flatten(), c=A_mean/dists, marker='o')
print("output values of size", A_mean.shape, A_var.shape, "and", dists.shape)
print(dists)
print(A_mean)

# multiple values of k:
print("\nfixed k =", [k, 2*k], end=" ")
[A_mean, A_var, dists] = ns.compute_local_stats(pos, val, y, k=np.array([k, 2*k], dtype=np.intc))
#P4.scatter(y[0,:].flatten(), y[1,:].flatten(), c=A_mean/dists[0], marker='o')
#P5.scatter(y[0,:].flatten(), y[1,:].flatten(), c=A_mean/dists[1], marker='o')
print("output values of size", A_mean.shape, A_var.shape, "and", dists.shape)
print(dists)
print(A_mean)


print("\nfixed R =", R, end=" ")
A_mean, A_var, nnn = ns.compute_local_stats(pos, val, y, R=np.array([R]))
print("output values of size", A_mean.shape, A_var.shape, "and", nnn.shape)

P4.scatter(y[0,:].flatten(), y[1,:].flatten(), c=A_mean*nnn, marker='o')
print(nnn)
Fig.savefig("essai.pdf")
