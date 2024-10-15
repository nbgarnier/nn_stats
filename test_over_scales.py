import numpy
import entropy.entropy as entropy
import entropy.tools as tools

x       = numpy.random.randn(1,1000)
tau_set = numpy.linspace(1,100,50, dtype=int)
print(tau_set)



#%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# load Modane data (pre-converted)
Npts   = 400000  # nb de pts en temps à conserver
#Npts   = 25690112 # all points

data = numpy.load("/Users/ngarnier/Documents/research/entropy/2019-01-16-ir/data_Modane.npy")
print("initial signal size :", data.size)
x=tools.reorder(data[0,0:Npts])
print("        I keep only :", x.size)
#t=numpy.arange(signal.size, dtype=float)
print("x has shape", x.shape)

fs=25000;            # sampling frequency
tau_eta = 1.8e-4*fs; # Kolmogorov scale (cf EPL)
tau_L   = 3.6e-2*fs; # large scale (cf EPL)

#tau_set = numpy.power(1.3, numpy.arange(25)+2).astype(int)  # for Modane
tau_set = numpy.arange(1,6,1)
tau_set = numpy.append(tau_set, numpy.arange(6,12,2));
tau_set = numpy.append(tau_set, numpy.arange(15,50,5));
tau_set = numpy.append(tau_set, numpy.arange(50,100,10));
tau_set = numpy.append(tau_set, numpy.arange(100,500,50));
tau_set = numpy.append(tau_set, numpy.arange(500,1000,100));
tau_set = numpy.append(tau_set, numpy.arange(1000,5000,500));
tau_set = numpy.append(tau_set, numpy.arange(5000,20000,2000));
print("tau", tau_set)
data_name="Modane"


entropy.set_sampling(Theiler=4, N_eff=4000, N_real=10)
H = tools.compute_over_scales(entropy.compute_entropy_increments, tau_set, x, k=4, verbosity_timing=1, get_samplings=1)


print(entropy.get_last_info())
entropy.get_last_sampling(1)

print(H)
