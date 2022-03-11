# Attempt enforcing alternate coefficients to have opposite sign 
# plus constrain their magnitudes based on higher order coefficients
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants
import random as rndm
import math as math
from math import exp, expm1, sqrt
from math import factorial as mf
from scipy.optimize import fmin
from scipy import interpolate
from scipy.optimize import curve_fit
import emcee
import acor
import scipy.signal as signal

PI=scipy.constants.pi
HH=scipy.constants.h
KK=scipy.constants.k
HbK=HH/KK
HbK9=HbK*1.0e9
neg_inf = -np.inf

np.set_printoptions(precision=20)

# # Read in the recombination line template and plot 

f = open('general_data/HI.HeI.HeII.dat','r')
#f = open('general_data/HI.dat','r')
lines = f.readlines()
f.close()
x1 = []
y1 = []
for line in lines:
    p = line.split()
    x1.append(float(p[0]))
    y1.append(float(p[1]))
   
x1 = np.asfarray(x1)
y1 = np.asfarray(y1)
# Limit all values to below 4 THz
x1 = x1[x1<=4000]*1.0e9 
y1 = y1[x1<=4000]
# Create new frequency array for interpolation later - 5 MHz to 4 THz
xnew = np.arange(50e6, 1e12,10e6)

#mini = np.array(signal.argrelmin(np.log10(y1)))
#maxi = np.array(signal.argrelmax(np.log10(y1)))
mini = np.array(signal.argrelmin(y1))
maxi = np.array(signal.argrelmax(y1))
mini = mini[0,:]
maxi = maxi[0,:]
maxi = maxi[0:-1]
print np.shape(mini)
print np.shape(maxi)
s_contrast1 = y1[maxi]-y1[mini]
#plt.plot(s_contrast1)
#plt.show()
plt.plot(x1, y1)
plt.scatter(x1[mini],y1[mini],c='r')
plt.scatter(x1[maxi],y1[maxi],c='g')
plt.xlim([40,300])
#plt.ylim([np.min(y1),1e-27])
plt.xlabel('Frequency (GHz)')
plt.ylabel('Brightness')
plt.show()

tck_min = interpolate.splrep(x1[mini],y1[mini], s=0)
tck_max = interpolate.splrep(x1[maxi],y1[maxi], s=0)
y_min = interpolate.splev(xnew, tck_min, der=0)
y_max = interpolate.splev(xnew, tck_max, der=0)
plt.plot(x1,y1,ls=':')
plt.scatter(x1[mini],y1[mini],c='r')
plt.scatter(x1[maxi],y1[maxi],c='r')
plt.plot(xnew,y_min,c='g',ls='--')
plt.plot(xnew,y_max,c='g',ls='--')
#plt.xlim(50,1000)
#plt.ylim(-0.5e-26,0.5e-26)
plt.show()
y_diff = y_max[xnew<1000] - y_min[xnew<1000]
plt.plot(np.log10(xnew[xnew<1000]),np.log10(y_diff))
plt.show()














# # b = np.r_[True, a[1:] < a[:-1]] & np.r_[a[:-1] < a[1:], True]
# # peakind = signal.find_peaks_cwt(y1, np.arange(1,15))
# # peakind, x1[peakind], y1[peakind]
# # print  'Plot the recombination line data'
# plt.plot(np.log10(x1), np.log10(y1))
# # plt.scatter(x1[peakind],y1[peakind],c='r')
# plt.xlim([-1,1])
# plt.ylim([-30,-25])
# # plt.ylim([np.min(y1),0.8e-27])
# # plt.xlabel('Frequency (GHz)')
# # plt.ylabel('Brightness')
# plt.show()
# # signal_contrast = np.diff(y1[peakind])
# # sig_cont_log = np.log10(signal_contrast)

# # plt.plot(sig_cont_log)
# # # plt.xlim([1,10])
# # # plt.ylim([np.min(signal_contrast),1.0e-27])
# # plt.show()
