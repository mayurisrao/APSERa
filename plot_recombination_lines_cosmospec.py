# plot_recomb_intensity_cosmospec.py
# Plot additive recombination spectrum
# Author: Mayuri S.Rao 
# Date: 11 March 2022

import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
import scipy.io
from scipy import interpolate

#path     = '/home/mayuris/workspace/RRI_DISTORTION_Lab/APSERa/Theory/'
freq, data     = np.loadtxt('total_spec_new.dat')
#freq, data     = data['Data2']/1e3  # Kelvin units
#fr       = np.loadtxt('freq_saras.txt') #Frequency in MHz
#signal   = interpolate.interp1d(fr, data)

#freq_custom = np.linspace(1, 10, 1001) #Use the GHz frequencies you want for the signals
#T21_full = signal(freq_custom) #This generates the atlas (nsignal, nfrequency)

plt.figure(figsize=(12,6))
plt.plot(freq, data)
#plt.plot(freq_custom, vanilla_model, lw=3.0, color='k', label='vanilla model')
plt.grid()
plt.legend()
plt.xlabel("Frequency (GHz)")
plt.ylabel(r"$\delta I_v [Jm^{-2}s^{-1}Hz^{-1}sr^{-1}]$")
plt.show()