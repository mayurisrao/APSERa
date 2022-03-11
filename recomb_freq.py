# Python code to read in the recombination line spectrum and plot along with
# markers at the expected line frequencies
#
# v0 19 Oct 2013

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants
from scipy.interpolate import interp1d

RYD = scipy.constants.Rydberg
cc = scipy.constants.c
print "Using Rydberg constant * velocity of light = ",RYD*cc," Hz "
zz = float(raw_input("Enter redshift (1+z) for the recombination lines: "))

f = open('HI.HeI.HeII.dat', 'r')  # read in the recombination line spectrum
npoints=0
Recl_F = []
Recl_I = []
for line1 in f:
	line1 = line1.strip()
	line1 = line1.split()
	line1 = np.asfarray(line1)
	npoints += 1
	Recl_F.append(line1[0])
	Recl_I.append(line1[1])
Recl_F = np.asfarray(Recl_F)
Recl_I = np.asfarray(Recl_I)
istart=400 
iend=4676
Recl_F = Recl_F[istart:iend]
Recl_I = Recl_I[istart:iend]
npoints=len(Recl_F)

logRecl_F = np.log10(Recl_F)
logRecl_I = np.log10(Recl_I)

# Create an array with the line frequencies (in log 10 Hz units)
linefreq=[]
for i in range (60,1,-1):
	t1=1.0/float(i*i)
	t2=1.0/float((i-1)*(i-1))
	linefreq.append(np.log10(RYD*(t2-t1)*cc*1.0e-9/zz))

# Create markers to identify the line locations on the plots
# Take nearest neighbour in frequency axis and make that intensity the y-axis value for the marker 
xx=[]
xx_linear=[]
yy=[]
for i in range (len(linefreq)):
	dist=3.0e10
	jmin=0
	for j in range (len(logRecl_I)):
		cdist = np.absolute(linefreq[i] - logRecl_F[j])
		if ( cdist < dist ): 
			jmin = j
			dist = cdist
	xx.append(logRecl_F[jmin])
	xx_linear.append(Recl_F[jmin])
	yy.append(logRecl_I[jmin])	
xx = np.asfarray(xx)
xx_linear = np.asfarray(xx_linear)
yy = np.asfarray(yy)

# define function for cubic spline interpolation between the markers.
# this spline fit is a baseline that may be subtracted to see the lines
yy2 = interp1d(xx, yy, kind='cubic') 

# define a x-axis array over which the interpolation will be done
index_lower=np.where(logRecl_F==np.min(xx))
index_upper=np.where(logRecl_F==np.max(xx))
logRecl_F_limited = logRecl_F[index_lower[0][0]:index_upper[0][0]]

# plot the spectrum plus markers plus interpolated baseline
ax1 = plt.subplot(1,1,1)
ax1.plot(logRecl_F,logRecl_I,linestyle='-')
ax1.scatter(xx,yy,color='r')
ax1.plot(logRecl_F_limited,yy2(logRecl_F_limited),linestyle='--')
plt.xlabel('log10 Frequency (GHz)')
plt.ylabel('log10 Intensity')
plt.show()

# Plot the residue
Recl_F_limited = Recl_F[index_lower[0][0]:index_upper[0][0]]
logRecl_I_limited = logRecl_I[index_lower[0][0]:index_upper[0][0]]
logRecl_I_residue = logRecl_I_limited - yy2(logRecl_F_limited)
logyy_residue = yy - yy2(xx)
yy_residue = 10.0**yy - 10.0**yy2(xx)
Recl_I_residue = 10.0**logRecl_I_limited - 10.0**yy2(logRecl_F_limited)
ax1 = plt.subplot(1,1,1)
ax1.plot(Recl_F_limited,Recl_I_residue,linestyle='-')
ax1.scatter(xx_linear,yy_residue,color='r')
plt.xlabel('Frequency (GHz)')
plt.ylabel('Intensity')
plt.show()

f.close()




