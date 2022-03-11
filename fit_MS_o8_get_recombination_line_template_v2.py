# Fit to recombination lines template plus foreground model plus CMB
#
# Fit smooth function in log-log space
#
# Compute derivatives and if any zero crossings are deemed present 
# in the range then blow up the chisq if so.
#
# Finally calculate a Markov chain for the parameters of the fit


# *****IMPORTANT*****
# The values of X1 and X2 need to be changed to be the end limits 
# of the frequency range in GHz.
# *****IMPORTANT*****

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants
import random as rndm
import math as math
from math import exp, expm1, sqrt
from math import factorial as mf
from scipy.optimize import fmin
from scipy import interpolate
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

PI=scipy.constants.pi
HH=scipy.constants.h
KK=scipy.constants.k
HbK=HH/KK
HbK9=HbK*1.0e9
neg_inf = -np.inf

X1 = 1.0
X2 = 7.0
print ' '
print 'CAUTION: this version is hard-coded for freq range: ',X1,' to ',X2,' GHz'

X1LOG = np.log10(X1)
X2LOG = np.log10(X2)
XDIF = X1LOG-X2LOG

np.set_printoptions(precision=20)

func1 = lambda p, x, yt: (HbK9*x)/( np.exp( HbK9*x/(10.0**p[0]) )-1.0 ) + yt + \
					10.0**( \
					(p[1]) + \
					(p[2])*(np.log10(x)-X1LOG)**1 \
						  )

func11 = lambda p, x: (HbK9*x)/( np.exp( HbK9*x/(10.0**p[0]) )-1.0 ) + \
					10.0**( \
					(p[1]) + \
					(p[2])*(np.log10(x)-X1LOG)**1 \
						  )

func2 = lambda p, x, yt: (HbK9*x)/( np.exp( HbK9*x/(10.0**p[0]) )-1.0 ) + yt + \
					10.0**( \
					(p[1]) + \
					(p[2])*(np.log10(x)-X1LOG)**1 + \
					(p[3])*(np.log10(x)-X1LOG)**2  \
						  )

func22 = lambda p, x: (HbK9*x)/( np.exp( HbK9*x/(10.0**p[0]) )-1.0 ) +  \
					10.0**( \
					(p[1]) + \
					(p[2])*(np.log10(x)-X1LOG)**1 + \
					(p[3])*(np.log10(x)-X1LOG)**2  \
						  )

func3 = lambda p, x, yt: (HbK9*x)/( np.exp( HbK9*x/(10.0**p[0]) )-1.0 ) + yt + \
					10.0**( \
					(p[1]) + \
					(p[2])*(np.log10(x)-X1LOG)**1 + \
					(p[3])*(np.log10(x)-X1LOG)**2 + \
					(p[4])*(np.log10(x)-X1LOG)**3  \
						  )

func33 = lambda p, x: (HbK9*x)/( np.exp( HbK9*x/(10.0**p[0]) )-1.0 ) + \
					10.0**( \
					(p[1]) + \
					(p[2])*(np.log10(x)-X1LOG)**1 + \
					(p[3])*(np.log10(x)-X1LOG)**2 + \
					(p[4])*(np.log10(x)-X1LOG)**3  \
						  )

func4 = lambda p, x, yt: (HbK9*x)/( np.exp( HbK9*x/(10.0**p[0]) )-1.0 ) + yt + \
					10.0**( \
					(p[1]) + \
					(p[2])*(np.log10(x)-X1LOG)**1 + \
					(p[3])*(np.log10(x)-X1LOG)**2 + \
					(p[4])*(np.log10(x)-X1LOG)**3 + \
					(p[5])*(np.log10(x)-X1LOG)**4  \
						  )

func44 = lambda p, x: (HbK9*x)/( np.exp( HbK9*x/(10.0**p[0]) )-1.0 ) + \
					10.0**( \
					(p[1]) + \
					(p[2])*(np.log10(x)-X1LOG)**1 + \
					(p[3])*(np.log10(x)-X1LOG)**2 + \
					(p[4])*(np.log10(x)-X1LOG)**3 + \
					(p[5])*(np.log10(x)-X1LOG)**4  \
						  )


# define the function to be minimized by scipy.optimize.fmin
# chisq1 = lambda p, x, y, yt: sqrt(((func1(p,x,yt)-y)**2).sum()/float(len(x)))

# chisq2 = lambda p, x, y, yt: sqrt(((func2(p,x,yt)-y)**2).sum()/float(len(x)))

# def chisq3 (p, x, y, yt):
# 	k2 = (mf(2)/mf(0))*p[3] + (mf(3)/mf(1))*p[4]*(np.log10(x)-X1LOG)
# 	for i in range (len(x)-1):
# 		if k2[i+1]*k2[i] < 0.0: 
# 			return 100.0+min(i,len(x)-i)
# 	return sqrt(((func3(p,x,yt)-y)**2).sum()/float(len(x)))

# def chisq4 (p, x, y ,yt):
# 	k2 = (mf(2)/mf(0))*p[3] + (mf(3)/mf(1))*p[4]*(np.log10(x)-X1LOG) + \
# 		(mf(4)/mf(2))*p[5]*(np.log10(x)-X1LOG)**2
# 	k3 = (mf(3)/mf(0))*p[4] + (mf(4)/mf(1))*p[5]*(np.log10(x)-X1LOG)
# 	for i in range (len(x)-1):
# 		if k2[i+1]*k2[i] < 0.0 or k3[i+1]*k3[i] < 0.0: 
# 			return 100.0+min(i,len(x)-i)
# 	return sqrt(((func4(p,x,yt)-y)**2).sum()/float(len(x)))
	
chisq11 = lambda p, x, y: sqrt(((func11(p,x)-y)**2).sum()/float(len(x)))

chisq22 = lambda p, x, y: sqrt(((func22(p,x)-y)**2).sum()/float(len(x)))

def chisq33 (p, x, y):
	k2 = (mf(2)/mf(0))*p[3] + (mf(3)/mf(1))*p[4]*(np.log10(x)-X1LOG)
	for i in range (len(x)-1):
		if k2[i+1]*k2[i] < 0.0: 
			return 100.0+min(i,len(x)-i)
	return sqrt(((func33(p,x)-y)**2).sum()/float(len(x)))

def chisq44 (p, x, y):
	k2 = (mf(2)/mf(0))*p[3] + (mf(3)/mf(1))*p[4]*(np.log10(x)-X1LOG) + \
		(mf(4)/mf(2))*p[5]*(np.log10(x)-X1LOG)**2
	k3 = (mf(3)/mf(0))*p[4] + (mf(4)/mf(1))*p[5]*(np.log10(x)-X1LOG)
	for i in range (len(x)-1):
		if k2[i+1]*k2[i] < 0.0 or k3[i+1]*k3[i] < 0.0: 
			return 100.0+min(i,len(x)-i)
	return sqrt(((func44(p,x)-y)**2).sum()/float(len(x)))
	

# Read in the recombination line template and plot 
# fr = open('general_data/recline_template.txt', 'r')
# linea = fr.readline()
# lineb = fr.readline()
# fr.close()
# linea = linea.strip()
# lineb = lineb.strip()
# xtemp = linea.split()
# ytemp = lineb.split()
# xtemp = np.asfarray(xtemp)
# ytemp = np.asfarray(ytemp)
# ftemplate = interpolate.interp1d(xtemp, ytemp)
f = open('spec_to_get_template_6nov15.txt','r')
#f = open('general_data/spec_to_fit_10may14_2pt7e10_30jul14_to_archive_30jul2014_00000.txt', 'r') 
line1 = f.readline()		# read the frequencies in GHz
line1 = line1.strip()
x1 = line1.split()
x1 = np.asfarray(x1)
x1log = np.log10(x1)

# yt = ftemplate(x1)
# yt = np.asfarray(yt)

#print ' '
#print  'Plot the interpolated template'
#plt.plot(x1,yt,c='b')
#plt.xlabel('(f_GHz)')
#plt.ylabel('(T_recline antenna temperature)')
#plt.show()						

# accumulated_scaling_factor = []
# BF = []
irecord=0
#pdf_pages = PdfPages('plot_5sig_o4_no_lines_set3.pdf')
for line2 in f:
#		line2 = f.readline()
		print " "
		irecord += 1
		print " processing record number ",irecord
		line2 = line2.strip()
		y1 = line2.split()
		y1 = np.asfarray(y1)
		y1log = np.log10(y1)

# # Initial guess for optimization for H2 (lines exist!)			
# 		p00 = [np.log10(3.0),1.0,-2.5] 

# 		p1 = fmin(chisq1, p00, args=(x1,y1,yt),ftol=1.0e-20,maxiter=50000, maxfun=100000) 
# #		p2 = fmin(chisq1, p1, args=(x1,y1,yt),ftol=1.0e-20,maxiter=50000, maxfun=100000) 
# #		p1 = fmin(chisq1, p2, args=(x1,y1,yt),ftol=1.0e-20,maxiter=50000, maxfun=100000) 
# 		current_chisq = chisq1(p1,x1,y1,yt)
# 		print " "
# 		print "chisq1 = ",current_chisq
# 		print 10.0**p1[0], p1[1], p1[2]
# 		print " "

# 		p00 = [p1[0],p1[1],p1[2],0.0] 
# 		p1 = fmin(chisq2, p00, args=(x1,y1,yt),ftol=1.0e-20,maxiter=50000, maxfun=100000) 
# #		p2 = fmin(chisq2, p1, args=(x1,y1,yt),ftol=1.0e-20,maxiter=50000, maxfun=100000) 
# #		p1 = fmin(chisq2, p2, args=(x1,y1,yt),ftol=1.0e-20,maxiter=50000, maxfun=100000)  
# 		final_chisq = chisq2(p1,x1,y1,yt)
# 		print " "
# 		print "chisq2 = ",final_chisq
# 		print 10.0**p1[0], p1[1], p1[2], p1[3]
# 		print " "

# 		p00 = [p1[0],p1[1],p1[2],p1[3],0.0]
# 		p1 = fmin(chisq3, p00, args=(x1,y1,yt),ftol=1.0e-20,maxiter=50000, maxfun=100000) 
# #		p2 = fmin(chisq3, p1, args=(x1,y1,yt),ftol=1.0e-20,maxiter=50000, maxfun=100000) 
# #		p1 = fmin(chisq3, p2, args=(x1,y1,yt),ftol=1.0e-20,maxiter=50000, maxfun=100000)  
# 		final_chisq = chisq3(p1,x1,y1,yt)
# 		print " "
# 		print "chisq3 = ",final_chisq
# 		print 10.0**p1[0], p1[1], p1[2], p1[3], p1[4]
# 		print " "

# 		p00 = [p1[0],p1[1],p1[2],p1[3],p1[4],0.0]
# 		p1 = fmin(chisq4, p00, args=(x1,y1,yt),ftol=1.0e-20,maxiter=50000, maxfun=100000) 
# #		p2 = fmin(chisq4, p1, args=(x1,y1,yt),ftol=1.0e-20,maxiter=50000, maxfun=100000) 
# #		p1 = fmin(chisq4, p2, args=(x1,y1,yt),ftol=1.0e-20,maxiter=50000, maxfun=100000)  
# 		final_chisq = chisq4(p1,x1,y1,yt)
# 		print " "
# 		print "chisq4 = ",final_chisq
# 		print 10.0**p1[0], p1[1], p1[2], p1[3], p1[4], p1[5]
# 		print " "

# 		yfit2 = func4(p1,x1,yt)  # This fit is of foreground plus template
# 		yres2 = y1 - yfit2         # This residual is for H2
# 		yres2 = np.asfarray(yres2)

# Initial guess for optimization for H0 (lines do not exist!)			
		p00 = [np.log10(3.0),1.0,-2.5] 

		p1 = fmin(chisq11, p00, args=(x1,y1),ftol=1.0e-20,maxiter=50000, maxfun=100000) 
#		p2 = fmin(chisq1, p1, args=(x1,y1),ftol=1.0e-20,maxiter=50000, maxfun=100000) 
#		p1 = fmin(chisq1, p2, args=(x1,y1),ftol=1.0e-20,maxiter=50000, maxfun=100000) 
		current_chisq = chisq11(p1,x1,y1)
		print " "
		print "chisq1 = ",current_chisq
		print 10.0**p1[0], p1[1], p1[2]
		print " "

		p00 = [p1[0],p1[1],p1[2],0.0] 
		p1 = fmin(chisq22, p00, args=(x1,y1),ftol=1.0e-20,maxiter=50000, maxfun=100000) 
#		p2 = fmin(chisq2, p1, args=(x1,y1),ftol=1.0e-20,maxiter=50000, maxfun=100000) 
#		p1 = fmin(chisq2, p2, args=(x1,y1),ftol=1.0e-20,maxiter=50000, maxfun=100000)  
		final_chisq = chisq22(p1,x1,y1)
		print " "
		print "chisq2 = ",final_chisq
		print 10.0**p1[0], p1[1], p1[2], p1[3]
		print " "

		p00 = [p1[0],p1[1],p1[2],p1[3],0.0]
		p1 = fmin(chisq33, p00, args=(x1,y1),ftol=1.0e-20,maxiter=50000, maxfun=100000) 
#		p2 = fmin(chisq3, p1, args=(x1,y1),ftol=1.0e-20,maxiter=50000, maxfun=100000) 
#		p1 = fmin(chisq3, p2, args=(x1,y1),ftol=1.0e-20,maxiter=50000, maxfun=100000)  
		final_chisq = chisq33(p1,x1,y1)
		print " "
		print "chisq3 = ",final_chisq
		print 10.0**p1[0], p1[1], p1[2], p1[3], p1[4]
		print " "

		p00 = [p1[0],p1[1],p1[2],p1[3],p1[4],0.0]
		p1 = fmin(chisq44, p00, args=(x1,y1),ftol=1.0e-20,maxiter=50000, maxfun=100000) 
#		p2 = fmin(chisq4, p1, args=(x1,y1),ftol=1.0e-20,maxiter=50000, maxfun=100000) 
#		p1 = fmin(chisq4, p2, args=(x1,y1),ftol=1.0e-20,maxiter=50000, maxfun=100000)  
		final_chisq = chisq44(p1,x1,y1)
		print " "
		print "chisq4 = ",final_chisq
		print 10.0**p1[0], p1[1], p1[2], p1[3], p1[4], p1[5]
		print " "

		yfit = func44(p1,x1)  # This fit is of foreground 
		yres = y1 - yfit         # This residual is for H0
		yres = np.asfarray(yres)

		#ydiff = yres - yres2
# Final plots in linear scale
		plt.figure()
		plt.grid(c=0.8)
		plt.plot(x1,yres*1.0e9,linestyle='-',lw=2)
		plt.tick_params(axis='both', which='major', labelsize=28)
		plt.tick_params(axis='both', which='minor', labelsize=28)
		#plt.plot(x1,yres2,linestyle=':')
		plt.xlabel('Frequency [GHz]',fontsize=30)
		plt.ylabel('Recombination line signal [nK]',fontsize=30)
		plt.show()

# 		fig = plt.figure()
# 		plt.plot(x1,ydiff)
# 		plt.xlabel('Frequency (GHz)')
# 		plt.ylabel('yres-yres2')
# 		pdf_pages.savefig(fig)
# # Compute the variance in the residue

# 		yres_variance = np.var(yres)
# 		yres_allan_variance = 0.0
# 		for i in range(len(yres)-1):
# 			yres_allan_variance += (yres[i+1]-yres[i])**2
# 		yres_allan_variance /= (len(yres)-1)
# 		print "Sigma of the residuals in yres = ",sqrt(yres_variance),sqrt(yres_allan_variance/2.0)
# 		yres_allan_variance2 = 0.0
# 		yres_variance2 = np.var(yres2)
# 		for i in range(len(yres2)-1):
# 			yres_allan_variance2 += (yres2[i+1]-yres2[i])**2
# 		yres_allan_variance2 /= (len(yres2)-1)
# 		print "Sigma of the residuals in yres2 = ",sqrt(yres_variance2),sqrt(yres_allan_variance2/2.0)

# # Compute odds ratio ratio or Bayes factor

		
# 		print "Difference between residuals :",ydiff.sum()


# 		H0 = 1.0
# 		H2 = 1.0
# 		NN = np.size(yres)
# 		for i in range (len(yres)):
# 			H0 *= exp(-yres[i]**2/(2.0*yres_allan_variance/2.0))
# 			H2 *= exp(-yres2[i]**2/(2.0*yres_allan_variance2/2.0))
# #			H2 *= exp(-yres2[i]**2/(2.0*yres_allan_variance2/2.0)) # keep the variance identiical in both cases
# #			print " i H0 H2 : ",i,H0,H2


# 		print "AV 2 and AV 0 :", yres_allan_variance2,yres_allan_variance
# 		bayes_factor = ((yres_allan_variance2/yres_allan_variance)**(-0.5*NN))*H2/H0 ##ORIGINAL!!
# 		#bayes_factor = ((yres_allan_variance/yres_allan_variance2)**(-0.5*NN))*H0/H2
# 		# bayes_factor = H2/H0
# #		print " bayes_factor (H2/H0) = ", bayes_factor
# 		BF.append(bayes_factor)
		
# #BF = np.asfarray(BF)
# thefile = open('bayes_factors_5sig_100spec_o4_no_lines_set3.txt', 'w')
# for item in BF:
#   thefile.write("%s\n" % item)
# thefile.close()
# BF = np.asfarray(BF)
# print "Median of the Bayes Factor from",irecord,"spectra noise realizations is",np.median(BF)
# print "Mean of the Bayes Factor from",irecord,"spectra noise realizations is",np.mean(BF)
# pdf_pages.close()
