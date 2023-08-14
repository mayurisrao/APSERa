import numpy as np
import matplotlib.pyplot as plt
import sys, getopt
import os
import time
from scipy import interpolate
import itertools
from scipy import signal
import copy
from scipy.optimize import fmin, basinhopping
from math import exp, expm1, sqrt, sin, cos
from math import factorial as mf
import scipy.optimize as op
import scipy.constants
import csv
from prettytable import PrettyTable
from scipy.optimize import leastsq
import warnings
import pickle

PI = np.pi
CVEL = scipy.constants.c

###############################################

#Rescales arr list values into range min1 to max1
def rescale(arr, min1, max1, log):
    arr = np.asfarray(arr)
    if log:
        arr = np.log10(arr)
    min_arr = np.amin(arr)
    max_arr = np.amax(arr)
    arr_sc  = ((max1 - min1))*((arr) - float(min_arr))/(max_arr - min_arr) + min1
    return arr_sc	

###############################################

def select_freq_1d(x1, low, high): 
    prl = 0
    prh = 0
    for i in range(len(x1)):
        if x1[i]<=low:
            i_low=i
            prl = 1
        if x1[i]>=high:
            i_high=i
            prh = 1
            break
    if prl==1 and prh==0:
        i_high = len(x1)
        warnings.warn("Could not find the higher limit. Using the miaximum")
    if prl==0 and prh==1:
        i_low = 0	
        warnings.warn("Could not find the lower limit. Using the minimum")
    if prl==0 and prh==0:
        warnings.warn("Could not find the limits. Going to use the full array")
        i_low = 0
        i_high = len(x1)
    return i_low, i_high

##########################################################################
func_poly = lambda p, x    : (np.polyval(p,x))


def chisq_poly (p, *args):
    x, y, e, flag, signal21, domain, additive, joint, do_opt, met, smooth = args
    wtt = 1/e**2
    
    if domain=='log_log':
        
        if additive and (not joint):
            model = 10.0**func_poly(p[:-1],x) + p[-1]
            c1=np.poly1d(p[:-1])
            nterms = len(p[:-1])
        
        elif additive and joint:
            model = 10.0**func_poly(p[:-2],x) + p[-2] + p[-1]*signal21 
            c1=np.poly1d(p[:-2])
            nterms = len(p[:-2])
        
        elif not(additive) and joint:
            model = 10.0**func_poly(p[:-1],x) + p[-1]*signal21
            c1=np.poly1d(p[:-1])
            nterms = len(p[:-1])
        
        else:
            model = 10.0**func_poly(p,x)
            c1=np.poly1d(p)
            nterms = len(p)
    
    elif joint:
        model = func_poly(p[:-1],x) + p[-1]*signal21
        c1=np.poly1d(p[:-1])
        nterms = len(p[:-1])
    
    else:
        model = func_poly(p,x)
        c1=np.poly1d(p)
        nterms = len(p)

    if not(do_opt):
        return model
    
    if nterms <= 3:
        return ( ((model[flag]-y[flag])**2  )*wtt[flag] ).sum() / wtt[flag].sum() 

    if smooth:
    ### Find the highest order polynomial with at most one zero-crossing
        MM = 2
        for j in range(MM, nterms-1):  # check derivative polynomials for each fit 
            c2 = np.polyder(c1,j)
            der = c2(x)
            nzdetect = 0
            for k in range((len(der)-1)):  # for each derivative, scan through the x axis
                if (der[k]*der[k+1])<0:    # if a ZC is found,
                    nzdetect += 1
                    if nzdetect >= 2:
                        if met=='fmin':
                            return 1000.0+min(j,len(x)-j)
                        elif met=='leastsq':
                            return np.ones(len(y)) * (1000.0+min(j,len(x)-j))
                        else:
                            raise ValueError("Specified method not implemeneted!")
    if met=='fmin':
        return ( ((model[flag]-y[flag])**2  )*wtt[flag] ).sum() / wtt[flag].sum() 
    elif met=='leastsq':
        return wtt[flag]*np.abs(y[flag]-model[flag])
    else:
        raise ValueError("Specified method not implemeneted!")
##########################################################################

def ms_fit_inf(frq0, dat0, **kwargs):

    '''

    Performs a maximally smooth function fit in the given domain and x-range 
    Has the capability for fitting an overall additive, and jointly fitting with 21-cm signal.

    Returns: best fit parameters, selected x, min index of x, max index of x, fit curve, residuals 

    Parameters:

    frq0 : x-values, 1-D array 
    dat0 : y-values, 1-D array
    flag : Flagged frequency channels. Following miriad convention: 1 is good, and 0 is bad
    error: error on y-values, optional, default is 1. (Weight used in fitting is 1/error**2)
    xtol : XTOL for minimization, default is 1e-6
    ftol : FTOL for minimization, default is 1e-6
    domain : Domain to fit in. Possible options are 'lin_lin', 'log_lin' and 'log_log'. Default is 'lin_lin'
    maxiter : Maximum iterations for minimization, default is 1e5
    maxfev  : Maximum function evaluations for minimization, default is 1e5
    algo : Algorithm for minimization. Default is 'Nelder-Mead'
    norder : Order of Maximally Smooth Function. Default is 6
    temp : Temperature parameter for Basinhopping. Default is 0.1
    stepsize : Stepsize parameter for Basinhopping. Deafult is 0.1
    niter : Number of iterations within Basinhopping. Default is 1
    basin_iter_loops : Number of times Basinhopping would be run. Default is 10
    xmin : Minimum frequency to fit for
    xmax : Maximum frequency to fit for
    rescale : Whether to rescale x values to -1 to 1. Sometimes that gives better fit. Deafult is True.
    seed : Seed for basinhopping. Default is 1. 
    additive: Fit for an overall constant. Considered only when log_log domain is used. Default is True.
    joint: Try forward modeling with scale factor approach. Default is False.
    signal: 21-cm signal for forward modeling. Default is False.
    add_ini: initial guess for additive. Considered only when domain is log_log
    op_file: If parameters need to be saved
    op_name: Name of the output file
    verbose: If table of all the parameters need to listed
    smooth: If it is MS fit. Default is True, else it is unconstrained fit
    ret_int: Return intermediate results in a dictionary. Keys are number of terms, and not the order

    '''

    pd = {}
    
    #Initialize default parameters
    pd['flag'], pd['error'], pd['xtol'], pd['ftol'], pd['domain'] = np.ones(len(frq0)), np.ones(len(frq0)), 1e-6, 1e-6, 'lin_lin'
    pd['maxiter'], pd['maxfev'], pd['algo'], pd['norder'] = 1e5, 1e5, 'Nelder-Mead', 6
    pd['temp'], pd['stepsize'], pd['niter'], pd['basin_iter_loops'] = 0.1, 0.1, 1, 10
    pd['xmin'], pd['xmax'], pd['rescale'], pd['seed'] , pd['additive']= frq0[0], frq0[-1], True, 1, True
    pd['joint'], pd['signal'], pd['add_ini'], pd['op_file'], pd['op_name'], pd['verbose'] , pd['smooth'] = False, [False], 0.0, True, "output_file.pickle", True, True
    pd['ret_int'] = False

    all_keys = ['flag','error', 'xtol', 'ftol', 'domain', 'maxiter', 'maxfev', 'algo', 'norder', 'temp',\
            'stepsize','niter','basin_iter_loops','xmin','xmax','rescale','seed','additive','joint','signal','add_ini', 'op_file', \
            'op_name','verbose', 'smooth', 'ret_int']

    #Update the parameters depedning upon the arguments
    for key in kwargs:
        if key in all_keys:
            pd[key] = kwargs[key]

    #Test for consistency
    if pd['domain']=='log_log' and not(pd['additive']):
        warnings.warn("log-log fitting generally does not fit to overall DC. You might want to enable additive")
    if pd['joint'] and not('signal' in [key for key  in kwargs]):
        raise ValueError ("Signal needs to be supplied while joint fitting!")
    
    #Print the final list of parameters
    t = PrettyTable(['Parameter', 'Value'])
    for key in pd:
        item = pd[key]
        if np.size(item)>1:
            item = str(len(item)) + " elements"
        t.add_row([key, item])
    
    if pd['verbose']: 
        print(t)

    imin, imax = select_freq_1d(frq0, pd['xmin'], pd['xmax'])
    frq0 = frq0[imin:imax]
    y1   = dat0[imin:imax]
    e1   = pd['error'][imin:imax]
    flag = np.array(pd['flag'][imin:imax], dtype='bool')

    if len(pd['signal'])>1:
        sig_ss = pd['signal'][imin:imax]
    else:
        sig_ss = False

    if pd['domain'] in {'log_log', 'log_lin'}:
        if pd['rescale']:
            x1 = np.asfarray(rescale(frq0,-1.0,+1.0, True))  # rescale x values to be from -1 to 1
        else:
            x1 = np.log10(frq0)

    if pd['domain'] == 'lin_lin':
        if pd['rescale']:
            x1 = np.asfarray(rescale(frq0,-1.0,+1.0, False))  # rescale x values to be from -1 to 1
        else:
            x1 = copy.deepcopy(frq0)

    NN=2
    if pd['domain']=='log_log':
        z = np.polyfit((x1), np.log10(y1), NN)			
        if pd['additive']:
            p00 = np.concatenate((z,np.asarray([pd['add_ini']])))
        else:
            p00 = copy.deepcopy(z)
    else:
        z = np.polyfit((x1), (y1), NN)			
        p00 = copy.deepcopy(z)
        
    #OPTIONS DECIDES OPTIMIZATION PARAMETER FOR MINIMIZATION
    do_opt = True
    met = 'fmin'
    args = (x1,y1,e1,flag,sig_ss,pd['domain'], pd['additive'], pd['joint'], do_opt, met, pd['smooth'])
    OPTIONS={'ftol':pd['ftol'], 'xtol': pd['xtol'], 'maxiter':pd['maxiter'], 'maxfev':pd['maxfev']}

    #Check for joint fitting mode
    if pd['joint']:
        p00 = np.append(p00,0.0)
    
    minimizer_kwargs={'method':pd['algo'], 'args':args, 'options':OPTIONS}

    param_inter = {}

    for i in range(0, pd['norder']-NN+1):

        if pd['verbose']:
            print("\nPolynomial Order : {}".format(i+2))

        for i in range (pd['basin_iter_loops']):
            pout = basinhopping(chisq_poly, p00,minimizer_kwargs=minimizer_kwargs,T=pd['temp'], \
                                  stepsize=pd['stepsize'], niter=pd['niter'], seed=pd['seed'])
            p00 = pout.x	

        pout = basinhopping(chisq_poly, p00,minimizer_kwargs=minimizer_kwargs,T=pd['temp'], \
                                  stepsize=pd['stepsize'], niter=pd['niter'], seed=pd['seed'])
        p1 = pout.x	

        param_inter[len(p1)] = p1 

        if pd['verbose']:
            print ("Coefficients : {} \n".format(p1)) 
            print("Residual rms : {}".format(np.sqrt(chisq_poly(p1, *args))))
        p00=np.array(np.concatenate(([0.0], p1)), dtype='float128')

    p1_basin = copy.deepcopy(p1)

    if pd['op_file']:
        with open(pd['op_name']+'.pickle', 'wb') as handle:
            pickle.dump(p1_basin, handle)
   
    ############################################### 
    #Change the optimization method, just for errors.
    #Results should NOT change

    args_list = list(args)
    args_list[-2]="leastsq"
    args = tuple(args_list)

    p_leastsq_full = leastsq(chisq_poly, p1_basin, args=args, \
            ftol=pd['ftol'], xtol=pd['xtol'], maxfev = int(pd['maxfev']), full_output=True)
    p_joint_leastsq        = p_leastsq_full[0]

    
    if pd['verbose']:
        print("Final set of coefficients (leastsq): {}".format(p_joint_leastsq))
        if(np.allclose(p_joint_leastsq,p1_basin)):
            print("\nLeastsq gives similar values to that of Basinhopping\n")
        else:
            warnings.warn("\nLeastsq gives different basinhopping. Fit might not have converged, don't trust the errors!\n")
        
    ############################################### 
    #Turn off optimization, change back to fmin, and return model
    args_list = list(args)
    args_list[-2]="fmin"
    args_list[-3]=False
    args = tuple(args_list)

    y_fit_basin = chisq_poly(p1_basin,*args)
    yres_basin  = y1 - y_fit_basin

    s_sq      	           = (len(e1)/(len(e1)-(pd['norder']+1)))* np.sum((1.0/e1**2) * (yres_basin)**2)/\
                             (np.sum(1.0/e1**2))
    err_p1 = False
    try:
        cov_x = s_sq * p_leastsq_full[1]
        cov_x = np.asarray(cov_x)
        err_p1 = np.sqrt(np.abs(np.diag(cov_x)))
    except:
        if pd['verbose']:
            print("Unable to compute covariance matrix!")

    if pd['ret_int']:
        return p1_basin, x1, imin, imax, y_fit_basin, yres_basin, err_p1, param_inter
    else:
        return p1_basin, x1, imin, imax, y_fit_basin, yres_basin, err_p1
