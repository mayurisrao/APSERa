import numpy as np
import scipy
from scipy.optimize import fmin, basinhopping
import copy
import sys 
from scipy import constants
sys.path.append('/home/pratush/Documents/saras3/utils')
from misc import *

def set_prior_BD(xx, tau, amp, cen, w, thresh):
    
    prior = 1e12
    
    if not((57<=cen<=83) and (1.4<=w<=30)):
        return prior

    #Compute Zero crossings in first derivative
    BD_signal      = T21(xx, tau, amp, cen, w)
    der1           = np.diff(BD_signal)
    zero_cross     = (der1[:-1] * der1[1:]) < 0
    no_zero_cross  = (zero_cross).sum()

    if no_zero_cross == 0:
        return prior

    zero_cross_ind = np.where(zero_cross==True)[0][0]
    zero_der_freq  = xx[zero_cross_ind+1]

    if no_zero_cross > 1:
        raise ValueError("More than 1 zero crossings found")
    if np.abs(zero_der_freq - cen)>1:
        raise ValueError("Zero derivative does not align with central frequency")

    signal_peak = BD_signal[zero_cross_ind+1]
    max_ind = len(xx) - (zero_cross_ind+1)

    # Test for increase in ampltitude from zero derivative
    for i in range(1, max_ind):
        pt1 = zero_cross_ind+1 - i
        pt2 = zero_cross_ind+1 + i

        if pt1 < 0 or pt1 >= len(xx):
            return prior

        val1 = np.abs(BD_signal[pt1] - signal_peak)
        val2 = np.abs(BD_signal[pt2] - signal_peak)

        per1 = np.abs(val1/signal_peak)
        per2 = np.abs(val2/signal_peak)

        if per1 >= 0.05 and per2 >= 0.05:
            prior = 1.0
            return prior

        #if val1 >= thresh and val2 >= thresh:
        #    prior = 1.0
        #    return prior

    return prior

def model_gauss(p, x, x_MHz, order, redshift, return_gauss=False):

    nterms_poly = order + 1

    p_poly  = p[0:nterms_poly]
    amp, cen, wid = p[nterms_poly:]

    if redshift:
        gaussian = gaussian_unity(amp, cen, wid, f_z(x_MHz))
    else:
        gaussian = gaussian_unity(amp, cen, wid, x_MHz)

    fg = 10**(np.polyval(p_poly,x))

    if return_gauss:
        model = gaussian
    else:
        model = fg + gaussian
    return model

def model_scale(p, x, order, BD):

    nterms_poly = order + 1

    p_poly  = p[0:nterms_poly]
    scale   = p[nterms_poly:]

    model = 10**(np.polyval(p_poly,x)) + scale*BD

    return model

def model_BD(p, x, x_MHz, order, return_BD=False):

    nterms_poly = order + 1

    p_poly  = p[0:nterms_poly]
    tau, amp, cen, w = p[nterms_poly:]
    BD = T21(x_MHz, tau, amp, cen, w)

    fg = 10**(np.polyval(p_poly,x))

    if return_BD:
        model = BD
    else:
        model = fg + BD
    return model

def model_ag(p, x, x_MHz, order, return_ag=False):

    nterms_poly = order + 1

    p_poly  = p[0:nterms_poly]
    amp, cen, wid1, wid2 = p[nterms_poly:]
    ag_gaussian = asym_gaussian(amp, cen, wid1, wid2, x_MHz)

    fg = 10**(np.polyval(p_poly,x))

    if return_ag:
        model = ag_gaussian
    else:
        model = fg + ag_gaussian
    return model

def model_gauss_1D(p, x, x_MHz, order, param_fits, return_gauss=False):

    nterms_poly = order + 1

    p_poly  = p[0:nterms_poly]

    if len(p) != (nterms_poly+1):
        raise ValueError("Wrong number of patameters.")

    if param_fits[0] == 'amp':
        cen, wid = param_fits[1:]
        amp = p[-1]

    elif param_fits[0] == 'cen':
        amp, wid = param_fits[1:] 
        cen = p[-1]

    elif param_fits[0] == 'wid':
        amp, cen = param_fits[1:]
        wid = p[-1]
    else:
        raise ValueError("param_fits not recognized. Choose between 'amp', 'cen' and 'wid'")

    gaussian = gaussian_unity(amp, cen, wid, x_MHz)
    fg = 10**(np.polyval(p_poly,x))

    if return_gauss:
        model = gaussian
    else:
        model = fg + gaussian

    return model

def model_BD_1D(p, x, x_MHz, order, param_fits, return_BD=False):

    nterms_poly = order + 1

    p_poly  = p[0:nterms_poly]

    if len(p) != (nterms_poly+1):
        raise ValueError("Wrong number of patameters.")

    tau, cen, w = param_fits
    amp = p[-1]

    BD = T21(x_MHz, tau, amp, cen, w)
    fg = 10**(np.polyval(p_poly,x))

    if return_BD:
        model = BD
    else:
        model = fg + BD
    return model

def model_log_log(p, x):

    model = 10**(np.polyval(p, x))
    return model 

def check_smoothness(p, x):

    nterms = len(p)
    c1 = np.poly1d(p)
    MM = 2

    if nterms <= 3:
        return True

    for i in range(MM, (nterms-1)):  # index 2 here implies quadratic is allowed

        c2 = np.polyder(c1,i)
        der = c2(x)
        
        no_zero_crossing = ((der[:-1] * der[1:]) < 0).sum()
        if no_zero_crossing > 0:
            return False
    return True

def model_sine(p, x, x_MHz, order=None, \
        return_sine=False, constraints=False, \
        constraints_dict = None, fit_sine_only=False,\
        smooth = False, fit_fg_only = False, \
        efficiency_weighted=False, efficiency=None):

    wt_sine = np.ones(len(x))

    if not(fit_sine_only):

        if order is None:
            raise ValueError("Order is needed")
    
        nterms_poly = order + 1
        
        p_poly = p[0:nterms_poly]
        fg = model_log_log(p_poly, x)
        poly_smooth = check_smoothness(p_poly, x)

        if smooth and not(poly_smooth):
            return 1e10
        
        if fit_fg_only:
            return fg

        if (len(p) != (nterms_poly+3)):#polynomial + 3 terms for sine
            raise ValueError("Wrong number of patameters")
        
        amp, length, phase = p[nterms_poly:]

    else:
        amp, length, phase = p
        fg = np.zeros(len(x))

    if constraints:
        if not((constraints_dict['amp'][0]    <= amp    <= constraints_dict['amp'][1]) and 
               (constraints_dict['length'][0] <= length <= constraints_dict['length'][1]) and
               (constraints_dict['phase'][0]  <= phase  <= constraints_dict['phase'][1])):
            return 1e10

    if efficiency_weighted:
        if efficiency is None:
            raise ValueError("efficiency is needed for the weight")
        wt_sine = efficiency 
    
    sine_sig = wt_sine * amp * np.sin(2*np.pi * (2*length/constants.c) * (x_MHz*1e6) + phase) 

    if return_sine:
        model = sine_sig
    else:
        model = fg + sine_sig
    return model

def chisq_poly_gauss(p, *args):

    x, x_MHz, y, e, order, redshift = args
    
    nterms_poly = order + 1
    p_poly  = p[0:nterms_poly]
    amp, cen, wid = p[nterms_poly:]
   
    if wid < 1.4:
        return 1e10
    
    wtt     = 1/e**2
    model_data = model_gauss(p, x, x_MHz, order, redshift)
    cc = np.sum(wtt*(y-model_data)**2)/np.sum(wtt)

    return cc

def chisq_poly_scale(p, *args):

    x, y, e, order, BD, ygchannel0, resol = args
    
    model_data = model_scale(p, x, order, BD)
    res        = y - model_data
    wtt        = 1/e**2
    
    rms_list, yres, _ = \
            smooth_res_MM([2, resol], ygchannel0, res, e)
    ac_sig = yres[-1]*rms_list[-1]/rms_list[0]

    cc = np.sum(wtt*(ac_sig)**2)/np.sum(wtt)

    return cc

def chisq_poly_scale_leastsq(p, *args):

    x, y, e, order, BD, ygchannel0, resol = args
    
    model_data = model_scale(p, x, order, BD)
    res        = y - model_data
    wtt        = 1/e
    
    rms_list, yres, _ = \
            smooth_res_MM([2, resol], ygchannel0, res, e)
    ac_sig = yres[-1]*rms_list[-1]/rms_list[0]

    cc = np.abs(wtt*(ac_sig))

    return cc

def chisq_poly_BD(p, *args):

    x, x_MHz, y, e, order = args
    
    nterms_poly = order + 1
    tau, amp, cen, w = p[nterms_poly:]
    prior = set_prior_BD(x_MHz, tau, amp, cen, w, 0.033)

    wtt        = 1/e**2
    model_data = model_BD(p, x, x_MHz, order)
    cc         = prior*np.sum(wtt*(y-model_data)**2)/np.sum(wtt)
    return cc

def chisq_poly_ag(p, *args):

    x, x_MHz, y, e, order = args
    
    nterms_poly = order + 1
    wtt        = 1/e**2
    model_data = model_ag(p, x, x_MHz, order)

    p_poly  = p[0:nterms_poly]
    amp, cen, wid1, wid2 = p[nterms_poly:]

    if not((57 <= cen <= 82) and (1.4 <= wid1) and (1.4 <= wid2)):
        return 1e10
    else:
        cc         = np.sum(wtt*(y-model_data)**2)/np.sum(wtt)
        return cc

def chisq_poly_gauss_1D(p, *args):

    x, x_MHz, y, e, order, param_fits = args # param_fits = ['amp'/'cen'/'wid', a, b]
    wtt     = 1/e**2
    model_data = model_gauss_1D(p, x, x_MHz, order, param_fits)
    cc = np.sum(wtt*(y-model_data)**2)/np.sum(wtt)
    return cc

def chisq_poly_BD_1D(p, *args):

    x, x_MHz, y, e, order, param_fits = args 
    wtt     = 1/e**2
    model_data = model_BD_1D(p, x, x_MHz, order, param_fits)
    cc = np.sum(wtt*(y-model_data)**2)/np.sum(wtt)
    return cc

def chisq_poly_sine(p, *args):

    x, x_MHz, y, e, order, constraints, constraints_dict, \
            fit_sine_only, smooth, fit_fg_only,\
            efficiency_weighted, efficiency= args 

    wtt     = 1/e**2
    model_data = model_sine(p, x, x_MHz, order, constraints=constraints, constraints_dict=constraints_dict, \
            fit_sine_only=fit_sine_only, smooth = smooth, fit_fg_only = fit_fg_only,\
            efficiency_weighted=efficiency_weighted, efficiency=efficiency)

    cc = np.sum(wtt*(y-model_data)**2)/np.sum(wtt)
    return cc

def do_fit_gauss(x1, x_MHz, y1, e1, order, redshift, p00):

    stepsize         = 1e-5
    T                = 1e-5
    niter            = 1
    seed             = 1

    args             = (x1, x_MHz, y1, e1, order, redshift)

    OPTIONS={'ftol':1e-10, 'xtol':1e-10, 'maxiter':1e5, 'maxfev':1e5}
    minimizer_kwargs={'method':'Nelder-Mead', 'args':args, 'options':OPTIONS}

    for i in range(100):
        pout = basinhopping(chisq_poly_gauss, p00,minimizer_kwargs=minimizer_kwargs,T=T, \
                              stepsize=stepsize, niter=niter, seed=seed)
        p00 = pout.x

    pout = basinhopping(chisq_poly_gauss, p00,minimizer_kwargs=minimizer_kwargs,T=T, \
                              stepsize=stepsize, niter=niter, seed=seed)
    p1 = pout.x

    print ("Coefficients : {} \n".format(p1))
    print("Residual rms : {}".format(np.sqrt(chisq_poly_gauss(p1, *args))))

    p1_basin    = copy.deepcopy(p1)
    y_fit_basin = model_gauss(p1_basin, x1, x_MHz, order, redshift)
    yres_basin  = y1 - y_fit_basin

    return p1_basin, yres_basin, y_fit_basin

def do_fit_scale(x1, y1, e1, order, BD,  ygchannel0, resol, p00):

    stepsize         = 1e-5
    T                = 1e-5
    niter            = 1
    seed             = 1

    args             = (x1, y1, e1, order, BD, ygchannel0, resol)

    OPTIONS={'ftol':1e-10, 'xtol':1e-10, 'maxiter':1e5, 'maxfev':1e5}
    minimizer_kwargs={'method':'Nelder-Mead', 'args':args, 'options':OPTIONS}

    for i in range(100):
        pout = basinhopping(chisq_poly_scale, p00,minimizer_kwargs=minimizer_kwargs,T=T, \
                              stepsize=stepsize, niter=niter, seed=seed)
        p00 = pout.x

    pout = basinhopping(chisq_poly_scale, p00,minimizer_kwargs=minimizer_kwargs,T=T, \
                              stepsize=stepsize, niter=niter, seed=seed)
    p1 = pout.x

    print ("Coefficients : {} \n".format(p1))
    print("Residual rms : {}".format(np.sqrt(chisq_poly_scale(p1, *args))))

    p1_basin    = copy.deepcopy(p1)
    y_fit_basin = model_scale(p1_basin, x1, order, BD)
    yres_basin  = y1 - y_fit_basin

    return p1_basin, yres_basin, y_fit_basin


def do_fit_BD(x1, x_MHz, y1, e1, order, p00):

    stepsize         = 1e-5
    T                = 1e-5
    niter            = 1
    seed             = 1

    args             = (x1, x_MHz, y1, e1, order)

    OPTIONS={'ftol':1e-10, 'xtol':1e-10, 'maxiter':1e5, 'maxfev':1e5}
    minimizer_kwargs={'method':'Nelder-Mead', 'args':args, 'options':OPTIONS}

    for i in range(100):
        pout = basinhopping(chisq_poly_BD, p00,minimizer_kwargs=minimizer_kwargs,T=T, \
                              stepsize=stepsize, niter=niter, seed=seed)
        p00 = pout.x

    pout = basinhopping(chisq_poly_BD, p00,minimizer_kwargs=minimizer_kwargs,T=T, \
                              stepsize=stepsize, niter=niter, seed=seed)
    p1 = pout.x

    print ("Coefficients : {} \n".format(p1))
    print("Residual rms : {}".format(np.sqrt(chisq_poly_BD(p1, *args))))

    p1_basin    = copy.deepcopy(p1)
    y_fit_basin = model_BD(p1_basin, x1, x_MHz, order)
    yres_basin  = y1 - y_fit_basin

    return p1_basin, yres_basin, y_fit_basin

def do_fit_ag(x1, x_MHz, y1, e1, order, p00):

    order            = 6
    stepsize         = 1e-5
    T                = 1e-5
    niter            = 1
    seed             = 1

    args             = (x1, x_MHz, y1, e1, order)

    OPTIONS={'ftol':1e-10, 'xtol':1e-10, 'maxiter':1e5, 'maxfev':1e5}
    minimizer_kwargs={'method':'Nelder-Mead', 'args':args, 'options':OPTIONS}

    for i in range(100):
        pout = basinhopping(chisq_poly_ag, p00,minimizer_kwargs=minimizer_kwargs,T=T, \
                              stepsize=stepsize, niter=niter, seed=seed)
        p00 = pout.x

    pout = basinhopping(chisq_poly_ag, p00,minimizer_kwargs=minimizer_kwargs,T=T, \
                              stepsize=stepsize, niter=niter, seed=seed)
    p1 = pout.x

    print ("Coefficients : {} \n".format(p1))
    print("Residual rms : {}".format(np.sqrt(chisq_poly_ag(p1, *args))))

    p1_basin    = copy.deepcopy(p1)
    y_fit_basin = model_ag(p1_basin, x1, x_MHz, order)
    yres_basin  = y1 - y_fit_basin

    return p1_basin, yres_basin, y_fit_basin

def do_fit_gauss_1D(x1, x_MHz, y1, e1, order, param_fits, p00):

    stepsize         = 1e-5
    T                = 1e-5
    niter            = 1
    seed             = 1

    args             = (x1, x_MHz, y1, e1, order, param_fits)

    OPTIONS={'ftol':1e-10, 'xtol':1e-10, 'maxiter':1e5, 'maxfev':1e5}
    minimizer_kwargs={'method':'Nelder-Mead', 'args':args, 'options':OPTIONS}

    for i in range(100):
        pout = basinhopping(chisq_poly_gauss_1D, p00,minimizer_kwargs=minimizer_kwargs,T=T, \
                              stepsize=stepsize, niter=niter, seed=seed)
        p00 = pout.x

    pout = basinhopping(chisq_poly_gauss_1D, p00,minimizer_kwargs=minimizer_kwargs,T=T, \
                              stepsize=stepsize, niter=niter, seed=seed)
    p1 = pout.x

    print ("Coefficients : {} \n".format(p1))
    print("Residual rms : {}".format(np.sqrt(chisq_poly_gauss_1D(p1, *args))))

    p1_basin    = copy.deepcopy(p1)
    y_fit_basin = model_gauss_1D(p1_basin, x1, x_MHz, order, param_fits)
    yres_basin  = y1 - y_fit_basin

    return p1_basin, yres_basin, y_fit_basin

def do_fit_BD_1D(x1, x_MHz, y1, e1, order, param_fits, p00):

    stepsize         = 1e-5
    T                = 1e-5
    niter            = 1
    seed             = 1

    args             = (x1, x_MHz, y1, e1, order, param_fits)

    OPTIONS={'ftol':1e-10, 'xtol':1e-10, 'maxiter':1e5, 'maxfev':1e5}
    minimizer_kwargs={'method':'Nelder-Mead', 'args':args, 'options':OPTIONS}

    for i in range(100):
        pout = basinhopping(chisq_poly_BD_1D, p00,minimizer_kwargs=minimizer_kwargs,T=T, \
                              stepsize=stepsize, niter=niter, seed=seed)
        p00 = pout.x

    pout = basinhopping(chisq_poly_BD_1D, p00,minimizer_kwargs=minimizer_kwargs,T=T, \
                              stepsize=stepsize, niter=niter, seed=seed)
    p1 = pout.x

    print ("Coefficients : {} \n".format(p1))
    print("Residual rms : {}".format(np.sqrt(chisq_poly_BD_1D(p1, *args))))

    p1_basin    = copy.deepcopy(p1)
    y_fit_basin = model_BD_1D(p1_basin, x1, x_MHz, order, param_fits)
    yres_basin  = y1 - y_fit_basin

    return p1_basin, yres_basin, y_fit_basin

def do_fit_sine(x1, x_MHz, y1, e1, p00, order=None, iterations=100, \
        constraints=False, constraints_dict = None, fit_sine_only=False,\
        smooth=False, fit_fg_only = False, efficiency_weighted=False, efficiency=None):

    stepsize         = 1e-5
    T                = 1e-5
    niter            = 1
    seed             = 1

    args             = (x1, x_MHz, y1, e1, order, constraints, constraints_dict, \
                        fit_sine_only, smooth, fit_fg_only,\
                        efficiency_weighted, efficiency)

    OPTIONS={'ftol':1e-10, 'xtol':1e-10, 'maxiter':1e5, 'maxfev':1e5}
    minimizer_kwargs={'method':'Nelder-Mead', 'args':args, 'options':OPTIONS}

    for i in range(iterations):
        pout = basinhopping(chisq_poly_sine, p00, minimizer_kwargs=minimizer_kwargs,T=T, \
                              stepsize=stepsize, niter=niter, seed=seed)
        p00 = pout.x

    pout = basinhopping(chisq_poly_sine, p00, minimizer_kwargs=minimizer_kwargs,T=T, \
                              stepsize=stepsize, niter=niter, seed=seed)
    p1 = pout.x

    print ("Coefficients : {} \n".format(p1))
    print("Residual rms : {}".format(np.sqrt(chisq_poly_sine(p1, *args))))

    p1_basin    = copy.deepcopy(p1)

    y_fit_basin = model_sine(p1_basin, x1, x_MHz, order, \
            constraints=constraints, constraints_dict=constraints_dict,\
            fit_sine_only=fit_sine_only,smooth=smooth, fit_fg_only=fit_fg_only,\
            efficiency_weighted=efficiency_weighted, efficiency=efficiency)

    yres_basin  = y1 - y_fit_basin

    return p1_basin, yres_basin, y_fit_basin
