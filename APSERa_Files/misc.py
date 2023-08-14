import numpy as np
from scipy import signal
import aipy as a
import warnings

def rescale(arr, min1, max1, log):
    arr = np.asfarray(arr)
    if log:
        arr = np.log10(arr)
    min_arr = np.amin(arr)
    max_arr = np.amax(arr)
    arr_sc  = ((max1 - min1))*((arr) - float(min_arr))/(max_arr - min_arr) + min1
    return arr_sc


get_FWHM = lambda sigma: 2*np.sqrt(2*np.log(2)) * sigma

get_sigma = lambda FWHM: FWHM/(2*np.sqrt(2*np.log(2)))

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, np.sqrt(variance))

def smooth_res(nsmooth, ygchannel0, y_res, e1):
    yres0_smooth_all = []
    rms_list         = []
    smoothing = []
    avg0, rms0       = weighted_avg_and_std(y_res, 1/e1**2)
    for mm in range (1,nsmooth):
        MM           = np.power(2,mm)  #  MM=32 gives approx 1 MHz FWHM smoothing
        win          = signal.hann(2*MM+1)
        yres0_smooth = []
        for i in range (len(ygchannel0)):
            i0     = ygchannel0[i]
            i1     = i0-MM
            i2     = i0+MM
            sumwt  = 0.0
            sumval = 0.0
            for j in range (len(ygchannel0)):
                ic      = ygchannel0[j]
                if ic < i1 or ic > i2 : continue
                jwin    = (ic-i0) + MM
                cwt     = win[jwin]/(e1[j]**2)
                sumwt  += cwt
                sumval += cwt*y_res[j]
            if sumwt != 0.0 : yres0_smooth.append(sumval/sumwt)
            else : yres0_smooth.append(0.0)

        yres0_smooth  = np.asfarray(yres0_smooth)
        avg1,rms1     = weighted_avg_and_std(yres0_smooth, 1/e1**2)
        yres0_smooth *= rms0/rms1
        rms_list.append(rms1)
        smoothing.append(float(MM)*250.0/8192.0)
        yres0_smooth_all.append(yres0_smooth)
    return np.asarray(rms_list), np.asarray(yres0_smooth_all), np.asarray(smoothing)

def smooth_res_MM(nchan_arr, ygchannel0, y_res, e1, return_error = False):
    yres0_smooth_all = []
    rms_list         = []
    smoothing = []
    avg0, rms0       = weighted_avg_and_std(y_res, 1/e1**2)
    for MM in nchan_arr:
        #  MM=32 gives approx 1 MHz FWHM smoothing
        win          = signal.hann(2*MM+1)
        yres0_smooth = []
        weights   = []
        for i in range (len(ygchannel0)):
            i0     = ygchannel0[i]
            i1     = i0-MM
            i2     = i0+MM
            sumwt  = 0.0
            sumval = 0.0
            for j in range (len(ygchannel0)):
                ic      = ygchannel0[j]
                if ic < i1 or ic > i2 : continue
                jwin    = (ic-i0) + MM
                cwt     = win[jwin]/(e1[j]**2)
                sumwt  += cwt
                sumval += cwt*y_res[j]
            if sumwt != 0.0 : 
                yres0_smooth.append(sumval/sumwt)
                weights.append(np.sqrt(1.0/sumwt))
            else : 
                yres0_smooth.append(0.0)
                weights.append(0.0)

        yres0_smooth  = np.asfarray(yres0_smooth)
        avg1,rms1     = weighted_avg_and_std(yres0_smooth, 1/e1**2)
        yres0_smooth *= rms0/rms1
        rms_list.append(rms1)
        smoothing.append(float(MM)*250.0/8192.0)
        yres0_smooth_all.append(yres0_smooth)
    if return_error:
        return np.asarray(rms_list), np.asarray(yres0_smooth_all), np.asarray(smoothing), weights
    else:
        return np.asarray(rms_list), np.asarray(yres0_smooth_all), np.asarray(smoothing)

def read_flagpyavg_file(file_in, PATH,  DATASIZE, chn_low, chn_hgh, channel): #taken from Ravi's MS func with 1 infl

    uvi = a.miriad.UV(PATH + file_in)
    #uvi.select('polarization', -3, 0,  include=True)
    # The first read is just to find the number of records to open appropriately sized arrays
    no_records=0
    for preamble, data in uvi.all():
            no_records+=1  # To get the number of records
    uvi.rewind()


    # Initialize 2D arrays to store data, fit and flags
    # These are for writing to MIRIAD

    p_sn=[]
    p_pol=[]
    p_data = np.zeros((no_records, DATASIZE), dtype=np.complex64)  # real+imag part
    p_flag = np.zeros((no_records, DATASIZE), dtype=np.int64) # integer flags

    i=0
    for preamble, data in uvi.all():
            p_sn.append(i) # Serial number for each record; starts with 0
            p_pol.append(uvi['pol']) # polarization list
            flags=np.logical_not(data.mask) # Miriad flags: 1 is gooddata, Mask: False is gooddata
            p_data[i]=data.data #Data array
            p_flag[i]=flags # Miriad type Flags array
            i += 1

    #To get the starting point right and set store_start to be the starting index
    for i in p_sn:
            if (p_pol[i]==-1 and p_pol[i+1]==-2 and p_pol[i+2]==-3 \
                            and p_pol[i+3]==-4):
                    store_start=i
                    break


    # Extract the indices for required frequency range chn_low to chn_hgh using first record
    ydata  = p_data[store_start].real # Real part of cross correlation
    yflag  = p_flag[store_start]      # Flag for each record
    yerror = p_data[store_start+3].real # STD for each record
    _, _, _, _, ilow, ihigh = select_freq_1d(channel, ydata, yerror, yflag, chn_low, chn_hgh)


    lengthfromstart = len(p_data) - store_start
    residuallength = (len(p_data) - store_start)%4
    SPEC_END   = (lengthfromstart - residuallength)//4


    ychannel  = channel[ilow:ihigh+1]
    ydata     = np.zeros((SPEC_END, len(ychannel)), dtype=np.float64)
    yres      = np.zeros((SPEC_END, len(ychannel)), dtype=np.float64)
    yfit      = np.zeros((SPEC_END, len(ychannel)), dtype=np.float64)
    yflag     = np.zeros((SPEC_END, len(ychannel)), dtype=np.int32)
    yerror    = np.zeros((SPEC_END, len(ychannel)), dtype=np.float64)

    j = 0
    for i in range(store_start, 4*SPEC_END, 4):
            ydata[j]  = p_data[i][ilow:ihigh+1].real
            yfit[j]   = p_data[i+1][ilow:ihigh+1].real
            yres[j]   = p_data[i+2][ilow:ihigh+1].real
            yerror[j] = p_data[i+3][ilow:ihigh+1].real
            yflag[j]  = p_flag[i][ilow:ihigh+1]
            j += 1

    del(uvi)

    return ychannel, ydata, yflag, yres, yfit, yerror

def select_freq_1d(x1, y1, w1, flag1, low, high):
    prl = 0
    prh = 0
    for i in range(0, len(x1)):
            if x1[i]<=low:
                    i_low=i
                    prl = 1
            if x1[i]>=high:
                    i_high=i
                    prh = 1
                    break
    if prl==1 and prh==0:
            i_high = len(x1)
    if prl==0 and prh==1:
        i_low = 0	
        warnings.warn("Could not find the lower limit. Using the minimum")
    if prl==0 and prh==0:
        warnings.warn("Could not find the limits. Going to use the full array")
        i_low = 0
        i_high = len(x1)

    x1    = x1[i_low:i_high]
    y1    = y1[i_low:i_high]
    w1    = w1[i_low:i_high]
    flag1 = flag1[i_low:i_high]

    return x1, y1, w1, flag1, i_low, i_high


f_c = lambda f : np.ceil((f * 8192.0/250.0) + 1)
c_f = lambda c :        (((c-1) * 250.0/8192.0))

gaussian  = lambda mu, sigma, x : (1.0/(np.sqrt(2*np.pi*sigma**2)))*np.exp(-(x-mu)**2/(2*sigma**2))

gaussian_unity  = lambda amp, cen, wid, x : amp * np.exp(-(x-cen)**2/(2*wid**2))

def asym_gaussian(amp, mu, sigma1, sigma2, x):
    x_less  = (x <= mu)
    x_great = (x >  mu)
    ymodel = np.zeros(len(x))
    ymodel[x_less]  = gaussian_unity(amp, mu, sigma1, x[x_less]) 
    ymodel[x_great] = gaussian_unity(amp, mu, sigma2, x[x_great])
    return ymodel

nu_21_cm = 1420.405751

z_f = lambda z : nu_21_cm/(z+1)
f_z = lambda f : nu_21_cm/f - 1

del_z = lambda nu, del_nu: nu_21_cm * (1/nu**2) * del_nu
del_n = lambda nu, del_z : del_z * nu **2 / nu_21_cm

T21 = lambda nu, tau, A, nu0, w: A*(1-np.exp(-tau*np.exp(B(nu, nu0, w, tau))))/(1-np.exp(-tau))
B   = lambda nu, nu0, w, tau   : 4*(nu-nu0)**2/(w**2)*np.log(   -(1.0/tau) * np.log(  (1+np.exp(-tau))/2.0  )   )

def RMS_data(ee, res, correction=False, norder=False, imin=None, imax=None):
    if imin is None:
        imin=0
    if imax is None:
        imax = len(ee)
    if correction and not(norder):
        raise ValueError ("norder is needed to apply correction")
    if not(correction):
        rms = np.sqrt( np.sum((1.0/ee[imin:imax]**2) *(res)**2)/  np.sum(1.0/ee[imin:imax]**2)  )
    else:
        rms = np.sqrt(len(ee[imin:imax])/(len(ee[imin:imax])-(norder+1)))*\
	np.sqrt( np.sum((1.0/ee[imin:imax]**2) *(res)**2)/  np.sum(1.0/ee[imin:imax]**2)  )
    return rms

def estimate_noise(ee, xx, yy, norder=None, imin=None, imax=None, return_res=False):
    if imin is None:
        imin=0
    if imax is None:
        imax = len(ee)
    if norder is None:
        warnings.warn("You did not supply order. Defaulting to 15th order.")
        norder=15
    _,_,rr_res = poly_fit(xx[imin:imax], yy[imin:imax], norder, e = ee[imin:imax] )
    rms = np.sqrt(len(ee[imin:imax])/(len(ee[imin:imax])-(norder+1)))*\
	  RMS_data(ee[imin:imax], rr_res[imin:imax])
    if return_res:
        return rms, rr_res[imin:imax]
    return rms

def poly_fit(_x, _y, norder, e=None, domain='lin_lin'):
    if domain == 'log_lin':
        x = np.log10(_x)
        y = _y
    elif domain == 'log_log':
        x = np.log10(_x)
        y = np.log10(_y)
    else:
        x = _x
        y = _y
    if e is None:
        c_poly = np.polyfit(x, y, norder)
    else:
        c_poly = np.polyfit(x, y, norder, w=1/e)
    if domain== 'log_log':
        fit_poly = 10**np.polyval(c_poly, x)
    else:
        fit_poly = np.polyval(c_poly, x)
    
    res_poly = _y - fit_poly
    return c_poly, fit_poly, res_poly

import collections

def get_iterable(x):
    if isinstance(x, collections.Iterable):
        return x.ravel()
    else:
        return (x,)

def get_dict(*args, keys_array = None):
    param_dict = {}
    if keys_array is None:
        keys_array = ['p', 'x', 'imin', 'imax', 'fit', 'res', 'err', 'param_list' , 'cases']
    
    for val, key in enumerate(keys_array):
        param_dict[key] = {}
        param_dict[key] = args[val]
        
    return param_dict

class Fitter(object):
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
