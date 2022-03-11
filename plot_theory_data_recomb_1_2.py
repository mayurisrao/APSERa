import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import matplotlib
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
import matplotlib.colors as mcolors
matplotlib.rc('font', **font)

f2 = 'total_spec_new.dat'

#x1,y1 = np.loadtxt(f1,unpack = True)
x2,y2 = np.loadtxt(f2,unpack = True)
# plt.figure()
# #plt.plot(np.log10(x1*1.0e9),y1,label='old')
# plt.plot(np.log10(x2*1.0e9),y2,label='new')
# plt.xlabel('log10[Freq GHz]')
# plt.ylabel('Intensity')
# plt.legend(loc='upper left')
# plt.grid()
# plt.show()
fig, ax = plt.subplots()
ax.plot(x2,y2,label='new')
ax.set_xlim([4.0,16])
ax.set_ylim([0.5e-27,1.5e-27])
line1 = 4.6
line2 = 8.5
line3 = 8.3 
line4 = 15.3
ax.axvspan(line1, line2, alpha=.35, color='orange')
ax.axvspan(line3, line4, alpha=.35, color='silver')
ax.text(12, 10, 'Band 5a', fontsize=25)
ax.set_xlabel(r'$\nu$ [GHz]')
ax.set_ylabel(r'$I_{\nu}$ [Jm$^{-2}$s$^{-1}$Hz$^{-1}$sr$^{-1}$]')
ax.grid()
plt.show()
