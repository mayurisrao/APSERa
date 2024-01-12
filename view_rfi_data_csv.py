import numpy as np
import matplotlib.pyplot as plt
from numpy import genfromtxt
import os
import matplotlib
import pandas as pd
from scipy import stats
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams["font.size"] = "18"
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib import colors
import sys
from pathlib import Path

print("Enter filename without extension")
filename = sys.argv[1]+".csv"
fout1 = filename + "_waterfall.png"
fout2 = filename +'_med.png'
rfi_df = pd.read_csv(filename, header=None)
rfi_data = np.array(rfi_df)

plt.figure(figsize=(9,15))
plt.imshow(rfi_data, aspect='auto', vmax=-60, vmin=-90, cmap='inferno', extent=(0.5,4.0,0,len(rfi_data)),interpolation='none')
plt.grid(c='w')
plt.colorbar(label='dBm')
plt.title(filename)
plt.xlabel('Frequency(GHz)')
plt.ylabel('time(seconds x 2)')
plt.savefig(fout1)
plt.show()

data_median = np.median(rfi_data[:,2:],axis=0)
plt.figure(figsize=(9,15))
plt.rcParams.update({'font.size': 10})
plt.plot(np.linspace(0.5,4.0,len(data_median)), data_median , linewidth=0.5)
plt.grid()
plt.xlabel('frequncy(GHz)')
plt.ylabel('intensity (dBm)')
plt.title(filename)
plt.legend(loc='upper right')
plt.savefig(fout2)
plt.show()
