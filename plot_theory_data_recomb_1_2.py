import numpy as np
import matplotlib.pyplot as plt

f2 = 'total_spec_new.dat'

#x1,y1 = np.loadtxt(f1,unpack = True)
x2,y2 = np.loadtxt(f2,unpack = True)
plt.figure()
#plt.plot(np.log10(x1*1.0e9),y1,label='old')
plt.plot(np.log10(x2*1.0e9),y2,label='new')
plt.xlabel('log10[Freq GHz]')
plt.ylabel('Intensity')
plt.legend(loc='upper left')
plt.grid()
plt.show()
