import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.size'] = 14
rcParams['font.family'] = 'serif'
rcParams['figure.autolayout'] = True
rcParams['lines.linewidth'] = 1.5
rcParams['lines.markersize'] = 6
rcParams['axes.titlesize'] = 14
rcParams['axes.labelsize'] = 14

af = np.loadtxt("error.txt")  #assumes af[:,0] = h, af[:,1] = error

fig = plt.figure()
plt.loglog(af[:,0], af[:,1], 'r-s', fillstyle='none', label='AF-3')


h_ref = af[:,0]
C = af[0,1] / h_ref[0]**3  
scale_factor = 0.7  # try 0.3, 0.5, 0.1 etc. to visually separate
plt.loglog(h_ref, scale_factor * C * h_ref**3, 'g--', label=r'Slope 3')

plt.xlabel('Mesh size')
plt.ylabel(r'$\mathrm{L}^1$ error')
plt.legend(loc='best')
plt.grid(True, linestyle='--')


from matplotlib.ticker import LogLocator, FormatStrFormatter

# Add more grid lines
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

from matplotlib.ticker import LogLocator

ax = plt.gca()

# Fewer minor ticks: only at [2, 5] between decades
ax.xaxis.set_major_locator(LogLocator(base=10.0))
ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=[2.0, 5.0], numticks=12))

ax.yaxis.set_major_locator(LogLocator(base=10.0))
ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=[2.0, 5.0], numticks=12))

plt.grid(True, which='both', linestyle='--', linewidth=0.6)




plt.savefig('order.pdf')
