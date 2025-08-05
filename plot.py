#from mpl import *
import numpy as np
import matplotlib.pyplot as plt



sol = np.loadtxt("sol.txt")
sol_int   = np.loadtxt("sol_int.txt")
ex = np.loadtxt("exact.txt")
#--------------------------density--------------------

fig = plt.figure()
plt.plot(sol[:,0][::2],sol[:,1][::2],'s',fillstyle='none',label='AF-cell average',c='b')
plt.plot(sol_int[:,0][::2], sol_int[:,1][::2],'ro',fillstyle='none',label='AF-cell interface')
plt.plot(ex[:,0],ex[:,1],'-',fillstyle='none',label='Exact',c='k')
#plt.xlim(-1,1)
plt.xlabel('x')
plt.ylabel('u')
#sfig.legend(loc=15)
plt.legend();
#plt.legend(bbox_to_anchor=(1,0), loc="lower right")
#plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
               # mode="expand", borderaxespad=0, ncol=2)
                
plt.grid(True, linestyle = '--')
#plt.axis('equal')
#plt.yticks(fontsize=12)
plt.savefig('u.pdf')

#plt.xlim(0.2,0.75)
#plt.savefig('uzoomed.pdf')