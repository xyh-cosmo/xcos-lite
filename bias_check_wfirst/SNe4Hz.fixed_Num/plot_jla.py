from pylab import *

jla = loadtxt('jlaxxx.txt')

fig = figure(figsize=(6,6))

ax = fig.add_axes([0.1,0.3,0.8,0.6])
ax.errorbar(jla[:,0],jla[:,1],yerr=jla[:,2],fmt='.',label=r'JLA mock SNe sample',alpha=0.5)
ID=argsort(jla[:,0])
ax.plot(jla[ID,0],jla[ID,3],'r--',label=r'fiducial')
ax.set_ylabel(r'$\mu(z)$',fontsize=15)

ax = fig.add_axes([0.1,0.1,0.8,0.2])
ax.errorbar(jla[:,0],jla[:,1]-jla[:,3],yerr=jla[:,2],fmt='.',label=r'residuals',alpha=0.5)
ax.hlines(0,xmin=jla[:,0].min(),xmax=jla[:,0].max(),linestyles='--')
ax.set_xlabel(r'z',fontsize=15)
ax.set_ylabel(r'$\Delta \mu(z)$',fontsize=15)
show()