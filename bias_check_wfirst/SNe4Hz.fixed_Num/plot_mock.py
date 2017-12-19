import sys
from pylab import *
from scipy.interpolate import interp1d

if len(sys.argv) < 2:
	print 'usage: %s mock_1.txt mock_2.txt ...'%(sys.argv[0])
	sys.exit(0)

# load fiducial mu(z)
mu_fid = loadtxt('data/zmu_fid.txt')

fmu = interp1d(mu_fid[:,0],mu_fid[:,2])

sn = []

for i in range(1,len(sys.argv)):
	tmp = loadtxt(sys.argv[i])
	sn.append(tmp)


# plot
for i in range(len(sn)):
	errorbar(sn[i][:,0],
			 (sn[i][:,1]-fmu(sn[i][:,0]))/sn[i][:,2],
			 yerr=sn[i][:,2],fmt='.',label=sys.argv[i])

legend(loc='best',frameon=False)

show()
