import sys
import numpy as np
import matplotlib.pylab as plt

# data format:
# SN_name z   mu   dmu   mu_ref

def read_sne_mock( mock_filename ):
	fp = open(mock_filename,'r')
	lines = fp.readlines()
	fp.close()

	snls3 = []
	for line in lines:
		sn = line.split()
		temp = []
		temp.append(float(sn[1]))
		temp.append(float(sn[2]))
		temp.append(float(sn[3]))
		temp.append(float(sn[4]))
		snls3.append(temp)

	return np.array(snls3)

def gethist(outfig_name=None):
	if len(sys.argv) < 2:
		print 'usage: %s snls3_mock.txt'%(sys.argv[0])
		sys.exit(0)

	mock_num = len(sys.argv) - 1
	snls3 = read_sne_mock(sys.argv[1])
	z = snls3[:,0]
	dmu = (snls3[:,1]-snls3[:,3])/snls3[:,2]
	# dmu = (snls3[:,1]-snls3[:,3])

	nbin_all = 15
	nbin_1 = 15
	nbin_2 = 15
	z1 = 0.4
	z2 = 0.6

	# sort redshifts
	zsort = np.argsort(z)
	n = len(z)/3
	z1 = z[zsort[n-1]]
	z2 = z[zsort[2*n-1]]

	print 'z1 = %g, z2 = %g'%(z1,z2)
	
	plt.figure()
	# plt.errorbar(snls3[:,0],snls3[:,1],yerr=snls3[:,2])
	ID1 = (z <  z1 )
	ID2 = (z >= z2 )
	
	plt.hist(dmu, bins=nbin_all, label=r'ALL', alpha=0.7, rwidth=0.75)
	plt.hist(dmu[ID1], bins=nbin_1, label=r'$z < 0.3$', alpha=0.7, rwidth=0.75)
	plt.hist(dmu[ID2], bins=nbin_2, label=r'$z >= 0.8$', alpha=0.7, rwidth=0.75)

	plt.legend(loc='upper left')
	if outfig_name is None:
		plt.show()
	else:
		plt.savefig(outfig_name)



if __name__ == '__main__':
	if len(sys.argv) == 2:
		gethist()
	elif len(sys.argv) == 3:
		gethist(sys.argv[2])
