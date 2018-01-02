import sys,os
import numpy as np
import matplotlib.pylab as plt
import scipy.stats as stats
from scipy.interpolate import interp1d

# load LCDM fiducial distance mudulus
z_dl_mu = np.loadtxt('data/zmu_fid.txt')
fun_mu = interp1d(z_dl_mu[:,0],z_dl_mu[:,2])

def mu_err(z,sigma_int=0.09,sigma_meas=0.08,sigma_lens=0.07,sigma_sys=0.02):
	sigma_tot = sigma_int**2 + sigma_meas**2 + (sigma_lens*z)**2 + (sigma_sys*(1.+z)/1.8)**2
	return sigma_tot**0.5
	
def get_wfirst_z(ntot=2725):
	bin_z = [ 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, \
			  0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65 ]	
	sn_num = [	69,	208,402,223,\
				327,136,136,136,\
				136,136,136,136,\
				136,136,136,136 ]

	sn_z = []
	for i in range(len(bin_z)):
		zmin = bin_z[i]-0.05
		zmax = bin_z[i]+0.05
		for j in range(sn_num[i]):
			z_rand = zmin + (zmax-zmin)*np.random.rand()
			sn_z.append(z_rand)
	
	return np.array(sn_z)

def get_mu_err(z):
	err = []
	for i in range(len(z)):
		err.append(mu_err(z[i]))

	return np.array(err)

# set the random number seed
# np.random.seed(1234567890)

num_of_mocks = 1000
P_VALUE = 0.

# test the existence of mock_WFIRST_1000
if os.path.isdir('mock_WFIRST_1000'):
	print 'directory: mock_WFIRST_1000 alreay exist!'
else:
	os.system('mkdir mock_WFIRST_1000')

cnt = 0
while cnt < num_of_mocks:
	z = get_wfirst_z()
	mu = fun_mu(z)
	mu_std = get_mu_err(z)
	dmu = np.zeros(mu.shape)
	for i in range(len(z)):
		dmu[i] = np.random.normal()*mu_std[i]
	
	ks_stats, pvalue = stats.kstest(dmu/mu_std,cdf='norm')
	if pvalue > P_VALUE:
		print 'yeap! got a good mock sample! pvalue = %g'%(pvalue)
		fname = 'mock_WFIRST_1000_B/WFIRST_SN_'+str(cnt+1)+'.txt'
		fp = open(fname,'w')
		for i in range(len(mu)):
			print >> fp, '%8.6f %10.8f %10.8f %10.8f'%(z[i],mu[i]+dmu[i],mu_std[i],mu[i])
		fp.close()
		cnt += 1
