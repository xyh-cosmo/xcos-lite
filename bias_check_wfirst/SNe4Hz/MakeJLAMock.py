import sys as sys
import numpy as np
from scipy.linalg import cholesky
import scipy.stats as stats
from scipy.interpolate import interp1d

# load LCDM fiducial distance mudulus
z_dl_mu = np.loadtxt('data/zmu_fid.txt')
fun_mu = interp1d(z_dl_mu[:,0],z_dl_mu[:,2])

# load JLA redshifts
z_jla = np.loadtxt('data/jla_z.txt')

# load covariance matrix
covmat = np.loadtxt('data/JLA_cov.txt')
var_jla= covmat.diagonal()**0.5


A = cholesky(covmat,lower=True)

def gen_err(covmat=None,use_full_cov=True):
	if covmat is None:
		print('covmat is not set, exit!')
		sys.exit(0)
	
	ndim = covmat.shape[0]
	mean = np.zeros(covmat.shape[0])
	if use_full_cov is False:
		covtemp = covmat.copy()
		covmat = np.zeros(covtemp.shape)
		for i in range(ndim):
			covmat[i,i] = covtemp[i,i]
	
	err = np.random.multivariate_normal(mean,covmat)
	return err

# set the random number seed
#np.random.seed(1234567890)

num_of_mocks = 1000
P_VALUE = 0.

cnt = 0
while cnt < num_of_mocks:
	mu = fun_mu(z_jla)
#	mu_err = gen_err(covmat,use_full_cov=True)
	mu_err = gen_err(covmat,use_full_cov=False)
	# r = np.random.randn(len(z_jla))
	# mu_err = np.matmul(A,r)
#	ks_stats, pvalue = stats.kstest(mu_err/var_jla,cdf='norm')
	# if pvalue >= P_VALUE:
#	print 'yeap! got a good mock sample! pvalue = %g'%(pvalue)
	fname = 'mock_JLA_1000/MOCK_JLA_'+str(cnt+1)+'.txt'
	fp = open(fname,'w')
	for i in range(len(mu)):
		print >> fp, '%s %8.6f %10.8f %10.8f %10.8f'%('mock_SN',z_jla[i],mu[i]+mu_err[i],var_jla[i],mu[i])
	fp.close()
	print 'saved %s'%(fname)
	cnt += 1
