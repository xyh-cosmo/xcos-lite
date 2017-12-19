import sys,os
import numpy as np
from scipy.stats import kstest

def read_mock_sne(mock):
    if os.path.isfile(mock):
        sn = []
        fp = open(mock)
        lines = fp.readlines()
        for line in lines:
            text = line.split()
            sn_i = []
            sn_i.append(np.float64(text[1]))
            sn_i.append(np.float64(text[2]))
            sn_i.append(np.float64(text[3]))
            sn_i.append(np.float64(text[4]))
            sn.append(sn_i)
        
        fp.close()
        return np.array(sn)
    else:
        print('can not open: %s'%(mock))
        sys.exit(0)



if __name__ == '__main__':
	if len(sys.argv) < 2:
		print 'useage: '
		print '\t%s mock_sne.txt'%(sys.argv[0])

	for i in range(len(sys.argv)-1):
		fname = sys.argv[i+1]
		sn = read_mock_sne(fname)
		s,p= kstest((sn[:,1]-sn[:,3])/sn[:,2], cdf='norm')
		print 'KStest result of: %s, S = %8.6f, P = %8.6f'%(fname,s,p)