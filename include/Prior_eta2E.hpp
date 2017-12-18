#include "Cosmology.hpp"

#ifndef __ETA2E__
#define __ETA2E__

using namespace Cosmology;
using namespace Data;
using namespace Simulation;
using namespace Likelihoods;
using namespace imcmc;

namespace Likelihoods{

	double Likelihood_eta2E( imcmc_double& 	param,
	                         double& 		lndet,
	                         double& 		chisq,
	                         void* 			model,
	                         void* 			data,
                             istate&        state );

}


#endif	//	__ETA2E__
