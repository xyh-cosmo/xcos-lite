/*
    This likelihood function is actually the transformed prior, defined in NEW parameter space
    {eta1, eta2, ..., eta_{N-1}, E_N}.

    the mathematical formula is given by:
    p(eta1, eta2, ..., eta_{N-1}, E_N) = (E_N-1)^{N-1} \times \prod_{i=2}^{N-1} \eta_i^{i-1}

    The Normalization constant is dropped.
*/

#include "Cosmology.hpp"
#include "GSL_Math.hpp"
#include "Likelihoods.hpp"
#include <imcmc/imcmc.hpp>
#include <imcmc/parser++.hpp>

using namespace std;
using namespace imcmc;
using namespace imcmc::parser;
using namespace Cosmology;
using namespace Data;
using namespace Likelihoods;

namespace Likelihoods{

//  ================================================================================
//  NOTE:   This so-called likelihood funtion is actually NOT a likelihood function,
//          but the log-form of the transforming function from {E_i} to {eta_i}
//  ================================================================================

	double Likelihood_eta2E( imcmc_double& 	param,
	                         double& 		lndet,
	                         double& 		chisq,
	                         void* 			model,
	                         void* 			data,
                             istate&        state ){
	//	===================================================
        state.this_like_is_ok = true;
	    lndet = chisq = 0;
	//	===================================================

	//	CosmologyWorkspace should have been updated !!!
        CosmologyWorkspace  *cw  = static_cast<CosmologyWorkspace*>(model);
        int E_idx_max = cw->background->E_table_size;

        std::string var_E, var_eta;
        var_E = "E" + Read::IntToString(E_idx_max);
        chisq += (E_idx_max-1.0)*log(param[var_E]-1.0);

        for( int i=2; i<E_idx_max; ++i ){
            var_eta = "eta" + Read::IntToString(i);
            chisq   += (i-1.0)*log(param[var_eta]);
        }

        return -lndet - 0.5*chisq;
    }
}
