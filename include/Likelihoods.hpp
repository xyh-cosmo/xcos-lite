#include "Cosmology.hpp"
#include "SNE.hpp"
#include "Prior_eta2E.hpp"

#ifndef __ICOSMO_LIKELIHOODS__
#define __ICOSMO_LIKELIHOODS__

using namespace Cosmology;
using namespace Data;
using namespace Likelihoods;

namespace Likelihoods{

/*
	struct Likelihood_Models contains all models needed by Likelihood_All(), which includes actually
	more than one likelihoods, each corresponding to a specific data set.
*/

	struct Likelihood_Models{

        int     my_rank;
		bool	use_cosmo_workspace;
	
        bool    use_eta2E;   //  @2016-2-25, add support for restricted ordering sampling of {E_i}
	//	this is used to output derived cosmic expansion rate, only used when {w(zi)} is used.
		bool	has_derived_params;		//	default: false.
		bool	derive_expansion_rate;	
		bool	derive_hubble_rate;		//	similar to derive_hubble_rate
		imcmc::imcmc_vector_string derived_params;

		ParamVector Params;		//	all model parameters

		CosmologyWorkspace	*cw;		//	first comes the cosmology workspace

		bool Init( std::string paramfile );
		bool Update( imcmc::imcmc_double param );

		Likelihood_Models();
		~Likelihood_Models();

		imcmc::imcmc_vector_string expansion_rate;		// E1, E2, ...
		imcmc::imcmc_vector_string hubble_rate;		     // H1, H2, ...
		imcmc::imcmc_double z_of_expansion_rates;		// z1, z2, ...
		imcmc::imcmc_double z_of_hubble_rates;			// z1, z2, ...
	};

	struct Likelihood_Data{
        int     my_rank;
		bool	use_sne;
        bool    use_cosmo_data;
		SNeIa	*sne;

		void Init( std::string paramfile );

		Likelihood_Data();
		~Likelihood_Data();
	};

//  ============================================================================
//  Likelihood_All returns the sum of likelihoods from all used data sets.
	double Likelihood_All(  imcmc::imcmc_double     &param,
                            double                  &lndet,
                            double                  &chisq,
                            void                    *model,
                            void                    *data,
                            imcmc::istate&          state );
};


#endif	//	__ICOSMO_LIKELIHOODS__
