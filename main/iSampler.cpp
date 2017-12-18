#include <iostream>
#include <fstream>
#include <iomanip>

#include "Cosmology.hpp"
#include "Likelihoods.hpp"
#include <imcmc/ensemble.hpp>
#include <imcmc/parser++.hpp>

using namespace std;
using namespace imcmc;
using namespace imcmc::parser;
using namespace Cosmology;
using namespace Data;
using namespace Likelihoods;

bool stop_on_gsl_err = false;

int main(int argc, char *argv[])
{
    MPI::Init(argc,argv);

    if( argc < 2 ){
        if( MPI::COMM_WORLD.Get_rank() == 0 ){
            cout << "usage :" << argv[0] << " *ini" << endl; 
        }
        MPI::Finalize();
        exit(0);
    }


    ensemble_workspace ew; // imcmc::ensemble sampler workspace

	Likelihood_Models	model;
	Likelihood_Data		data;

	model.Init(argv[1]);
	data.Init(argv[1]);

    // stop on GSL error?
    stop_on_gsl_err = Read::Read_Bool_from_File(argv[1],"stop_on_gsl_err");
    if( stop_on_gsl_err != true )
        gsl_set_error_handler_off();

//  TODO: add a new API inside Likelihood_Models, so one can add nuisance parameters in a simpler way such as:
	if( data.use_sne == true ){
		// MPI_cout("copying SNe Ia Light-Curve parameters into Params ...", MPI::COMM_WORLD.Get_rank());
		data.sne->Params.CopyInto( model.Params );
	}

    if( model.has_derived_params )
	    ew.add_likelihood( Likelihood_All, model.Params.MCMC_Params, model.derived_params, &model, &data );
    else
        ew.add_likelihood( Likelihood_All, model.Params.MCMC_Params, &model, &data );

    ew.init(argv[1]);

    ew.do_sampling();

    MPI::Finalize();
}
