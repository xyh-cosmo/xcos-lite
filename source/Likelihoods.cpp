/*
 *	This is a unified interface to all the likelihood functions, because Cosmological Model
 *	is shared by more than one likelihood functions, thus put them together is more convenient.
 */

#include "Cosmology.hpp"
#include "Likelihoods.hpp"
#include <imcmc/imcmc.hpp>
#include <imcmc/ensemble.hpp>
#include <imcmc/parser++.hpp>

using namespace std;
using namespace imcmc;
using namespace imcmc::parser;
using namespace Cosmology;
using namespace Data;
using namespace Likelihoods;

namespace Likelihoods{

//	===========================
//	Likelihood Models
//	===========================

	Likelihood_Models::Likelihood_Models(){
        my_rank = MPI::COMM_WORLD.Get_rank();
		use_cosmo_workspace 	= false;
        use_eta2E               = false;
		has_derived_params		= false;
        derive_hubble_rate      = false;
		derive_expansion_rate	= false;
        cw = NULL;
	}

	Likelihood_Models::~Likelihood_Models(){
		if ( cw != NULL ){
			delete cw;
			cw = NULL;
		}
	}

	bool Likelihood_Models::Init( std::string paramfile ){

        MPI_cout("Initializing Cosmology Workspace ...", my_rank);

		cw	= new CosmologyWorkspace;

		if( cw->Init(paramfile) == false )
			return false;

		use_cosmo_workspace	= true;

		cw->Params.CopyInto(this->Params);

	//	========================================================================
	//	add derived parameters
	//	========================================================================
		has_derived_params	= Read::Read_Bool_from_File(paramfile, "has_derived_params");

		if( has_derived_params ){
			derive_expansion_rate 	= Read::Read_Bool_from_File(paramfile, "derive_expansion_rate");
			derive_hubble_rate		= Read::Read_Bool_from_File(paramfile, "derive_hubble_rate");
            use_eta2E               = Read::Read_Bool_from_File(paramfile, "use_eta2E");
		}
        else{
            MPI_cout("*** no derived parameters ***", my_rank);
        }

	//	------------------------------------------------------------------------
	//	read zi from *.ini file, keyword is 'derive_E_H_at_z'
	//	NOTE: the redshifts MUST be given in increasing order!!!
	//	------------------------------------------------------------------------
		if( derive_expansion_rate || derive_hubble_rate ){

			int nz_E=-1;
			int nz_H=-1;
			double *z_E = NULL;
			double *z_H = NULL;

			if( Read::Has_Value(paramfile, "derive_E_at_z", "double") ){
				nz_E = Read::Num_of_Value_for_Key(paramfile, "derive_E_at_z");
				z_E = new double[nz_E];
				Read::Read_Array_of_Double_from_File(paramfile, "derive_E_at_z", z_E, nz_E);
			}

			if( Read::Has_Value(paramfile, "derive_H_at_z", "double") ){
				nz_H = Read::Num_of_Value_for_Key(paramfile, "derive_H_at_z");
				z_H = new double[nz_H];
				Read::Read_Array_of_Double_from_File(paramfile, "derive_H_at_z", z_H, nz_H);
			}

			if( derive_expansion_rate && (nz_E > 0) ){
				for( int i=1; i<=nz_E; ++i ){
					string pname = "E_"+Read::IntToString(i);
					expansion_rate.push_back(pname);
					z_of_expansion_rates[pname] = z_E[i-1];
					derived_params.push_back(pname);
				}

				delete[] z_E;
			}

			if( derive_hubble_rate && (nz_H > 0)  ){
				for( int i=1; i<=nz_H; ++i ){
					string pname = "H_"+Read::IntToString(i);
					hubble_rate.push_back(pname);
					z_of_hubble_rates[pname] = z_H[i-1];
					derived_params.push_back(pname);
				}

				delete[] z_H;
			}

		}	//

        if( use_eta2E ){
        //  add E_1, E_2, ..., E_{N-1} into derived_params.
        //  NOTE: derived {E_i} CAN NOT have the same name as {E_i} used inside in class background
		//	Note (added @Apr-5-2017): no need to specify "derive_E_at_z" in the input *.ini parameter file.
            std::string var_E, var_eta;
            for( int i=1; i<=cw->background->E_table_size-1; ++i ){
                std::string idx = Read::IntToString(i);
                var_E   = "derived_E" + idx;
                var_eta = "eta" + idx;
                derived_params.push_back(var_E);
                Params.AddParam(var_eta, true);
            }
        }

		return true;
	}

//	===========================
//	Likelihood Data sets
//	===========================

	Likelihood_Data::Likelihood_Data(){
        my_rank = MPI::COMM_WORLD.Get_rank();
		sne 	= new SNeIa;
        use_cosmo_data = false;
	}

	void Likelihood_Data::Init( std::string paramfile ){
		use_sne 	= Read::Read_Bool_from_File(paramfile, "use_sne");
		if( use_sne ){
			sne->Init(paramfile);
			MPI_cout("using SNeIa", my_rank);
		}
        use_cosmo_data = use_sne;
		MPI_cout("finishing setting likelihood data sets", my_rank);
	}

	Likelihood_Data::~Likelihood_Data(){
		if( use_sne ){
			delete sne;
			sne = NULL;
		}
	}

//	===========================
//	Likelihood functions
//	===========================

	double Likelihood_All(  imcmc_double&   param,
                            double&         lndet,
                            double&         chisq,
                            void*           model,
                            void*           data,
                            istate&         state ){

        state.this_like_is_ok=true;		//	this passed to imcmc::ensemble

		Likelihood_Models*	LM = static_cast<Likelihood_Models*>(model);
		Likelihood_Data*	LD = static_cast<Likelihood_Data*>(data);

    //  Need to solve {E_i} from {eta_i} and E_N before updating Models
        if( LM->use_eta2E ){

            double E_N, Exx;
            int E_index_max = LM->cw->background->E_table_size;
            std::string var_E, var_derE, var_eta;

            var_E   = "E" + Read::IntToString(E_index_max);
            E_N     = param[var_E];
            Exx     = E_N;

            for( int i=E_index_max-1; i>=1; --i ){
                std::string idx = Read::IntToString(i);
                var_E           = "E" + idx;
                var_derE        = "derived_E" + idx;
                var_eta         = "eta" + idx;
                param[var_derE] = param[var_eta] * (Exx-1.0) + 1.0; //  update derived_Ei
                param[var_E]    = param[var_derE];                  //  update E_i
                Exx             = param[var_derE];                  //  update Exx for next round.
            }
        }

        if( LD->use_cosmo_data == true ){
	        LM->cw->background->Params.UpdateParams(param);			//	Update cosmological parameters

			if( LM->cw->background->Update() == false ){ //	update cosmological background quantities.
				state.this_like_is_ok = false;
				chisq = 1e50;
				return -0.5*chisq;
			}
        }

		if( LM->has_derived_params ){

			imcmc_vector_string_iterator it;

			if( LM->derive_expansion_rate ){
				it = LM->expansion_rate.begin();
				double H0 = LM->cw->background->Hubble_z(0.0);
				while( it != LM->expansion_rate.end() ){
					param[*it] = LM->cw->background->Hubble_z(LM->z_of_expansion_rates[*it]) / H0;
					++it;
				}
			}

			if( LM->derive_hubble_rate ){
				it = LM->hubble_rate.begin();
				while( it != LM->hubble_rate.end() ){
					param[*it] = LM->cw->background->Hubble_z(LM->z_of_hubble_rates[*it]);
					++it;
				}
			}
		}

		double lndet_temp = 0;
		double chisq_temp = 0;

		lndet = chisq = 0;

		if( LD->use_sne ){
			Likelihood_SNE( param,
                            lndet_temp,
                            chisq_temp,
                            LM->cw,
                            LD->sne,
                            state );

			lndet += lndet_temp;
			chisq += chisq_temp;
		}

        if( LM->use_eta2E ){
			Likelihood_eta2E( param,
                              lndet_temp,
                              chisq_temp,
                              LM->cw,
                              NULL,
                              state );

			lndet += lndet_temp;
			chisq += chisq_temp;
        }

        return -lndet - 0.5*chisq;
	}

}
