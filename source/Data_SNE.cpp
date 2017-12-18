//  =======
//   SNeIa
//  =======

#include <armadillo>
#include <imcmc/imcmc.hpp>
#include <imcmc/parser++.hpp>
#include "Cosmology.hpp"
#include "SNE.hpp"

using namespace std;
using namespace imcmc;
using namespace imcmc::parser;
using namespace arma;

namespace Data{

	SNeIa::SNeIa(){
		initialized = false;
		sne_union	= NULL;
		sne_snls	= NULL;
		sne_jla		= NULL;
		sne_lsst    = NULL;
		sne_wfirst_mock		= NULL;
		sne_snls_jla_mock 	= NULL;
		mock_marg_MB 		= false;
		my_rank		= MPI::COMM_WORLD.Get_rank();
	}

	SNeIa::~SNeIa(){
		if( (sne_data_name == _MOCK_) && initialized ){
			if( sne_wfirst_mock != NULL )
				delete sne_wfirst_mock;
			if( sne_snls_jla_mock != NULL )
				delete sne_snls_jla_mock;
		}

		if( ((sne_data_name == _UNION_) && sne_union != NULL) && initialized )
			delete sne_union;

		if( ((sne_data_name == _SNLS_) && sne_snls != NULL) && initialized )
			delete sne_snls;

		if( ((sne_data_name == _JLA_) && sne_jla != NULL) && initialized )
			delete sne_jla;

		if( ((sne_data_name == _LSST_) && sne_lsst != NULL) && initialized )
			delete sne_lsst;
	}

	void SNeIa::Init( std::string paramfile ){

		if( Read::Read_Bool_from_File(paramfile, "use_mock_sne") ){
		    MPI_cout("using mock SNe sample", my_rank);
		    sne_data_name   = _MOCK_;

			if( Read::Read_Bool_from_File(paramfile,"use_wfirst_mock") ){
				MPI_cout("using WFIRST SNe sample",my_rank);
				sne_wfirst_mock = new WFIRST_Mock;
			    sne_wfirst_mock->ReadData(paramfile);
			}

			if( Read::Read_Bool_from_File(paramfile,"use_snls_jla_mock") ){
				MPI_cout("using SNLS or JLA mock SNe sample",my_rank);
				sne_snls_jla_mock = new SNLS_JLA_mock;
				sne_snls_jla_mock->ReadData(paramfile);
			}

			mock_marg_MB = Read::Read_Bool_from_File(paramfile,"mock_marg_MB");

			if( mock_marg_MB == false ){
				if( Read::Has_Key_in_File(paramfile,"MB") ){
					Params.AddParam("MB",true);
					if( Read::Has_Key_in_File(paramfile,"MB_std") )
						MB_std = Read::Read_Double_from_File(paramfile,"MB_std");
					else
						throw runtime_error("==> can not find MB_std in: "+paramfile);
				}
				else{
					throw runtime_error("==> can not find MB in: "+paramfile);
				}

				mock_MB_use_Gaussian_Prior = Read::Read_Bool_from_File(paramfile,"mock_MB_use_Gaussian_Prior");
			}

			initialized     = true;
		}
		else if( Read::Read_Bool_from_File(paramfile, "use_union") ){
		    MPI_cout("using Union2.1 SNe sample", my_rank);
		    sne_data_name   = _UNION_;
		    sne_union       = new UNION;
		    sne_union->ReadData(paramfile);
		    //  no nuisance parameter is used

		    initialized     = true;
		}
		else if( Read::Read_Bool_from_File(paramfile, "use_snls") ){
		    MPI_cout("using SNLS SNe sample", my_rank);
		    sne_data_name   = _SNLS_;
		    sne_snls        = new SNLS;
		    sne_snls->ReadData(paramfile);

		    //  nuisance parameters: alpha, beta, scriptMB
		    Params.AddParam("alpha_snls", true);
		    Params.AddParam("beta_snls", true);
		    Params.AddParam("scriptMB_snls", true);

		    initialized     = true;
		}
		else if( Read::Read_Bool_from_File(paramfile, "use_jla") ){
		    MPI_cout("using JLA SNe sample", my_rank);
		    sne_data_name   = _JLA_;
		    sne_jla         = new JLA;
		    sne_jla->ReadData(paramfile);

		    //  nuisance parameters: alpha, beta, scriptMB
		    Params.AddParam("alpha_jla", true);
		    Params.AddParam("beta_jla", true);
		    Params.AddParam("MB_jla", true);
		    Params.AddParam("DeltaMB_jla", true);

		    initialized     = true;
		}
		else if( Read::Read_Bool_from_File(paramfile, "use_lsst_sn") ){
			MPI_cout("using LSST photoz SNe mock sample", my_rank);
			sne_data_name	= _LSST_;
			sne_lsst		= new SNeIa_LSST;
			sne_lsst->ReadData(paramfile);
			Params.AddParam("MB",true);
			Params.AddParam("MB_a",true);
			Params.AddParam("MB_b",true);
		}
		else{
		    string err = "*** SNeIa::Init() ==> no SNeIa data available!";
		    throw runtime_error(err);
		}

	}

}
