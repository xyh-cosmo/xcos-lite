//  ========================================================
//  ==================  Mock samples  ======================
//  ========================================================
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

	WFIRST_Mock::WFIRST_Mock(){
		my_rank			= MPI::COMM_WORLD.Get_rank();
		initialized 	= false;
	}

	WFIRST_Mock::~WFIRST_Mock(){
		if( initialized ){
			delete[] sne_z;
			delete[] sne_mu;
			delete[] sne_dmu;
		}
	}

	void WFIRST_Mock::ReadData( string paramfile ){

		vector<double> z, mu, dmu;

		string fname;
		fname = Read::Read_String_from_File(paramfile, "mock_wfirst_datafile");
		ifstream infile(fname.c_str());

	    bool has_err = true;	//	deaut include errors
	    has_err = Read::Read_Bool_from_File(paramfile, "mock_wfirst_has_error");

		if( !infile.good() )
		    throw runtime_error("*** void SNe_Mock::ReadData() ==> failed to open: " + fname);
		else{
		    string line;

		    int count=0;

		    MPI_cout("reading mock sne sample: " + fname, my_rank);

		    while( getline(infile, line) ){

		        if( !Read::Is_Commented(line) && !Read::Is_Empty(line) ){
		            stringstream stream(line);
		            double zi, mui, dmui, mui_theory;
		            stream >> zi >> mui >> dmui >> mui_theory;  //  mui_thoery is NOT used, it is kept for debug.
		            z.push_back(zi);

		            if( has_err )
						mu.push_back(mui);
		            else
						mu.push_back(mui_theory);

		            dmu.push_back(dmui);
		            ++count;
		        }
		    }

		    this->sne_num   = count;
		    this->sne_z     = new double[count];
		    this->sne_mu    = new double[count];
		    this->sne_dmu   = new double[count];

		    for(int i=0; i<this->sne_num; ++i){
		        sne_z[i]    = z[i];
		        sne_mu[i]   = mu[i];
		        sne_dmu[i]  = dmu[i];
		    }

		    MPI_cout("finished reading: " + Read::IntToString(count) + " mock SNeIa are loaded", my_rank);
		}

		infile.close();
	}

	SNLS_JLA_mock::SNLS_JLA_mock(){
		my_rank		= MPI::COMM_WORLD.Get_rank();
		initialized	= false;
		use_full_covmat = true;
	}

	SNLS_JLA_mock::~SNLS_JLA_mock(){
		// do nothing ...
	}

	void SNLS_JLA_mock::ReadData( string paramfile ){
		vector<double> vec_z, vec_mu, vec_dmu, vec_mu_no_err;

		if( !Read::Has_Key_in_File(paramfile, "mock_snls_jla_datafile") ){
			throw runtime_error("*** cannot find key: mock_snls_jla_datafile in file:"+paramfile+", stop!");
		}

		string fname;
		fname = Read::Read_String_from_File(paramfile,"mock_snls_jla_datafile");
		ifstream infile(fname.c_str());

		bool has_err = true;
		has_err = Read::Read_Bool_from_File(paramfile,"mock_snls_jla_has_err");

		if( Read::Has_Key_in_File(paramfile,"mock_snls_jla_use_full_covmat") ){
			use_full_covmat = Read::Read_Bool_from_File(paramfile,"mock_snls_jla_use_full_covmat");
		}

		if( !infile.good() ){
			throw runtime_error("*** failed to open: "+fname);
		}
		else{
			string line;
			int count = 0;
			while( getline(infile,line) ){
				if( !Read::Is_Commented(line) && !Read::Is_Empty(line) ){
					stringstream stream(line);
					string name;
					double zi, mui, dmui, muixx;
					stream >> name >> zi >> mui >> dmui >> muixx;

					vec_z.push_back(zi);
					vec_dmu.push_back(dmui);
					if( has_err )
						vec_mu.push_back(mui);
					else
						vec_mu.push_back(muixx);

					++count;
				}
			}

			sne_num = count;

			z 	= arma::zeros(count);
			mu	= arma::zeros(count);
			mB	= arma::zeros(count);
			dmu	= arma::zeros(count);

			for( int i=0; i<count; ++i ){
				z(i)	= vec_z[i];
				mu(i)	= vec_mu[i];
				dmu(i)	= vec_dmu[i];
			}

			mB = mu + (-19.3);
		}

		// read covariance
		if( !Read::Has_Key_in_File(paramfile,"mock_snls_jla_covmat") ){
			throw runtime_error("*** cannot find covariance matrix file in: "+paramfile);
		}

		string covmat_file = Read::Read_String_from_File(paramfile,"mock_snls_jla_covmat");
		cov.load(covmat_file,arma::raw_ascii);

		if( use_full_covmat == false ){
			cout << "==> use_full_covmat == false, so only the diagonal of the covariance matrix will be used.\n";
			arma::mat cov_temp = cov;
			cov = arma::diagmat(cov_temp);
		}
		
		icov = cov.i();
	}

}
