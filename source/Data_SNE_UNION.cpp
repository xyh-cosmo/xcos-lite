//  ========================================================
//  ==================  Union2.1    ========================
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

	void UNION::ReadData( string paramfile ){

#if defined(_DEBUG_UNION_)
		cout << "--->   reading union2.1 data\n";
#endif

		std::string dataset = Read::Read_String_from_File( paramfile, "union2.1_dataset" );
		std::string data 	= Read::Read_String_from_File( dataset, "data_file" );
		int size 			= Read::Read_Int_from_File( dataset, "data_size" );
		sne_num = size;

		z   = arma::zeros(size, 1);
		mu  = arma::zeros(size, 1);
		dmu = arma::zeros(size, 1);
		P   = arma::zeros(size, 1);
		icov= arma::zeros(size, size);

		systematic = Read::Read_Bool_from_File( paramfile, "union2.1_with_sys" );
		string covmat;
		if( systematic == true )
			covmat = Read::Read_String_from_File( dataset, "cov_sys" );
		else if( systematic == false )
			covmat = Read::Read_String_from_File( dataset, "cov_nosys" );

	//	now read in data
		std::ifstream infile_data(data.c_str());
		std::string line, sn_name;

		int i=0;
		while( std::getline(infile_data, line) ){
		    if( !Read::Is_Commented(line) && !Read::Is_Empty(line) ){
				std::stringstream stream(line);
				stream >> sn_name >> z(i) >> mu(i) >> dmu(i) >> P(i);
				++i;
		    }
		}

		infile_data.close();

		if( i != size )
		    throw runtime_error("*** UNION::ReadData() ==> wrong number of UNION2.1 distance moduli data !");

	//  now let's read in covariance matrix

		std::ifstream infile_cov(covmat.c_str());	// size is 580*580

		if( !infile_cov )
		    throw runtime_error("*** UNION::ReadData() ==> cannot open: " + covmat);
		else{
			i = 0;
			while( std::getline(infile_cov, line) ){
		        if( !Read::Is_Commented(line) && !Read::Is_Empty(line) ){
					std::stringstream stream(line);
					double covx;
					for(int j=0; j<size; ++j){
						stream >> covx;
		                icov(i,j) = covx;
					}
					++i;
		        }
			}
		}

		infile_cov.close();

		if( i != size )
		    throw runtime_error("*** UNION::ReadData() ==> wrong number of UNION2.1 distance moduli cov data !");
	}


}
