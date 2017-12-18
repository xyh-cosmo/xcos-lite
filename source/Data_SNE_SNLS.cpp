//  ========================================================
//  ==================   SNLS   ============================
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

	void SNLS::ReadData( string paramfile ){
		string dataset      = Read::Read_String_from_File(paramfile, "snls_dataset");
		string data_file    = Read::Read_String_from_File(dataset, "data_file");

		int size            = Read::Read_Int_from_File(dataset, "size");
		sne_num = size;

	//  allocate memory
	//	1)	light-curve data
		zcmb        = arma::zeros(size,1);
		zhel        = arma::zeros(size,1);
		dz          = arma::zeros(size,1);
		mb          = arma::zeros(size,1);
		dmb         = arma::zeros(size,1);
		s           = arma::zeros(size,1);
		ds          = arma::zeros(size,1);
		c           = arma::zeros(size,1);
		dc          = arma::zeros(size,1);
		var3        = arma::zeros(size,1);
		dvar3       = arma::zeros(size,1);
		cov_m_s     = arma::zeros(size,1);
		cov_m_c     = arma::zeros(size,1);
		cov_s_c     = arma::zeros(size,1);

	//	2)	covariance matrix
		cov_mB_mB       = arma::zeros(size,size);
		cov_mB_alpha    = arma::zeros(size,size);
		cov_mB_beta     = arma::zeros(size,size);
		cov_alpha_alpha = arma::zeros(size,size);
		cov_alpha_beta  = arma::zeros(size,size);
		cov_beta_beta   = arma::zeros(size,size);

		cov_tot         = arma::zeros(size,size);
		icov_tot        = arma::zeros(size,size);

		pecz            = Read::Read_Double_from_File(dataset, "pecz");

		intrinsicdisp   = Read::Read_Double_from_File(dataset, "intrinsicdisp");

		use_four_disp	= Read::Read_Bool_from_File(dataset, "use_four_disp");
		intrinsicdisp0  = Read::Read_Double_from_File(dataset, "intrinsicdisp0");
		intrinsicdisp1  = Read::Read_Double_from_File(dataset, "intrinsicdisp1");
		intrinsicdisp2  = Read::Read_Double_from_File(dataset, "intrinsicdisp2");
		intrinsicdisp3  = Read::Read_Double_from_File(dataset, "intrinsicdisp3");

		twoscriptmfit   = Read::Read_Bool_from_File(dataset, "twoscriptmfit");
		scriptmcut      = Read::Read_Double_from_File(dataset, "scriptmcut");

		string mag_covmat_file;
		string mag_stretch_covmat_file;
		string mag_colour_covmat_file;
		string stretch_covmat_file;
		string stretch_colour_covmat_file;
		string colour_covmat_file;

		mag_covmat_file             = Read::Read_String_from_File(dataset, "mag_covmat_file");
		stretch_covmat_file         = Read::Read_String_from_File(dataset, "stretch_covmat_file");
		colour_covmat_file          = Read::Read_String_from_File(dataset, "colour_covmat_file");
		mag_stretch_covmat_file     = Read::Read_String_from_File(dataset, "mag_stretch_covmat_file");
		mag_colour_covmat_file      = Read::Read_String_from_File(dataset, "mag_colour_covmat_file");
		stretch_colour_covmat_file  = Read::Read_String_from_File(dataset, "stretch_colour_covmat_file");

	//	read data
		std::string line;
		int sn_num=0;

#ifdef _DEBUG_SNLS_
	std::cout << "---->	reading SNLS light-curve data\n";
#endif

		std::ifstream infile( data_file.c_str() );

		if( !infile )
		    throw runtime_error("*** SNLS::ReadData() ==> cannot open: " + data_file);

		std::getline(infile, line);
		while( std::getline(infile, line) ){
		    if( !Read::Is_Commented(line) && !Read::Is_Empty(line) ){
			    std::string sn_namex;
			    double zcmbx;
			    double zhelx;
			    double dzx;
			    double mbx;
			    double dmbx;
			    double sx;
			    double dsx;
			    double cx;
			    double dcx;
			    double var3x;
			    double dvar3x;
			    double cov_m_sx;
			    double cov_m_cx;
			    double cov_s_cx;
			    int setx;

				std::stringstream stream(line);
				stream 	>> sn_namex
						>> zcmbx >> zhelx >> dzx
						>> mbx >> dmbx
						>> sx >> dsx
						>> cx >> dcx
						>> var3x >> dvar3x
						>> cov_m_sx
						>> cov_m_cx
						>> cov_s_cx
						>> setx;

				sne_name.push_back(sn_namex);
				set.push_back(setx);

		        zcmb(sn_num)    = zcmbx;
		        zhel(sn_num)    = zhelx;
		        dz(sn_num)      = dzx;
		        mb(sn_num)      = mbx;
		        dmb(sn_num)     = dmbx;
		        s(sn_num)       = sx;
		        ds(sn_num)      = dsx;
		        c(sn_num)       = cx;
		        dc(sn_num)      = dcx;
		        var3(sn_num)    = var3x;
		        dvar3(sn_num)   = dvar3x;
		        cov_m_s(sn_num) = cov_m_sx;
		        cov_m_c(sn_num) = cov_m_cx;
		        cov_s_c(sn_num) = cov_s_cx;

				++sn_num;
		    }
		}

		infile.close();

		if( sn_num != size )
		    throw runtime_error("*** SNLS::ReadData() ==> Fatal error: wrong number of SNe light-curve data!");

		int I, J;
		double covx;
	//  mag_covmat
		std::ifstream file1(mag_covmat_file.c_str());

		if( !file1 )
		    throw runtime_error("*** SNLS::ReadData() ==> cannot open: " + mag_covmat_file);


#ifdef _DEBUG_SNLS_
	std::cout << "---->	SNLS: reading SNLS cov_mB_mB\n";
#endif

		I = 0;
		J = 0;
 	    std::getline(file1, line);  //  skip the comments
		while( std::getline(file1, line) ){
		    if( !Read::Is_Commented(line) && !Read::Is_Empty(line) ){
		        std::stringstream stream(line);

		        while( stream >> covx ){
		            if( J < size ){
		                cov_mB_mB(I,J) = covx;
		                ++J;
		            }
		        }

		        if( J == size ){
		            ++I;
		            J = 0;
		        }
		    }
		}

		file1.close();

	//  stretch_covmat
		std::ifstream file2(stretch_covmat_file.c_str());

		if( !file2 )
			throw runtime_error("*** SNLS::ReadData() ==> cannot open: " + stretch_covmat_file);


#ifdef _DEBUG_SNLS_
	std::cout << "---->	SNLS: reading SNLS cov_alpha_alpha\n";
#endif

		I = 0;
		J = 0;
	    std::getline(file2, line);  //  skip the comments
		while( std::getline(file2, line) ){
		    if( !Read::Is_Commented(line) && !Read::Is_Empty(line) ){
		        std::stringstream stream(line);

		        while( stream >> covx ){
		            if( J < size ){
		                cov_alpha_alpha(I,J) = covx;
		                ++J;
		            }
		        }

		        if( J == size ){
		            ++I;
		            J = 0;
		        }
		    }
		}

		file2.close();

	//  colour_covmat
		std::ifstream file3(colour_covmat_file.c_str());

		if( !file3 )
		    throw runtime_error("*** SNLS::ReadData() ==> cannot open: " + colour_covmat_file);


#ifdef _DEBUG_SNLS_
	std::cout << "---->	SNLS: reading SNLS cov_beta_beta\n";
#endif

		I = 0;
		J = 0;
	    std::getline(file3, line);  //  skip the comments
		while( std::getline(file3, line) ){
		    if( !Read::Is_Commented(line) && !Read::Is_Empty(line) ){
		        std::stringstream stream(line);

		        while( stream >> covx ){
		            if( J < size ){
		                cov_beta_beta(I,J) = covx;
		                ++J;
		            }
		        }

		        if( J == size ){
		            ++I;
		            J = 0;
		        }
		    }
		}

		file3.close();

	//  mag_stretch
		std::ifstream file4(mag_stretch_covmat_file.c_str());

		if( !file4 )
		    throw runtime_error("*** SNLS::ReadData() ==> cannot open: " + mag_stretch_covmat_file);


#ifdef _DEBUG_SNLS_
	std::cout << "---->	SNLS: reading SNLS cov_mB_alpha\n";
#endif

		I = 0;
		J = 0;
	    std::getline(file4, line);  //  skip the comments
		while( std::getline(file4, line) ){
		    if( !Read::Is_Commented(line) && !Read::Is_Empty(line) ){
		        std::stringstream stream(line);

		        while( stream >> covx ){
		            if( J < size ){
		                cov_mB_alpha(I,J) = covx;
		                ++J;
		            }
		        }

		        if( J == size ){
		            ++I;
		            J = 0;
		        }
		    }
		}

		file4.close();

	//  mag_colour
		std::ifstream file5(mag_colour_covmat_file.c_str());

		if( !file5 )
		    throw runtime_error("*** SNLS::ReadData() ==> cannot open: " + mag_colour_covmat_file);


#ifdef _DEBUG_SNLS_
	std::cout << "---->	SNLS: reading SNLS cov_mB_beta\n";
#endif

		I = 0;
		J = 0;
	    std::getline(file5, line);  //  skip the comments
		while( std::getline(file5, line) ){
		    if( !Read::Is_Commented(line) && !Read::Is_Empty(line) ){
		        std::stringstream stream(line);

		        while( stream >> covx ){
		            if( J < size ){
		                cov_mB_beta(I,J) = covx;
		                ++J;
		            }
		        }

		        if( J == size ){
		            ++I;
		            J = 0;
		        }
		    }
		}

		file5.close();

	//  stretch_colour
		std::ifstream file6(stretch_colour_covmat_file.c_str());

		if( !file6 )
		    throw runtime_error("*** SNLS::ReadData() ==> cannot open: " + stretch_colour_covmat_file);


#ifdef _DEBUG_SNLS_
	std::cout << "---->	SNLS: reading SNLS cov_alpha_beta\n";
#endif

		I = 0;
		J = 0;
	    std::getline(file6, line);  //  skip the comments
		while( std::getline(file6, line) ){
		    if( !Read::Is_Commented(line) && !Read::Is_Empty(line) ){
		        std::stringstream stream(line);

		        while( stream >> covx ){
		            if( J < size ){
		                cov_alpha_beta(I,J) = covx;
		                ++J;
		            }
		        }

		        if( J == size ){
		            ++I;
		            J = 0;
		        }
		    }
		}

		file6.close();
	}

	void SNLS::UpdateNuisance( imcmc_double par ){
		alpha   = par["alpha_snls"];
		beta    = par["beta_snls"];
		scriptMB= par["scriptMB_snls"];
	}

	void SNLS::UpdateCov(){

#ifdef _DEBUG_SNLS_
		std::cout << "---> SNLS: updating covariance matrix ..\n";
#endif

		cov_tot.fill(0);

		cov_tot += (cov_mB_mB
		        +   alpha*alpha*cov_alpha_alpha
		        +   beta*beta*cov_beta_beta
		        +   2*alpha*cov_mB_alpha
		        -   2*beta*cov_mB_beta
		        -   2*alpha*beta*cov_alpha_beta );


	//  now add the diagonal terms
		arma::mat cov_temp = arma::zeros(sne_num, sne_num);

		for(int i=0; i<sne_num; ++i){

		//  light curve fitting uncertainties
		    double cov_ii   =   dmb(i) * dmb(i)
		                    +   alpha*alpha * ds(i) * ds(i)
		                    +   beta*beta * dc(i) * dc(i)
		                    +   2*alpha * cov_m_s(i)
		                    -   2*beta * cov_m_c(i)
		                    -   2*alpha*beta * cov_s_c(i);

		//  peculiar velocity and redshift error...
		    double zfacsq   = ( 5.0/log(10.0) )*( 5.0/log(10.0) );
		    double emptyfac = ( 1.0 + zcmb(i) ) / ( zcmb(i) * (1 + 0.5*zcmb(i)));
		    double zerrsq   = dz(i)*dz(i) + pecz*pecz;
		    zerrsq *= zfacsq*emptyfac*emptyfac;
		    cov_ii += zerrsq;

		//  intrinsic dispersion
		    if( use_four_disp == true ){
		        if( set[i] == 0 )
		            cov_ii += (intrinsicdisp0 * intrinsicdisp0);
		        else if( set[i] == 1 )
		            cov_ii += (intrinsicdisp1 * intrinsicdisp1);
		        else if( set[i] == 2 )
		            cov_ii += (intrinsicdisp2 * intrinsicdisp2);
		        else if( set[i] == 3 )
		            cov_ii += (intrinsicdisp3 * intrinsicdisp3);
		    }
		    else
		        cov_ii += (intrinsicdisp * intrinsicdisp);

		//  lensing: sigma_lens = 0.055*z
		    cov_ii += pow(0.055*zcmb(i), 2);

		    cov_temp(i, i) = cov_ii;
		}

		cov_tot += cov_temp;

	//	now invert the covariance matrix
		icov_tot = cov_tot.i();

    //  ###############################################################
    //  The following is used to extract the full covariance matrix
	// 	save the diagonal terms
/*
		arma::mat diag_term = arma::zeros(sne_num,sne_num);
		double alpha_best = 1.43;
		double beta_best = 3.26;

		for(int i=0; i<sne_num; ++i){

		//  light curve fitting uncertainties
		    double cov_ii   =   dmb(i) * dmb(i)
		                    +   alpha_best*alpha_best * ds(i) * ds(i)
		                    +   beta_best*beta_best * dc(i) * dc(i)
		                    +   2*alpha_best * cov_m_s(i)
		                    -   2*beta_best * cov_m_c(i)
		                    -   2*alpha_best*beta_best * cov_s_c(i);

		//  peculiar velocity and redshift error...
		    double zfacsq   = ( 5.0/log(10.0) )*( 5.0/log(10.0) );
		    double emptyfac = ( 1.0 + zcmb(i) ) / ( zcmb(i) * (1 + 0.5*zcmb(i)));
		    double zerrsq   = dz(i)*dz(i) + pecz*pecz;
		    zerrsq *= zfacsq*emptyfac*emptyfac;
		    cov_ii += zerrsq;

		//  intrinsic dispersion
		    if( use_four_disp == true ){
		        if( set[i] == 0 )
		            cov_ii += (intrinsicdisp0 * intrinsicdisp0);
		        else if( set[i] == 1 )
		            cov_ii += (intrinsicdisp1 * intrinsicdisp1);
		        else if( set[i] == 2 )
		            cov_ii += (intrinsicdisp2 * intrinsicdisp2);
		        else if( set[i] == 3 )
		            cov_ii += (intrinsicdisp3 * intrinsicdisp3);
		    }
		    else
		        cov_ii += (intrinsicdisp * intrinsicdisp);

		//  lensing: sigma_lens = 0.055*z
		    cov_ii += pow(0.055*zcmb(i), 2);

			diag_term(i,i) = cov_ii;

			cout << "I = " << i << " Dii = " << cov_ii << endl;
		}


		diag_term += (cov_mB_mB
				+   alpha_best*alpha_best*cov_alpha_alpha
				+   beta_best*beta_best*cov_beta_beta
				+   2*alpha_best*cov_mB_alpha
				-   2*beta_best*cov_mB_beta
				-   2*alpha_best*beta_best*cov_alpha_beta );

        cout <<"==> saving SNLS3 diagonal cov terms ....\n";
        diag_term.save("snls_diag_dmu.txt",arma::raw_ascii);
        exit(0);
*/
	}

}
