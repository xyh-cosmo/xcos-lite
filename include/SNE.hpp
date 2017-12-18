#include "Cosmology.hpp"
#include <armadillo>
#include "jla.h"
#include <ctime>

using namespace std;
using namespace arma;
using namespace Cosmology;

#ifndef __SNE__
#define __SNE__

//  ============================================================================
//  Interface to cosmology
//  ============================================================================

namespace Cosmology{

//  distance modulus & apparent magnitudes:
//  NOTE: these two DO NOT include the peak-light corrections.

    double SNE_DistanceModulus( double z,
                                CosmologyWorkspace& cw  );

    double SNE_DistanceModulus( double z,
                                Background& bg  );

    void SNE_DistanceModulus(   double z[],
                                double mu[],
                                int sne_num,
                                CosmologyWorkspace& cw  );

    void SNE_DistanceModulus(   double z[],
                                double mu[],
                                int sne_num,
                                Background& bg  );

//  for using armadillo
    void SNE_DistanceModulus(   arma::vec& z,
                                arma::vec& mu,
                                int sne_num,
                                CosmologyWorkspace& cw  );

    void SNE_DistanceModulus(   arma::vec& z,
                                arma::vec& mu,
                                int sne_num,
                                Background& bg  );
}

//  ============================================================================
//  Simultion Part
//  ============================================================================
namespace Simulation{

//  generating mock SNe samples
//  use WFIRST_AFTA assumed redshift number density & error model.

	struct SNe_Mock_Nz{
	//	this structure is designed to store user provided redshifts of SNeIa

		std::string	sn_z_fname;			//	filename of the user provided redshifts
	};

	struct SNe_Mock_Nz_MagErr{
	//	this structure is designed to store user provided redshifts and errors on the magnitudes

		std::string sn_z_magerr_fname;	//	filename of the user provided redhsifts & mu-error
	};

    struct SNe_Mock_Generator{

    //  use_flux_as_obs: is true, then assume Gauusian errors on flux instead of distance modulus
        bool SN_flux_as_obs;

    //  SN factory can be combined with wfirst OR deep SN to form a new mock sample, so that
    //  there will be plenty of SNe at both low and high redshifts.
        bool    gen_snfactory;
        bool    gen_wfirst;
        bool    gen_deepsn;

        bool    fix_sne_num_in_each_bin_wfirst;

        int     sne_num_wanted_snfactory;   //  ~300 nearby low-z SNe
        int     sne_num_wanted_wfirst;      //  ~2750
        int     sne_num_wanted_deepsn;      //  ~13000

    //  these z_zmax set the maximum redshifts of the SNe in each mock samples.
        double  z_max_wanted_snfactory;     //  ~0.08
        double  z_max_wanted_wfirst;        //  ~1.7
        double  z_max_wanted_deepsn;        //  ~1.5

    //  Since the first step of generating mocks is to randomly sample the redshifts
    //  from given distribution, we use the following vector<double> to hold all
    //  the redshifts drawn from the given z-distribution(s)
        std::vector<double> sn_z;

    //  ========================================================================
    //  these are reserved for generating union2.1/snls3/jla -like mock samples
		bool	use_union_z;
		bool	use_snls_z;
		bool	use_jla_z;

		std::vector<double>	union_z;
		std::vector<double> snls_z;
		std::vector<double> jla_z;

    //  ========================================================================
        gsl_rng  *seed;
        unsigned int seed_of_seed;

    //  ========================================================================
    //  here we use WFIRST error model only, should be sufficient for most purpose
        double  err_int_fac0, err_int_fac1;     //  default values: 0.1, 0.033
        double  err_meas_fac;                   //  default value:  0.08
        double  err_lens_fac;                   //  default value:  0.07
        double  err_sys_fac;                    //  default value:  0.02

        string  outfile;
        bool    use_default_fname;
        int     mock_num;

        double  sigma_mB_int(double z);
        double  sigma_mB_meas(double z);
        double  sigma_mB_lens(double z);
        double  sigma_mB_sys(double z);
        double  sigma_mB_Tot(double z);

        SNe_Mock_Generator();
        ~SNe_Mock_Generator();

        bool gen_mock;   //  might be false if not configured properly.
        void Config( std::string mock_sne_conf );

        void Init_SNFactoy( std::string& paramfile );
        void Update_sn_z_SNFactory();

        bool include_wide_survey;
        void Init_DEEPSN( std::string& paramfile );
        void Update_sn_z_DEEPSN();

        void Init_WFIRST( std::string& paramfile );
        void Update_sn_z_WIFRST();

        void SampleZs();
        void GenMocks( CosmologyWorkspace& cw );
    };

}

//  ============================================================================
//  Data Structures
//  ============================================================================

namespace Data{

    enum SNE_DataName { _MOCK_=0, _UNION_, _SNLS_, _JLA_, _LSST_ };

    struct WFIRST_Mock{ //  simulated SNe data
    	int 		my_rank;

        int         sne_num;
        bool        initialized;
        bool        perturbed;  //  if false, then use theoretical distance modulus

        double      *sne_z;
        double      *sne_mu;
        double      *sne_dmu;

        WFIRST_Mock();
        ~WFIRST_Mock();

        void        ReadData( string paramfile );
    };

    struct SNLS_JLA_mock{
    //  this structure is written for SNLS3 and JLA mock SN samples
        int my_rank;
        int sne_num;

        bool initialized;
        bool use_full_covmat;   // if false, then use only the diagonal terms of covmat

        arma::vec z;
        arma::vec mu;
        arma::vec mB; // use mB as observable is more appropriate
        arma::vec dmu;

		arma::mat cov;	// covariance matrix
        arma::mat icov;  // inverse of the covariance matrix

        SNLS_JLA_mock();
        ~SNLS_JLA_mock();

        void ReadData( string paramfile );
        // double Chisq( CosmologyWorkspace& cw );
    };

    struct UNION{
    	int 		my_rank;
        int         sne_num;    //  580
        arma::vec   z;          //  redshifts of SNe
        arma::vec   mu, dmu;    //  distance modulus and uncertainties
        arma::vec   P;          //  probability of low mass
        arma::mat   icov;       //  covariance matrix, two are provided, nosys / sys, has aready been inverted
        bool        systematic; //  true = sys, false = nosys

        UNION(){
        	my_rank	= MPI::COMM_WORLD.Get_rank();
        }

        double      max_z();
        void        ReadData( string paramfile );
    };

    struct SNLS{
    	int 		my_rank;
        int         sne_num;            //  472

        double 	    pecz;               // = 0.001 --> 300 km/sec

        bool        use_four_disp;      // if true, then use four dispersion for each subset
        double      intrinsicdisp;      //  0.13
        double      intrinsicdisp0;     //  0.0675
        double      intrinsicdisp1;     //  0.1133
        double      intrinsicdisp2;     //  0.0815
        double      intrinsicdisp3;     //  0.0989

        bool        twoscriptmfit;      //  default is set to false
        double      scriptmcut;         //  10

        SNLS(){
        	my_rank	= MPI::COMM_WORLD.Get_rank();
        }

    //  lightcurvae data
        vector<string>  sne_name;
        vector<int>     set;

        arma::vec   zcmb, zhel, dz;
        arma::vec   mb, dmb;
        arma::vec   s, ds;
        arma::vec   c, dc;
        arma::vec   var3, dvar3;
        arma::vec   cov_m_s;
        arma::vec   cov_m_c;
        arma::vec   cov_s_c;

    //  covariance matrix,
        arma::mat   cov_mB_mB;
        arma::mat   cov_mB_alpha;
        arma::mat   cov_mB_beta;
        arma::mat   cov_alpha_alpha;
        arma::mat   cov_alpha_beta;
        arma::mat   cov_beta_beta;

    //  every time alpha, beta, scriptMB been changed, the cov_tot will be updated
        double      alpha, beta, scriptMB;
        arma::mat   cov_tot, icov_tot;
        void        UpdateCov();
        void        UpdateNuisance( imcmc::imcmc_double par );

        arma::vec   get_zcmb();
        arma::vec   get_zhel();
        double      max_z();

        void        ReadData( string paramfile );
    };

    //  ============================================================================
    //  LSST SNeIa's redshifts are measured photometrically, there will be extrac
    //  constributions to the total errors in the distance modulus.  Besides, theoretically
    //  one should take into account the uncertainties in the redshifts.
    struct SNeIa_LSST{

        int sn_num;
        double *sn_z;
        double *sn_mu;
        double *sn_err;
        double *sn_zerr;

        SNeIa_LSST();
        ~SNeIa_LSST();

        bool ReadData( std::string& paramfile );

        gsl_interp_accel *acc_ns;
        gsl_spline *spline_ns;

        arma::mat A;   // 100 x 200
        int A_nrow, A_ncol;
        double zs_min, zs_max;
        double zp_min, zp_max;
        double *zs,*zp;
        void init_matrix_A( std::string& paramfile );
        double zp_current;

        double muzperr_zmin, muzperr_zmax;
        double *muzperr_z, *muzperr;
        gsl_interp_accel *acc_muzperr;
        gsl_spline *spline_muzperr;
    };

    double photoz_mu_integrand( double z, void *param );

//  JLA is defined in another independent source file

    struct SNeIa{

    	int 			my_rank;

        ParamVector     Params;
        SNE_DataName    sne_data_name;

        UNION           *sne_union;
        SNLS            *sne_snls;
        JLA             *sne_jla;
        SNeIa_LSST      *sne_lsst;
        WFIRST_Mock     *sne_wfirst_mock;
        SNLS_JLA_mock   *sne_snls_jla_mock;

        void            Init( std::string paramfile );
        bool            initialized;

		double			MB_std;
        bool            mock_marg_MB;   // if true, then use likelihood with MB marginalized
        bool            mock_MB_use_Gaussian_Prior; // if true, use gaussian prior on MB

        SNeIa();
        ~SNeIa();
    };

}


namespace Likelihoods{

//  ========================
//  Individual likelihoods
//  ========================

    double Likelihood_SNe_Mock( imcmc::imcmc_double&    param,
                                double&                 lndet,
                                double&                 chisq,
                                void*                   model,
                                void*                   data,
                                imcmc::istate&          state );

    double Likelihood_MOCK_SNLS_JLA( imcmc::imcmc_double&    param,
                                double&                 lndet,
                                double&                 chisq,
                                void*                   model,
                                void*                   data,
                                imcmc::istate&          state );

    double Likelihood_SNE_UNION(imcmc::imcmc_double&    param,
                                double&                 lndet,
                                double&                 chisq,
                                void*                   model,
                                void*                   data,
                                imcmc::istate&          state );

    double Likelihood_SNE_SNLS( imcmc::imcmc_double&    param,
                                double&                 lndet,
                                double&                 chisq,
                                void*                   model,
                                void*                   data,
                                imcmc::istate&          state );

    double Likelihood_SNE_JLA(  imcmc::imcmc_double&    param,
                                double&                 lndet,
                                double&                 chisq,
                                void*                   model,
                                void*                   data,
                                imcmc::istate&          state );

    double likelihood_SNeIa_LSST(   imcmc::imcmc_double&    full_params,
                                    double&                 lndet,
                                    double&                 chisq,
                                    void*                   model,
                                    void*                   data,
                                    imcmc::istate&          state );
//  ========================
//  Unified likelihood
//  ========================

    double Likelihood_SNE(  imcmc::imcmc_double&    param,
                            double&                 lndet,
                            double&                 chisq,
                            void*                   model,
                            void*                   data,
                            imcmc::istate&          state );

}


#endif  //  __SNE__
