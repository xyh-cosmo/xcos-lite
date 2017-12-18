#include <armadillo>
#include "Cosmology.hpp"
#include "SNE.hpp"

using namespace std;
using namespace imcmc;
using namespace Cosmology;
using namespace Data;

namespace Likelihoods{

	double Likelihood_SNE_MOCK( imcmc_double&   param,
                                double&         lndet,
                                double&         chisq,
                                void*           model,
                                void*           data,
                                istate&         state ){

        state.this_like_is_ok = true;
		lndet = chisq = 0;
		CosmologyWorkspace  *cw = static_cast<CosmologyWorkspace*>(model);

        //  ============================================================================================================
	    //  update model parameters, in this case no LC-nuisance parameters are present, and now we need only background
        //  ============================================================================================================

		SNeIa *sne = static_cast<SNeIa*>(data);

		if( sne->mock_marg_MB == false ){
			double *mu = new double[sne->sne_wfirst_mock->sne_num];
			SNE_DistanceModulus( sne->sne_wfirst_mock->sne_z, mu, sne->sne_wfirst_mock->sne_num, *cw );

			double MB = param["MB"];
			for(int i=0; i<sne->sne_wfirst_mock->sne_num; ++i){	//	mB = mu + MB
			    double dmB = (mu[i]+MB) - (sne->sne_wfirst_mock->sne_mu[i]+(-19.3));
			    chisq += pow(dmB/sne->sne_wfirst_mock->sne_dmu[i], 2);
			}

			delete[] mu;
		}
		else{	//	the following calculation has marginalized MB
			double a = 0.0;
			double b = 0.0;
			double H0 = cw->background->H0;
			double H0_c = H0/(1e-3*_c_);
			for( int i=0; i<sne->sne_wfirst_mock->sne_num; ++i ){
			//	when using only SN data, H0=70 is assumed, be careful!
				double dL = cw->background->ComovingDistance_z(sne->sne_wfirst_mock->sne_z[i])
							* (1+sne->sne_wfirst_mock->sne_z[i]) * H0_c; // remove H0,c dependence
				double mB_th  = 5.0*log10(dL);   // see appendix-B of A. Conley for the details
				double mB_obs = (sne->sne_wfirst_mock->sne_mu[i]+(-19.3));	// convert distance modulus to apparent magnitude, assumed MB=-19.3
				double dmB = mB_obs - mB_th;
				chisq += pow(dmB/sne->sne_wfirst_mock->sne_dmu[i], 2);
				a += 1.0/pow(sne->sne_wfirst_mock->sne_dmu[i], 2);
				b += dmB/pow(sne->sne_wfirst_mock->sne_dmu[i], 2);
			}

			chisq += ( log(0.5*a/_PI_) - b*b/a );
		}

		return -lndet -0.5*chisq;
	}

	double Likelihood_MOCK_SNLS_JLA( imcmc::imcmc_double&    param,
	                                double&                 lndet,
	                                double&                 chisq,
	                                void*                   model,
	                                void*                   data,
	                                imcmc::istate&          state ){
	//  @4-25-2017
	//  This likelihood function is for SNLS3 & JLA mock SN sample, but with covariance
	//	between different SNe, unlike the case for WFIRST mock SN sample.
		state.this_like_is_ok = true;
		lndet = chisq = 0;
		CosmologyWorkspace  *cw = static_cast<CosmologyWorkspace*>(model);
		SNeIa *sne = static_cast<SNeIa*>(data);

		double MB = param["MB"];
		arma::vec mu = arma::zeros(sne->sne_snls_jla_mock->sne_num);
		SNE_DistanceModulus( sne->sne_snls_jla_mock->z, mu, sne->sne_snls_jla_mock->sne_num, *cw );
		arma::vec mB = mu + MB;
		arma::vec dmB = mB - sne->sne_snls_jla_mock->mB;

		chisq = arma::as_scalar( dmB.t() * sne->sne_snls_jla_mock->icov * dmB );

		if( sne->mock_MB_use_Gaussian_Prior )
			chisq += pow((MB+19.3)/sne->MB_std,2);

		return -lndet - 0.5*chisq;
	}

	double Likelihood_SNE_UNION( imcmc_double&  param,
                                 double&        lndet,
                                 double&        chisq,
                                 void*          model,
                                 void*          data,
                                 istate&        state ){

        state.this_like_is_ok = true;
		lndet = chisq = 0;

        CosmologyWorkspace  *cw  = static_cast<CosmologyWorkspace*>(model);
		SNeIa               *sne = static_cast<SNeIa*>(data);

		arma::vec mu = arma::zeros(sne->sne_union->sne_num);
		SNE_DistanceModulus( sne->sne_union->z, mu, sne->sne_union->sne_num, *cw );
		arma::vec dmu = mu - sne->sne_union->mu;

        //  ============================================================================
        //	remember that Union2.1 data has assumed H0=70 km/s/Mpc, so the theoretical
		//	calculation should also taking into account of this by adding 5*log10(bg->H0/70)
		//	to mu Union2.1 data is taken from cosmomc, in which H0 has been marginalized,
		//	you may check that change H0 will not change chi2 .., so I comment out the
		//	following 're-scaling relation' mu = mu + 5.0*log10(cw->background->H0/70.);
        //  =============================================================================

		chisq = arma::as_scalar( dmu.t() * sne->sne_union->icov * dmu );

		return -lndet - 0.5*chisq;
	}

	double Likelihood_SNE_SNLS( imcmc_double&   param,
                                double&         lndet,
                                double&         chisq,
                                void*           model,
                                void*           data,
                                istate&         state ){

        state.this_like_is_ok = true;

		lndet = chisq = 0;

        //  ================================
	    //  update cosmological parameters
        //  ================================
		CosmologyWorkspace  *cw  = static_cast<CosmologyWorkspace*>(model);
		SNeIa               *sne = static_cast<SNeIa*>(data);

		sne->sne_snls->UpdateNuisance(param);   // update LC nuisance parameters, now we need alpha, beta and MB
		sne->sne_snls->UpdateCov();

		arma::vec mb   = arma::zeros(sne->sne_snls->sne_num);;

        //  =========================================================================
	    //  DONOT forget to add the corrections
	    //  here we do not use "SNE_DistanceModulus( z, mu, sne_snls->size(), *cw, 70. )",
		//	since we need to use zcmb and zhel, and H0 dependence has been removed.
        //  =========================================================================

		for(int i=0; i<sne->sne_snls->sne_num; ++i){
		    double dl_hel = cw->background->ComovingDistance_z(sne->sne_snls->zcmb(i)) * (1+sne->sne_snls->zhel(i));
	//        dl_hel *= (cw->background->H0/(1e-3*_c_));   //  remove H0 dependence.
		    mb(i)   = 5.*log10(dl_hel)  // normal mu(z)
		            - sne->sne_snls->alpha*(sne->sne_snls->s(i)-1)
		            + sne->sne_snls->beta*sne->sne_snls->c(i)
		            + sne->sne_snls->scriptMB
		            + 5.*log10(70./(1e-3*_c_)); //  adding this terms allows to fit H0 by fixing scriptMB.

	//  scriptMB = MB + 5log10(c/H0), the extra factor 25 cancles each other.
	 //
	//        mb(i)   = 5.*log10(dl_hel) + 25  // normal mu(z)
	//                - sne->sne_snls->alpha*(sne->sne_snls->s(i)-1)
	//                + sne->sne_snls->beta*sne->sne_snls->c(i)
	//                + sne->sne_snls->scriptMB
	//                + (5.*log10(70./(1e-3*_c_)) -25); //  adding this terms allows to fit H0 by fixing scriptMB
		}

		arma::vec dmb = mb - sne->sne_snls->mb;
		chisq = arma::as_scalar( dmb.t() * sne->sne_snls->icov_tot * dmb );

		return -lndet - 0.5*chisq;
	}

    //  ================================================================================================================
	//	JLA likelihood function wrapper (v4)
    //  NOTE: JLA handles intrinsic dispersions (as well as sigma_z and sigma_lens) differently from SNLS3 ...
    //  more details can be found in the head of supernovae_JLA.f90 of cosmomc
    //  ================================================================================================================
	double Likelihood_SNE_JLA(  imcmc_double&   param,
                                double&         lndet,
                                double&         chisq,
                                void*           model,
                                void*           data,
                                istate&         state ){

        state.this_like_is_ok = true;

		lndet = chisq = 0;

        //  ============================================================================================================
	    //  update model parameters
        //  by default, pass a CosmologyWorkspace object
        //  ============================================================================================================
		CosmologyWorkspace  *cw  = static_cast<CosmologyWorkspace*>(model);

#if defined( _UPDATE_IN_EACH_LIKELIHOOD_ )
		cw->Update(param);
#endif

		SNeIa               *sne = static_cast<SNeIa*>(data);

        //  ============================================================================================================
	    //  update LC nuisance parameters
        //  NOTE: JLA::UpdateCov() actually doing nothing, all real calculation is done inside JLA
        //        code provided by "Marc Betoule"
        //  ============================================================================================================
		sne->sne_jla->UpdateNuisance(param);
		sne->sne_jla->UpdateCov();

		double alpha    = sne->sne_jla->alpha;
		double beta     = sne->sne_jla->beta;
		double MB       = sne->sne_jla->MB;
		double DeltaMB  = sne->sne_jla->DeltaMB;

		double nuisance[4] = { alpha, beta, MB, DeltaMB };

		double *z   = sne->sne_jla->getZ();
		double *mu	= new double[sne->sne_jla->size()];

		SNE_DistanceModulus( z, mu, sne->sne_jla->size(), *cw);

		chisq = sne->sne_jla->computeLikelihood(mu, nuisance);

		delete[] mu;
		delete[] z;

#if defined(_DEBUG_JLA_)
		cout << "_DEBUG_JLA_: chisq = " << chisq << endl;
#endif

		return -lndet - 0.5*chisq;
	}

	double Likelihood_SNE_LSST(	imcmc_double&   full_params,
                                double&         lndet,
                                double&         chisq,
                                void*           model,
                                void*           data,
                                istate&         state ){
	    lndet = chisq = 0.0;

		CosmologyWorkspace  *cw  = static_cast<CosmologyWorkspace*>(model);
		SNeIa               *sne = static_cast<SNeIa*>(data);

	    // CosmoTheory*    theory = static_cast<CosmoTheory*>(model);
	    // SNeIa_LSST*     lsst = static_cast<SNeIa_LSST*>(data);

	    double dmu;
	    double *muzs = new double[sne->sne_lsst->A_ncol];
	    double *muzp = new double[sne->sne_lsst->A_nrow];

	    gsl_spline *spline_muzp = gsl_spline_alloc(gsl_interp_cspline, sne->sne_lsst->A_nrow);
	    gsl_interp_accel *acc_muzp = gsl_interp_accel_alloc();

	    // for( int i=0; i<lsst->A_ncol; ++i ){
	    //     muzs[i] = 5.0*log10( theory->engine->get_Dl(lsst->zs[i]) ) + 25.0;
	    // }

		SNE_DistanceModulus( sne->sne_lsst->zs, muzs, sne->sne_lsst->A_ncol, *cw );

	    for( int i=0; i<sne->sne_lsst->A_nrow; ++i ){
	        muzp[i] = 0.0;
	        for( int j=0; j<sne->sne_lsst->A_ncol; ++j ){
	            muzp[i] += sne->sne_lsst->A(i,j)*muzs[j];
	        }
	    }

	    gsl_spline_init(spline_muzp, sne->sne_lsst->zp, muzp, sne->sne_lsst->A_nrow);

	    double MB = full_params["MB"];
	    double MB_a = full_params["MB_a"];
	    double MB_b = full_params["MB_b"];

	    double mu, muerr;
	    double MB_z;

	    for( int i=0; i<sne->sne_lsst->sn_num; ++i ){

	        MB_z = MB*( 1.0+ sne->sne_lsst->sn_z[i]*(MB_a + MB_b*sne->sne_lsst->sn_z[i]) );
	        mu = gsl_spline_eval(spline_muzp, sne->sne_lsst->sn_z[i], acc_muzp);

	        if( sne->sne_lsst->sn_z[i] <= sne->sne_lsst->muzperr_zmin ){
				muerr = gsl_spline_eval(sne->sne_lsst->spline_muzperr,
										sne->sne_lsst->muzperr_zmin,
										sne->sne_lsst->acc_muzperr);
			}
	        else if( sne->sne_lsst->sn_z[i] >= sne->sne_lsst->muzperr_zmax ){
				muerr = gsl_spline_eval(sne->sne_lsst->spline_muzperr,
										sne->sne_lsst->muzperr_zmax,
										sne->sne_lsst->acc_muzperr);
			}
	        else{
				muerr = gsl_spline_eval(sne->sne_lsst->spline_muzperr,
										sne->sne_lsst->sn_z[i],
										sne->sne_lsst->acc_muzperr);
			}

	//      NOTE: the mock SNe data provides {zi,mu(zi)}, while the true observable is mb, the apparent magnitude.
	//      cosmological theory predicts only distance modulus, mu(z), we add to it the UNKNOWN absolute magnitude
	//      MB(z), so we can predict mb(z).

	        double mb = mu + MB_z;

	        dmu = mb - ( sne->sne_lsst->sn_mu[i] + (-19.3) );
	        chisq += pow(dmu/muerr,2);
	    }

	    delete[] muzs;
	    delete[] muzp;
	    gsl_spline_free(spline_muzp);
	    gsl_interp_accel_free(acc_muzp);

	    return -lndet -0.5*chisq;
	}

	//  ================================================================================================================
	//  a unified interface
	//  void *data is casted into SNeIa pointer sne*, then one can determine from sne_ia->sne_data_name
	//  that which one of the four likelihood functions will be used.
    //  ================================================================================================================

	double Likelihood_SNE(  imcmc_double&   param,
                            double&         lndet,
                            double&         chisq,
                            void*           model,
                            void*           data,
                            istate&         state ){

        state.this_like_is_ok = true;
		chisq = 0.0;
		lndet = 0.0;
		double lnpost=0;

		SNeIa *sne = static_cast<SNeIa*>(data);

		if( sne->sne_data_name == _MOCK_ ){

			double chisq_temp, lndet_temp;

			if( sne->sne_wfirst_mock != NULL ){
				//cout << "calling wfirst SNe likelihood ...\n";
				lnpost = Likelihood_SNE_MOCK( param,
											  lndet_temp,
											  chisq_temp,
											  model,
											  data,
											  state );
				chisq += chisq_temp;
				lndet += lndet_temp;
			}

			// combine mock JLA/SNLS3 with WFIRST mock
			if( sne->sne_snls_jla_mock != NULL ){
				//cout << "calling snls_jla mock SNe likelihood ...\n";
	  			lnpost = Likelihood_MOCK_SNLS_JLA(	param,
	  												lndet_temp,
	  												chisq_temp,
	  												model,
	  												data,
	  												state );
				chisq += chisq_temp;
				lndet += chisq_temp;
	  		}

			if( sne->mock_MB_use_Gaussian_Prior ){	// we always assume MB = -19.3 as the fiducial value.
				// cout << "adding MB GP...\n";
				double MB = param["MB"];
				chisq += pow((MB+19.3)/sne->MB_std,2);
			}

//			exit(0);
        }
		else if( sne->sne_data_name == _UNION_ ){
		    lnpost = Likelihood_SNE_UNION( param,
                                           lndet,
                                           chisq,
                                           model,
                                           data,
                                           state );
        }
		else if( sne->sne_data_name == _SNLS_ ){
		    lnpost = Likelihood_SNE_SNLS( param,
                                          lndet,
                                          chisq,
                                          model,
                                          data,
                                          state );
        }
		else if( sne->sne_data_name == _JLA_ ){
		    lnpost = Likelihood_SNE_JLA( param,
                                         lndet,
                                         chisq,
                                         model,
                                         data,
                                         state );
		}
		else if( sne->sne_data_name == _LSST_ ){
		    lnpost = Likelihood_SNE_LSST( param,
                                          lndet,
                                          chisq,
                                          model,
                                          data,
                                          state );
		}
		else{
		    string err = "\n*** double Likelihood_SNE() ==> error happened when just to call one of the\n";
		          err += "*** likelihood function:\n";
		          err += "--> 1) Likelihood_SNE_MOCK()\n";
				  err += "--> 2) Likelihood_MOCK_SNLS_JLA()\n";
		          err += "--> 3) Likelihood_SNE_UNION()\n";
		          err += "--> 4) Likelihood_SNE_SNLS()\n";
		          err += "--> 5) Likelihood_SNE_JLA()\n";
				  err += "--> 6) Likelihood_SNE_LSST()\n\n";
		    throw runtime_error(err);
		}

#if defined(_DEBUG_SNE_LIKE_)
		cout << "_DEBUG_SNE_LIKE_: chisq = " << chisq << endl;
#endif

//		return lnpost;
		return -lndet - 0.5*chisq;
	}

}
