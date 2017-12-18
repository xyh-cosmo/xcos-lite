#include <armadillo>

#include "Cosmology.hpp"

using namespace std;

namespace Cosmology{

//  functions to compute distance modulus

	double SNE_DistanceModulus( double z,
		                        CosmologyWorkspace& cw ){
		return 5.*log10( cw.background->LuminosityDistance_z(z) ) + 25.;
	}

	double SNE_DistanceModulus( double z,
		                        Background& bg ){
		return 5.*log10( bg.LuminosityDistance_z(z) ) + 25.;
	}

//	=================
//	using reference
//	=================

	void SNE_DistanceModulus(   double z[],
		                        double mu[],
		                        int sne_num,
		                        CosmologyWorkspace& cw ){
		for(int i=0; i<sne_num; ++i)
		    mu[i] = 5.*log10( cw.background->LuminosityDistance_z(z[i]) ) + 25.;
	}

	void SNE_DistanceModulus(   double z[],
		                        double mu[],
		                        int sne_num,
		                        Background& bg ){
		for(int i=0; i<sne_num; ++i)
		    mu[i] = 5.*log10( bg.LuminosityDistance_z(z[i]) ) + 25.;
	}

//	==================
//  using armadillo
//	==================
	void SNE_DistanceModulus(   arma::vec& z,
		                        arma::vec& mu,
		                        int sne_num,
		                        CosmologyWorkspace& cw ){

		for(int i=0; i<sne_num; ++i)
		    mu(i) = 5.*log10( cw.background->LuminosityDistance_z(z(i)) ) + 25.;
	}

	void SNE_DistanceModulus(   arma::vec& z,
		                        arma::vec& mu,
		                        int sne_num,
		                        Background& bg ){
		for(int i=0; i<sne_num; ++i)
		    mu(i) = 5.*log10( bg.LuminosityDistance_z(z(i)) ) + 25.;
	}


}
