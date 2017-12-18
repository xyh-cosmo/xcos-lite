#include "Cosmology.hpp"
#include "imcmc/parser++.hpp"

using namespace std;
using namespace Constants;
using namespace imcmc::parser;

namespace Cosmology{

//	Hubble parameter
	double Background::Hubble_z( double z ){
		double lnE	 = bg_lnE_at_z.Get_Value(z);
		return exp(lnE)*H0;
	}

//	Comoving distance
	double Background::ComovingDistance_z( double z ){
		double lnz = log(1.+z);
		return bg_dc_at_lnz.Get_Value(lnz);
	}

//	Angular distance
	double Background::AngularDiameterDistance_z( double z ){
		double lnz = log(1.+z);
		return bg_dc_at_lnz.Get_Value(lnz)/(1.0+z);
	}

//	Luminosity distance	
	double Background::LuminosityDistance_z( double z ){
		double lnz = log(1.+z);
		return bg_dc_at_lnz.Get_Value(lnz)*(1.0+z);
	}

//	Distance modulus
	double Background::DistanceModulus_z( double z ){
		double dl = (1.+z)*bg_dc_at_lnz.Get_Value( log(1.+z) );
		return 5.0*log10(dl) + 25.0;
	}
}
