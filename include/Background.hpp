#include "Verbose.hpp"
#include "GSL_Math.hpp"
#include "Constants.hpp"
#include "CosmoParams.hpp"
#include "Error.hpp"
#include "Precision.hpp"

#ifndef __BACKGROUND__
#define __BACKGROUND__

using namespace std;
using namespace GSL_Math;

namespace Cosmology{

//	default size of array for spline, in the future different background quantities might have
//	different defualt size
#ifndef _BACKGROUND_SPLINE_SIZE_
	#define _BACKGROUND_SPLINE_SIZE_ 100
#endif


#ifndef _K_ABS_MIN_
    #define _K_ABS_MIN_  1E-6
#endif

class Background{
	public:
        int     my_rank;                //  openmpi rank id

		double	H0, h;              	//  Hubble constant, H0=100*h
		double	Omega_k;            	//  curvature fraction: Omega_k = -K/H0^2
		double	K;                  	//  curvature

		int		E_table_size;
		double*	E_table;				//	used when use_tabulated_E is true
		double*	E_table_z;
		double	E_table_zmax;			//	max redshift , min value is always 0
		
		string 	E_table_interp_method;	//	"linear"; "cspline"

    //  This object holds all cosmological parameters, to communicate with "imcmc"
	    ParamVector Params;

		bool Init( std::string paramfile );		//	update init using param values in paramfile
		bool Init_to_default();					//	init to default values
		bool Init_background_array();			//	init the array_interps

		// void UpdateParam( std::string pname, double pvalue );	// NOT implemented
		bool Update();

        void PrintParams();
        void PrintParams( std::ofstream& os );

		Background();
		~Background();

		int array_size;
        array_interp    bg_lnE_at_z;    // log(E(z))
        array_interp    bg_dc_at_lnz;   // comoving distance at log(1+z)

		double Hubble_z( double z );
		double ComovingDistance_z( double z );
		double LuminosityDistance_z( double z );
		double AngularDiameterDistance_z( double z );
		double DistanceModulus_z( double z );
};

	double background_invHz( double z, void *param );
}

#endif  //  __BACKGROUND__
