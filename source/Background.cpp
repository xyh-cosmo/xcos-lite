#include "Cosmology.hpp"
#include "imcmc/parser++.hpp"
#include <iomanip>

using namespace std;
using namespace imcmc::parser;
using namespace GSL_Math;

namespace Cosmology{

	Background::Background(){
        my_rank = MPI::COMM_WORLD.Get_rank();
		Init_to_default();
	}

	Background::~Background(){
		MPI_cout("*** Background::~Background() ==> cleaning background (tabulated_E)...",my_rank);
		delete[] E_table;
		delete[] E_table_z;
	}

	bool Background::Init_to_default(){

        bg_lnE_at_z.Set_Name("bg_lnE_at_z");
        bg_dc_at_lnz.Set_Name("bg_dc_at_lnz");

        E_table = NULL;
        E_table_z = NULL;

		return true;
	}

	bool Background::Init( std::string paramfile ){

        MPI_cout("initializing Hubble constant ...", my_rank);

		if( Read::Has_Value(paramfile, "H0", "double") && Read::Has_Value(paramfile, "h", "double") ){
			std::string err = "\n*** Background::Init(std::string paramfile) ==> only one of H0, ";
			err += "h can be set in: " + paramfile;
		    throw std::runtime_error(err);
		}
		else if( (!Read::Has_Value(paramfile, "H0", "double")) && (!Read::Has_Value(paramfile, "h", "double")) ){
			std::string err = "\n*** Background::Init(std::string paramfile) ==> must have H0 or ";
			err += "h been set in: " + paramfile;
		    throw std::runtime_error(err);
		}

		if( Read::Has_Value(paramfile, "H0", "double") ){
		    H0  = Read::Read_Double_from_File(paramfile, "H0");
		    h   = H0/100.;
		    Params.AddParam_with_Value("H0", H0, __do_mcmc__);
		}
		else if( Read::Has_Value(paramfile, "h", "double") ){
		    h   = Read::Read_Double_from_File(paramfile, "h");
		    H0  = h*100.;
		    Params.AddParam_with_Value("h", h, __do_mcmc__);
		}

		E_table_size	= Read::Read_Int_from_File(paramfile, "E_bin_num");
		if( E_table_size < 2 ){
			print_icosmo_error;
			std::string err = "\n*** Background::Init(std::string paramfile) ==> bin number of tabulated E is two small,";
			err += "*** increase that number please!";
			throw std::runtime_error(err);
		}

		if( Read::Has_Key_in_File(paramfile,"E_table_interp_method") ){
			E_table_interp_method = Read::Read_String_from_File(paramfile,"E_table_interp_method");

			if( (Read::SameStrings(E_table_interp_method,"linear") == false) && (Read::SameStrings(E_table_interp_method,"cspline") == false) ){
				print_icosmo_error;
				std::string err = "\n*** Background::Init(std::string paramfile) ==> E_table_interp_method has to be linear or cspline!";
				throw std::runtime_error(err);
			}
		}
		else{
			print_icosmo_error;
			std::string err = "\n*** Background::Init(std::string paramfile) ==> E_table_interp_method is not found";
			throw std::runtime_error(err);
		}

		E_table_z 	= new double[E_table_size+1];
		E_table		= new double[E_table_size+1];

	//	now read in {zi, Ei}, note that there is no need to read {z0, E0}
		std::string var_E;
		double Et, zmax=0.0;

		E_table_z[0]	= 0.0;
		E_table[0]		= 1.0;

		bool E_table_z_is_ready = true;

		double *ztemp = NULL;

		if( E_table_size != Read::Num_of_Value_for_Key(paramfile, "E_table_z", "double") ){
			E_table_z_is_ready = false;
		}
		else{
			ztemp = new double[E_table_size];
			Read::Read_Array_of_Double_from_File(paramfile, "E_table_z", ztemp, E_table_size);

			//	check the ordering
			for( int i=0; i<E_table_size-1; ++i ){
				if( ztemp[i] >= ztemp[i+1] ){
					std::cout << "\n*** Background::Init(std::string paramfile) --> wrong ordering of input redshifts\n";
					E_table_z_is_ready = false;
				}
			}
		}

		if( E_table_z_is_ready ){
			for( int i=1; i<=E_table_size; ++i ){

				E_table_z[i] = ztemp[i-1];

				if( ztemp[i-1] >= zmax )
					zmax = ztemp[i-1];
			}
			E_table_zmax = zmax;
		}
		else{
			throw std::runtime_error("*** error happens in setting E_table_z, check your *.ini file\n");
		}

		if( ztemp != NULL )
			delete[] ztemp;

    //  add model parameters
		for( int i=1; i<=E_table_size; ++i ){
			var_E 			= "E" + Read::IntToString(i);
			Et 				= Read::Read_Double_from_File(paramfile, var_E);
			E_table[i]		= Et;
			Params.AddParam_with_Value(var_E, Et, __do_mcmc__);
		}

		//	Now get all parameters, time to initialize the array_interp's
		if( Read::Has_Value(paramfile, "background_array_size", "int") )
			array_size = Read::Read_Int_from_File(paramfile, "background_array_size");
		else
			array_size = _BACKGROUND_SPLINE_SIZE_;

	    MPI_cout("initializing background interpolating array ...", my_rank);
		return  Init_background_array();
	}


//	=============================================================
//	Note that here we assume all arrays share the same size !!!!
//	=============================================================
/*	orders:
	0) dark energy equation of state, both w(a) and weff(a)
	1) aH
	2) conformal time, age
	3) comoving distance ( this can be derived from conformal time )
	4) sound speed and sound horizon
 */

	bool Background::Init_background_array(){	//	update interpolating arrays
	
		double *ztemp	= new double[E_table_size+1];
		double *lnEtemp	= new double[E_table_size+1];

		ztemp[0]	= 0.0;
		lnEtemp[0]	= 0.0;	//	log(1.) = 0

		for( int i=0; i<=E_table_size; ++i ){
			ztemp[i]	= E_table_z[i];
			lnEtemp[i]	= log(E_table[i]);
		}

    //  log(E(z))
		if( bg_lnE_at_z.Init(ztemp, lnEtemp, E_table_size+1, E_table_interp_method) == false ) {
			cout << "@Error: failed to initialize bg_lnE_at_z using " << E_table_interp_method << " interpolation!\n";
			return false;
		}

	//	comoving distance
		double *z		= new double[array_size];
		double *lnz 	= new double[array_size];   // this is actually log(1+z)
		double *dc		= new double[array_size];

        double dlnz = (log(1.+E_table_zmax) - log(1.))/(array_size-1.);
		for( int i=0; i<array_size; ++i ){
			lnz[i] 	= i*dlnz;
			z[i]    = exp(lnz[i])-1.0;
		}
		
		z[array_size-1] = E_table_zmax; // to avoid numerical error
		
		for( int i=0; i<array_size; ++i ){
			dc[i] = _c_in_km_s_*Romberg_Integrator( background_invHz,
                                                    0.0,
								                    z[i],
								                    this,
								                    1e-6);
//            double mu = 5.*log10(dc[i]*(1+z[i])) + 25;
//            cout.precision(15);
//			cout << " lnz[" << i << "] = " << lnz[i]
//			     << " dc[" << i << "] = " << dc[i] 
//			     << " mu = " << mu << endl;
		}
//exit(0);

		if( bg_dc_at_lnz.Init(lnz, dc, array_size) == false ){
			cout << "@Error: failed to initialize bg_dc_at_lna !\n";
			return false;
		}

		delete[] z;
		delete[] lnz;
		delete[] dc;
		delete[] ztemp;
		delete[] lnEtemp;

		return true;
	}

	//  =======================================================================
	//  void Background::Update()
	//  update all cosmological paramters, especially those derived parameters
	//  =======================================================================

	bool Background::Update(){

		if( Params.Do_MCMC("H0") ){     //  if H0 is read from parameter file
		    H0  = Params["H0"];
		    h   = 0.01*H0;
		}
		else if( Params.Do_MCMC("h") ){ //  if h is read from parameter file
		    h   = Params["h"];
		    H0  = 100.*h;
		}

	//	tabulated Hubble parameter .or. cosmic expansion rate: {zi, Hi} or {zi, Ei}
		std::string var_E;
		for( int i=1; i<=E_table_size; ++i ){
			var_E		= "E" + Read::IntToString(i);
			E_table[i]	= Params[var_E];
		}

	//	after updating the parameters, now initialize the interpolating arrays.
		return Init_background_array();
	}

	double background_invHz( double z, void *param ){
		Background *bg = static_cast<Background*>(param);
//		cout << "1.0/bg->Hubble_z(z) ==> " << 1.0/bg->Hubble_z(z) << endl;
		return 1.0/bg->Hubble_z(z);
	}
}
