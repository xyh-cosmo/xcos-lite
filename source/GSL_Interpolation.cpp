#include <iostream>
#include <imcmc/parser++.hpp>
#include "GSL_Math.hpp"
#include "Error.hpp"

using namespace std;
using namespace imcmc::parser;

namespace GSL_Math{


	Array_for_Interpolation::Array_for_Interpolation(){
		array_size	= -1;
		initialized	= false;
	}

	Array_for_Interpolation::~Array_for_Interpolation(){
		if( initialized ){
			gsl_spline_free(spline);
			gsl_interp_accel_free(acc);
		}
	}

	int Array_for_Interpolation::Get_Size(){
		return array_size;
	}

	bool Array_for_Interpolation::Init( double x[], double f[], int size, string mtd ){

		if( array_size == -1 ){	//	not allocated
			array_size	= size;
			acc			= gsl_interp_accel_alloc();

			if( Read::SameStrings(mtd, "cspline") ){
				spline	= gsl_spline_alloc(gsl_interp_cspline, array_size);
			}
			else if( Read::SameStrings(mtd, "linear") ){
				spline	= gsl_spline_alloc(gsl_interp_linear, array_size);
			}
			else{
				print_icosmo_error;
				string err = "\n*** Array_for_Interpolation::Init( double x[], double f[], int size, string name, string method )\n";
				err += "*** currently only scpline and linear are supported!!";
				throw runtime_error(err);
			}
		}
		else if( size != array_size ){	//	allocated, by check size
			string err = "\n*** Array_for_Interpolation::Init( double x[], double f[], int size, string name, string method )\n";
			err += "*** the array to be interpolated should has the same size as array_size";
			throw std::runtime_error(err);
		}

		int status = gsl_spline_init( spline, x, f, size );

		if( status != GSL_SUCCESS ){
			cout << "failed to initialize GSL interpolatiion.\n";
			initialized = false;
			return false;
		}

		return true;
	}

	void Array_for_Interpolation::Set_Name( std::string name ){
		array_name = name;
	}

	void Array_for_Interpolation::Get_Value( double x, double& value ){
		value = gsl_spline_eval( spline, x, acc );
	}

	double Array_for_Interpolation::Get_Value( double x ){
		double value = gsl_spline_eval( spline, x, acc );
		return value;
	}


	void Array_for_Interpolation::Write_Array_to( std::string ofile ){

		ofstream of(ofile.c_str());

		if( !of.good() ){
			std::string err = "\n*** void Array_for_Interpolation::Write_Array_to( std::string ofile )\n";
			err += "*** failed to open output file: " + ofile;
			throw std::runtime_error(err);
		}

		for( int i=0; i<int(spline->size); ++i ){
			of 	<< std::setw(15) << std::setprecision(12) << spline->x[i] << " "
				<< std::setw(15) << std::setprecision(12) << spline->y[i] << endl;
		}

		of.close();
	}

}
