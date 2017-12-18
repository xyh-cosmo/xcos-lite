#include <iostream>
#include <fstream>
#include <stdexcept>
#include <iomanip>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_heapsort.h>
#include <gsl/gsl_statistics_double.h>

/********************************************************
    Some GSL wrappers
********************************************************/

#ifndef _MATH_HPP_
#define _MATH_HPP_

//  ================================================================================================
//  some GSL settings, for examples maximum iteration number or precisions
#define     GSL_QAGS_ITER       1000
#define     GSL_QAGS_ABSTOT     1.E-10
#define     GSL_QAGS_RELTOT     1.E-10

#define		_ROMBERG_ITER_MAX_	100
//  ================================================================================================

namespace GSL_Math{

//	================================================================================================
//	Multi-Normal random number generater
//	use gsl_linalg_cholesky_decomp( gsl_matrix *A )
//	if use armadillo, this is to be done by calling chol( X, layout ), more details see armadillo
//	documentation

	void GSL_Cholesky_Decompose( gsl_matrix *A, gsl_matrix *B );

	void GSL_Gaussian_ND(	gsl_rng		*r,
							gsl_vector	*mu,
							gsl_matrix	*sigma_chol,	//	cholesky-decomposed already, and the upper triangle is removed(set to 0)
							double		*gaussian_nd,
							bool		check=false );

//  ================================================================================================
//  Integration
    double Simple_Integrator( 	double 	(*f)(double, void *),
								double 	x1,
								double 	x2,
								void 	*param	);

    double GSL_Integrator(	double 	(*f)(double, void *),
							double 	x1,
							double 	x2,
							void 	*param );

    double GSL_Integrator( 	double 	(*f)(double, void *),
							double 	x1,
							double 	x2,
							double 	abs_eps,
							double 	rel_eps,
							void 	*param	);

//	added @ 4-16-2015
	double Romberg_Integrator( 	double 	(*f)(double, void *),
								double 	a,
								double 	b,
								void 	*param,
								double 	eps=1E-5 );

	double Romberg_Integrator( 	double 	xi[],
								double 	yi[],
								int 	size,
								double 	a,
								double 	b,
								double 	eps=1E-5 );

//  ================================================================================================
//  Linear algebra
    double GSL_LnDet( gsl_matrix *A, int N );
    void   GSL_Inv_Matrix(gsl_matrix *A, gsl_matrix *IA, int N);

//  ================================================================================================
//  Interpolation
    // class Interpolator{
    //     private:
	// 		bool				initialized;
    //         std::string         Method;     //  linear or cubic?
    //         gsl_interp_accel    *acc;
    //         gsl_spline          *spline;
    //     public:
    //         Interpolator();
    //         Interpolator( double x[], double y[], int size, std::string method="cubic" );
    //         bool Init( double x[], double y[], int size, std::string method="cubic" );
    //         double Eval(double xi);
    //         ~Interpolator();
    // };


//	new designed interpolation
	struct Array_for_Interpolation{

		gsl_interp_accel	*acc;
		gsl_spline			*spline;
		std::string			method;
		std::string			array_name;
		int					array_size;
		bool				initialized;

		Array_for_Interpolation();
		~Array_for_Interpolation();

		int  Get_Size();
		bool Check_x(double x[], int size, int& loc);	//	check whether x is increasing order! enabled
														//	when -D_DEBUG_ARRAY_INTERP_ is added to compiler-flag

		bool Init( double x[], double f[], int size, std::string mtd = "cspline" );
		void Set_Name( std::string name );
		void Get_Value( double x, double& value );
		void Write_Array_to( std::string ofile );
		double Get_Value( double x );
	};

	typedef Array_for_Interpolation  array_interp;
	typedef Array_for_Interpolation  array_1d_interp;

	// struct GSL_Sort_Pairs{
	// //	this class can be used to sort pairs of vectors into ascending or descending order,
	// //	default order: ascending
	// 	gsl_vector	*vec;
	//
	// 	GSL_Sort_Pairs();
	// 	~GSL_Sort_Pairs();
	// 	void sort_to_ascend();
	// 	void sort_to_descend();
	// }

}

#endif
