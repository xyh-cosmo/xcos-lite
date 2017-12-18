/*	
 *	algorithm is taken from: http://en.wikipedia.org/wiki/Multivariate_normal_distribution
 *
 *	1) Find any real matrix A such that A AT = Σ. When Σ is positive-definite, the Cholesky 
 *	decomposition is typically used, and the extended form of this decomposition can always be used 
 *	(as the covariance matrix may be only positive semi-definite) in both cases a suitable matrix A 
 *	is obtained. An alternative is to use the matrix A = UΛ½ obtained from a spectral decomposition 
 *	Σ = UΛUT of Σ. The former approach is more computationally straightforward but the matrices A 
 *	change for different orderings of the elements of the random vector, while the latter approach 
 *	gives matrices that are related by simple re-orderings. In theory both approaches give equally 
 *	good ways of determining a suitable matrix A, but there are differences in computation time.
 *
 *	2) Let z = (z1, ..., zN)^T be a vector whose components are N independent standard normal variates 
 *	(which can be generated, for example, by using the Box–Muller transform).
 *
 *	3) Let x be μ + Az. This has the desired distribution due to the affine transformation property.
 */

#include <iostream>
#include <stdexcept>
#include "GSL_Math.hpp"

using namespace std;

namespace GSL_Math{

//	================================================================================================
//	Multi-Normal random number generater
//	use gsl_linalg_cholesky_decomp( gsl_matrix *A )
//	if use armadillo, this is to be done by calling chol( X, layout ), more details see armadillo
//	documentation.

	void GSL_Cholesky_Decompose( gsl_matrix *A, gsl_matrix *B ){
		if( (B->size1 != A->size1) || (B->size2 != A->size2) ){
			std::string err = "\n*** void GSL_Cholesky_Decompose( gsl_matrix *A, gsl_matrix *B ) ==>";
			err += "*** B must have the same size with A!";
			throw std::runtime_error(err);
		}
		else if( B->size1 != B->size2 ){
			std::string err = "\n*** void GSL_Cholesky_Decompose( gsl_matrix *A, gsl_matrix *B ) ==>";
			err += "*** A and B must be squre matrix!";
			throw std::runtime_error(err);
		}

		gsl_matrix_memcpy(B, A);
		gsl_linalg_cholesky_decomp(B);

		for( unsigned int i=0; i<B->size1; ++i ){
			for( unsigned int j=i+1; j<B->size2; ++j ){
				gsl_matrix_set(B, i, j, 0.0);
			}
		}
	}

	void GSL_Gaussian_ND(	gsl_rng		*r,
							gsl_vector	*mu,
							gsl_matrix	*sigma_chol,	//	cholesky-decomposed already, upper triangle is removed
							double		*gaussian_nd,
							bool		check ){
		if( check ){
			if( (mu->size != sigma_chol->size1) || (mu->size != sigma_chol->size2) ){
				std::string err = "\n*** void GSL_Gaussian_ND(...) ==> check the dimension of mu and";
				err += "*** sigma_chol!!!";
				throw std::runtime_error(err);
			}

			for( unsigned int i=0; i<sigma_chol->size1; ++i ){
				for( unsigned int j=i+1; j<sigma_chol->size2; ++j ){
					if( fabs(gsl_matrix_get(sigma_chol, i, j)) > 1E-15 ){
						std::string err = "\n*** void GSL_Gaussian_ND(...) ==> check your input ";
						err += "*** Cholesky decomposed matrix!";
						throw std::runtime_error(err);
					}
				}
			}
		}

		int size = mu->size;

		gsl_vector *z = gsl_vector_alloc(mu->size);
		gsl_vector *x = gsl_vector_alloc(mu->size);

		for( int i=0; i<size; ++i )
			gsl_vector_set( z, i, gsl_ran_gaussian(r, 1.0) );

		gsl_blas_dgemv( CblasNoTrans, 1.0, sigma_chol, z, 0, x );

		gsl_vector_add( x, mu );

		for( int i=0; i<size; ++i )
			gaussian_nd[i] = gsl_vector_get(x, i);

		gsl_vector_free(z);
		gsl_vector_free(x);
	}

}
