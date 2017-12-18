#include <iostream>
#include <string>
#include "GSL_Math.hpp"
#include "Error.hpp"

using namespace std;

namespace GSL_Math{

    double GSL_LnDet( gsl_matrix *A, int N ){
		int signum;
		double lndetA;
		gsl_matrix *A_Copy	= gsl_matrix_alloc(N, N);
		gsl_matrix_memcpy(A_Copy, A);

		gsl_permutation *p	= gsl_permutation_alloc(N);
		gsl_linalg_LU_decomp(A_Copy, p, &signum);

	//	Not calcuate ln|A| from LU decompsition of A
		lndetA	= gsl_linalg_LU_lndet(A_Copy);

		return lndetA;
	}

    void   GSL_Inv_Matrix(gsl_matrix *A, gsl_matrix *IA, int N){
		int signum;

		gsl_matrix *A_Copy  = gsl_matrix_alloc(N, N);   //  make a copy of matrix A
		gsl_matrix_memcpy(A_Copy, A);

		gsl_permutation *p = gsl_permutation_alloc(N);  
		gsl_linalg_LU_decomp(A_Copy, p, &signum);       //  LU decomposition
		gsl_linalg_LU_invert(A_Copy, p, IA);            //  invert

		gsl_matrix_free(A_Copy);
		gsl_permutation_free(p);
	}

}
