//============================================================================
// Name        : gmres_ca.cpp
// Author      : Robert Ernstbrunner
// Version     :
// Copyright   : Your copyright notice
// Description : ca-gmres in C++
//============================================================================

/**********/
/*  TODO  */
/**********/

/* check sparse_status_t in every class/function */
/* make template functions virtual, or leave 'em */
/* remove static functions and initialize objects */
/* change allocation of 'R' back to malloc */
/* modify givens rotations to H */

/**************/
/*  END TODO  */
/**************/

#ifndef GMRES_CA_HPP_

#include "gmres.hpp"
#include "matrix_reader.hpp"
#include <papi.h>
#include <omp.h>
#include "_hypre_utilities.h"
#include "HYPRE_struct_ls.h"

#define GMRES_CA_HPP_

#define s 3
#define ScalarType double

int main() {
	
	std::vector<ScalarType, mkl_allocator<ScalarType>> test;
	test.reserve(2);
	
	std::cout.setf(std::ios_base::fixed);
	std::cout.precision(3);
	
	float                         rtime, ptime, mflops;
	long long                     flpops;
	ScalarType                    *V;                           // basis for Krylov subspace K(A, v) = span{v, (A^1)v, (A^2)v,...,(A^s-1)v} (in col-major)
	ScalarType                    *Q;                           // orthonormal basis for Krylov subspace K(A,v), Arnoldi output (in col-major)
	ScalarType                    *R;                           // upper triangular matrix containing construction info for H (in col-major)
																															// (at each outer iteration only the last s columns of R are needed)
	ScalarType                    *H;                           // Hessenberg matrix, Arnoldi output (in row-major)
	ScalarType                    *r;                           // residual vector b - Ax
	ScalarType                    *x;                           // solution vector x
	ScalarType                    *tx;                          // true solution
	ScalarType 						        *zeta;                        // rotated beta*e1
	ScalarType                    lambda;                       // Ritz value for Newton Basis
	ScalarType                    lambda_imag;                  // imaginary part of Ritz value for Newton Basis
	// ScalarType                    rTol;                         //
	// size_t                        itTol;                        //
	bool                          restart = true;               // restart as long as value is true

	std::vector<pair_t, mkl_allocator<pair_t>>  theta_vals;     // ritz values in modified leja ordering	

	sparse_status_t               stat;
	sparse_matrix_t               A;                            // n x n matrix A
	MatrixInfo<ScalarType>        minfo;                        //
	size_t                        n;                            // dim(A)
	const size_t                  t = 3;                        // number of 'outer iterations' before restart
	const size_t                  m = s*t;                      // restart length
	const size_t                  ritz_num = s;                 // number of Ritz-values
	size_t                        i, j, k;                      // indices used in for-loops	
	
	if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
		exit(1);

	stat = matrix_reader<ScalarType>::read_matrix_from_file("../matrix_market/sparse9x9.mtx", &A, &minfo);
	// stat = matrix_reader<ScalarType>::read_matrix_from_file("../matrix_market/matlab_example.mtx", &A, &minfo);
	// stat = matrix_reader<ScalarType>::read_matrix_from_file("../matrix_market/mini_test.mtx", &A, &minfo);
	// stat = matrix_reader<ScalarType>::read_matrix_from_file("../matrix_market/goodwin.mtx", &A, &minfo);
	// stat = matrix_reader<ScalarType>::read_matrix_from_file("../matrix_market/nasa4704.mtx", &A, &minfo);
	// stat = matrix_reader<ScalarType>::read_matrix_from_file("../matrix_market/bmw7st_1.mtx", &A, &minfo);

	if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
			exit(1);
	
	PAPI_shutdown();

	printf("runtime readfile: %f\n", rtime);

	if(stat != SPARSE_STATUS_SUCCESS) throw std::invalid_argument("MatCreate failed");

	n = minfo.n;

	V = (ScalarType *)mkl_malloc(n*(s + 1) * sizeof(ScalarType), 64);if(V == NULL){return 1;}
	Q = (ScalarType *)mkl_calloc(n*(m + 1), sizeof(ScalarType), 64);if(Q == NULL){return 1;}
	R = (ScalarType *)mkl_calloc((m + s + 1)*s, sizeof(ScalarType), 64);if(R == NULL){return 1;}
	H = (ScalarType *)mkl_calloc((m + 1)*m, sizeof(ScalarType), 64);if(H == NULL){return 1;}
	r = (ScalarType *)mkl_malloc(n * sizeof(ScalarType), 64);if(r == NULL){return 1;}
	x = (ScalarType *)mkl_malloc(n * sizeof(ScalarType), 64);if(x == NULL){return 1;}
	tx = (ScalarType *)mkl_malloc(n * sizeof(ScalarType), 64);if(tx == NULL){return 1;}
	zeta = (ScalarType *)mkl_calloc(m + 1, sizeof(ScalarType), 64);if(zeta == NULL){return 1;}
	
	/* set true solution */
	for (i = 0; i < n; ++i)
		tx[i] = 1;
	
	/********************************************/
	/*  compute RHS b                           */
	/*  compute residual vector (r = b - Ax_0)  */
	/*  therefore, b == r                       */
	/********************************************/
	gmres<ScalarType>::mv(A, tx, r, 1);
	
	zeta[0] = cblas_dnrm2 (n, r, 1);

	// Q is in col-major!
	for(i = 0; i < n; ++i)
		Q[i] = r[i] / zeta[0];

	stat = gmres<ScalarType>::gmres_init(n, A, H, Q, theta_vals, ritz_num, m);
	
	std::cout << std:: endl << "theta_vals (FINAL):" << std:: endl;
	for (auto p: theta_vals) {
		std::cout << p.second << "   \toutlist: " << p.first << std:: endl;
	}
	
	printf("\n\n============= H:\n");
	for(i = 0; i < ritz_num + 1; ++i) {
		for(j = 0; j < ritz_num; ++j) {
			printf("%2.2f ", H[i*m + j]);
		}
		std::cout << std::endl;
	}	
	std::cout << std::endl;

	/**************/
  /*  reduce H  */
	/**************/
	
	stat = gmres<ScalarType>::reduce_H(H, s, m, 0, zeta);	

	std::cout << "\n\n============= zeta:\n";
	for(i = 0; i < m+1; ++i)
		std::cout << zeta[i] << std::endl;	

	std::cout << "\n\n============= H reduced:\n";
	for(i = 0; i < m+1; ++i) {
		for(j = 0; j < m; ++j) {
			printf("%2.2f ", H[i*m + j]);
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;	
	do{
		// outer iterations before restart :: k < {'t' resp. 't+1'} ::
		for (k = ritz_num / s; k < 2; k++) {
			
			// compute A*q_(s+1) and store result in V[:,0]
			lambda = theta_vals.at(0).second.real();
			gmres<ScalarType>::mv(A, &Q[k*s*n], &V[0], s);
			cblas_daxpy(n, -lambda, &Q[k*s*n], 1, &V[0], 1);		

			for (i = 0; i < s; ++i) {

				lambda = theta_vals.at(i).second.real();
				lambda_imag = theta_vals.at(i).second.imag();

				// leave in for debug, matrix may not be sparse enough
				// stat = mkl_sparse_d_create_csr (&A, SPARSE_INDEX_BASE_ZERO, n, n, minfo.rows_start, minfo.rows_end, minfo.col_indx, minfo.values);
				
				gmres<ScalarType>::mv(A, &V[i*n], &V[(i+1)*n], s);
				cblas_daxpy(n, -lambda, &V[i*n], 1, &V[(i+1)*n], 1);

				if (lambda_imag < 0) {
					cblas_daxpy(n, lambda_imag*lambda_imag, &V[(i-1)*n], 1, &V[(i+1)*n], 1);
				}
			}
			
// void cblas_dgemm (const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa,
									// const CBLAS_TRANSPOSE transb, const MKL_INT m, 
                  // const MKL_INT n, const MKL_INT k,
									// const double alpha, const double *a,
									// const MKL_INT lda, const double *b,
									// const MKL_INT ldb, const double beta,
									// double *c, const MKL_INT ldc);

			/**************************/
			/*  "BLOCK" GRAM SCHMIDT  */
			/**************************/

			cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, (s*k)+1, s, n, 1, &Q[0], n, V, n, 0, R, m + s + 1);
			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, s, (s*k)+1, -1, &Q[0], n, R, m + s + 1, 1, V, n);

			/**********/
			/*  TSQR  */
			/**********/

			gmres<ScalarType>::tsqr(V, &Q[n*(k*s + 1)], &R[s*k+1], n, s, m);
			
			
			/************************************************************************************************************************/
			/*************************************  inverse of upper triangular matrix  *********************************************/
			/*************************  https://software.intel.com/en-us/mkl-developer-reference-c-trtri ****************************/
			/* lapack_int LAPACKE_dtrtri (int matrix_layout , char uplo , char diag , lapack_int n , double * a , lapack_int lda ); */
			/************************************************************************************************************************/
			
			
		} // end for "outer iteration"

		restart = false;

	} while(restart);

	/**************************************/
	/*  solve the system (standard GMRES) */
	/**************************************/

// void cblas_dtrsv (const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const CBLAS_DIAG diag, 
                  // const MKL_INT n, const double *a, const MKL_INT lda, double *x, const MKL_INT incx);
	cblas_dtrsv (CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit, ritz_num, H, m, zeta, 1);

// void cblas_dgemv (const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE trans, const MKL_INT m, const MKL_INT n, const double alpha, const double *a, const MKL_INT lda, const double *x, const MKL_INT incx, const double beta, double *y, const MKL_INT incy);	
	cblas_dgemv (CblasColMajor, CblasNoTrans, n, ritz_num, 1, Q, n, zeta, 1, 0, x, 1);	
	
	if (n < 15) {
		std::cout << "\n\n============= y:\n";
		for(i = 0; i < m+1; ++i)
			std::cout << zeta[i] << std::endl;

		std::cout << std::endl;

		std::cout << "\n\n============= sol x:\n";
		for(i = 0; i < n; ++i)
			std::cout << x[i] << std::endl;

		std::cout << std::endl;
	}
	stat = mkl_sparse_destroy(A);
	mkl_free(V);
	mkl_free(Q);
	mkl_free(R);
	mkl_free(H);
	mkl_free(r);
	mkl_free(x);
	mkl_free(tx);
	mkl_free(zeta);

	mkl_free_buffers();

	return 0;
}

#endif

/*
	transpose matrix in place
*/
// void mkl_dimatcopy (const char ordering, const char trans, size_t rows, size_t cols, const double alpha, double * AB, size_t lda, size_t ldb);
// mkl_dimatcopy('R', 'T', s+1, n, 1, Q, n, s+1);

/*****************************************/
/*  cout 'V', 'Q' and result matrix 'R'  */
/*****************************************/

			// if (n < 15) {
				// printf("\n============= final Q (col major):\n");
				// for(size_t k = 0; k < n; ++k) {
					// for(j = 0; j < m+1; ++j) {
						// std::cout << Q[j*n + k] << " ";
					// }
					// std::cout << std::endl;
				// }
			// }
			
			// if (n < 15) {
				// printf("\n============= V (col major):\n");
				// for(size_t k = 0; k < n; ++k) {
					// for(j = 0; j < s; ++j) {
						// std::cout << V[j*n + k] << " ";
					// }
					// std::cout << std::endl;
				// }
			// }

			// printf("\n============= R (col major):\n");
			// for(size_t o = 0; o < m + s + 1; ++o) {
				// for(j = 0; j < s; ++j) {
					// std::cout << R[(m + s + 1)*j + o] << " ";
				// }
				// std::cout << std::endl;
			// }

			// printf("\n============= R_ (col major):\n");
			// for(size_t o = 0; o < s; ++o) {
				// for(j = 0; j < s; ++j) {
					// std::cout << R[(m + s + 1)*j + (o + (s*k + 1))] << " ";
				// }
				// std::cout << std::endl;
			// }
			
			// printf("\n============= Q row major:\n");
			// for(size_t j = 0; j < n*(s+1); ++j) {
				// std::cout << Q[j] << ", ";
			// }
			// std::cout << std::endl;