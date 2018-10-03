//============================================================================
// Name        : gmres_ca.cpp
// Author      : Robert Ernstbrunner
// Version     :
// Copyright   : Your copyright notice
// Description : ca-gmres in C++
//============================================================================

#ifndef GMRES_CA_HPP_

// #include "arnoldi_ca.hpp"
// #include "tsqr.hpp"
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
	
  // HYPRE_StructSolver   solver;
	
	std::vector<ScalarType, mkl_allocator<ScalarType>> test;
	test.reserve(2);
	
	std::cout.setf(std::ios_base::fixed);
	std::cout.precision(3);
	
	float                         rtime, ptime, mflops;
	long long                     flpops;
	ScalarType                    *V;                           // basis for Krylov subspace K(A, v) = span{v, (A^1)v, (A^2)v,...,(A^s-1)v}
	ScalarType                    *Q;                           // orthonormal basis for Krylov subspace K(A,v), Arnoldi output
	ScalarType                    *H;                           // Hessenberg matrix, Arnoldi output
	ScalarType                    *r;                           // residual vector b - Ax
	ScalarType                    *x;                           // solution vector x
	ScalarType                    *tx;                          // true solution
	ScalarType 						        *zeta;                        // rotated beta*e1
	ScalarType                    lambda;                       // Ritz value for Newton Basis
	ScalarType                    lambda_old;                   // old Ritz value for Newton Basis
	ScalarType                    lambda_imag;                  // imaginary part of Ritz value for Newton Basis
	// ScalarType                    rTol;                         //
	// size_t                        itTol;                        //

	std::vector<pair_t, mkl_allocator<pair_t>>  theta_vals;     // ritz values in modified leja ordering	

	sparse_status_t               stat;
	sparse_matrix_t               A;                            // n x n matrix A
	MatrixInfo<ScalarType>        minfo;                        //
	size_t                        n;                            // dim(A)
	const size_t                  t = 3;                        // number of 'outer iterations' before restart
	const size_t                  m = s*t;                      // restart length
	const size_t                  ritz_num = s;                 // number of Ritz-values
	size_t                        i, j;                         // indices used in for-loops

	ScalarType c_giv = 0;
	ScalarType s_giv = 0;
	ScalarType x1 = 0;
	ScalarType x2 = 0;

	
	
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

	tx = (ScalarType *)mkl_calloc(n, sizeof(ScalarType), 64);if(tx == NULL){return 1;}
	r = (ScalarType *)mkl_calloc(n, sizeof(ScalarType), 64);if(r == NULL){return 1;}
	x = (ScalarType *)mkl_calloc(n, sizeof(ScalarType), 64);if(x == NULL){return 1;}
	zeta = (ScalarType *)mkl_calloc(m+1, sizeof(ScalarType), 64);if(zeta == NULL){return 1;}
	V = (ScalarType *)mkl_calloc(n * (s), sizeof(ScalarType), 64);if(V == NULL){return 1;}
	Q = (ScalarType *)mkl_calloc(n * (m+1), sizeof(ScalarType), 64);if(Q == NULL){return 1;}
	H = (ScalarType *)mkl_calloc((m+1) * m, sizeof(ScalarType), 64);if(H == NULL){return 1;}		
	
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


	/*****************************/
  /*  classic givens_rotation  */
	/*****************************/
	
	for (i = 0; i < ritz_num; ++i) {
			x1 = H[i*m + i];
			x2 = H[(i + 1)*m + i];
			
// void cblas_drotg (double *x1, double *x2, double *c, double *s);
			cblas_drotg(&x1, &x2, &c_giv, &s_giv);

// void cblas_drot (const size_t n, double *x, const size_t incx, double *y, const size_t incy, const double c, const double s);
			cblas_drot(m - i, &H[i*m + i], 1, &H[(i + 1)*m + i], 1, c_giv, s_giv);

// void cblas_drot (const size_t n, double *x, const size_t incx, double *y, const size_t incy, const double c, const double s);
			cblas_drot(1, &zeta[i], 1, &zeta[i + 1], 1, c_giv, s_giv);
	}


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
	
	
	// void cblas_dtrsv (const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const CBLAS_DIAG diag, 
                  // const MKL_INT n, const double *a, const MKL_INT lda, double *x, const MKL_INT incx);
	cblas_dtrsv (CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit, ritz_num, H, m, zeta, 1);
	

// void cblas_dgemv (const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE trans, const MKL_INT m, const MKL_INT n, const double alpha, const double *a, const MKL_INT lda, const double *x, const MKL_INT incx, const double beta, double *y, const MKL_INT incy);	
	cblas_dgemv (CblasColMajor, CblasNoTrans, n, ritz_num, 1, Q, n, zeta, 1, 0, x, 1);	
	
	
	std::cout << "\n\n============= y:\n";
	for(i = 0; i < m+1; ++i)
		std::cout << zeta[i] << std::endl;

	std::cout << std::endl;

	std::cout << "\n\n============= sol x:\n";
	for(i = 0; i < n; ++i)
		std::cout << x[i] << std::endl;

	std::cout << std::endl;
	
	
	
	// ScalarType **diag_ptr = minfo.diag_ptr;
	lambda_old = 0;
	
	for (i = ritz_num - 1; i < ritz_num + s - 1 && i < m; ++i) {

		lambda = theta_vals.at((i - ritz_num + 1)%s).second.real();
		lambda_imag = theta_vals.at((i - ritz_num + 1)%s).second.imag();
				
		// if (n < 15) {
			// printf("\n============= Q col major:\n");
			// for(size_t k = 0; k < n; ++k) {
				// for(j = 0; j < m+1; ++j) {
					// std::cout << Q[j*n + k] << " ";
				// }
				// std::cout << std::endl;
			// }
		// }		

		// #pragma omp parallel for if(n > 100000000000000) default(none) shared(A, diag_ptr, n, lambda) schedule(auto)
		for (j = 0; j < n; ++j) {
			*minfo.diag_ptr[j] = *minfo.diag_ptr[j] - lambda + lambda_old;
		}

		// stat = mkl_sparse_d_create_csr (&A, SPARSE_INDEX_BASE_ZERO, n, n, minfo.rows_start, minfo.rows_end, minfo.col_indx, minfo.values);
		
		gmres<ScalarType>::mv(A, &Q[i*n], &Q[(i+1)*n], s);
		
		if (lambda_imag < 0) {
			
			cblas_daxpy(n, lambda_imag*lambda_imag, &Q[(i-1)*n], 1, &Q[(i+1)*n], 1);	
		}
		lambda_old = lambda;
	}
	
	for (j = 0; j < n; ++j) {
		*minfo.diag_ptr[j] = *minfo.diag_ptr[j] + lambda_old;
	}

	if (n < 15) {
		printf("\n============= Q col major after:\n");
		for(i = 0; i < n; ++i) {
			for(j = 0; j < m+1; ++j) {
				std::cout << Q[j*n + i] << " ";
			}
			std::cout << std::endl;
		}
	}

	stat = mkl_sparse_destroy(A);
	mkl_free(V);
	mkl_free(H);
	mkl_free(Q);
	mkl_free(r);
	mkl_free(tx);
	
	mkl_free_buffers();

	return 0;
}

#endif

/*
	transpose matrix in place
*/
// void mkl_dimatcopy (const char ordering, const char trans, size_t rows, size_t cols, const double alpha, double * AB, size_t lda, size_t ldb);
// mkl_dimatcopy('R', 'T', s+1, n, 1, Q, n, s+1);

// printf("\n============= Q row major:\n");
// for(size_t j = 0; j < n*(s+1); ++j) {
	// std::cout << Q[j] << ", ";
// }
// std::cout << std::endl;