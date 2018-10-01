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
	ScalarType                    *V;                           // matrix with vectors {A^0v, A^1v, ..., A^kv}
	ScalarType                    *r;                           // residual vector b - Ax
	ScalarType                    *tx;                          // true solution
	ScalarType                    lambda;                       // Ritz value for Newton Basis
	ScalarType                    lambda_old;                   // old Ritz value for Newton Basis
	ScalarType                    lambda_imag;                  // imaginary part of Ritz value for Newton Basis
	ScalarType                    *H;                           //
	// ScalarType                    rTol;                         //
	// size_t                        itTol;                        //

	std::vector<pair_t, mkl_allocator<pair_t>>  theta_vals;     // ritz values in modified leja ordering	

	sparse_status_t               stat;
	sparse_matrix_t               A;                            // n x n matrix A
	MatrixInfo<ScalarType>        minfo;                        //
	size_t                        n;                            // dim(A)
	const size_t                  t = 3;                        // number of 'outer iterations' before restart
	const size_t                  m = s*t;                      // restart length
	size_t                        i, j, k;                      // indices used in for-loops

	if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
		exit(1);

	// stat = matrix_reader<ScalarType>::read_matrix_from_file("../matrix_market/matlab_example.mtx", &A, &minfo);
	// stat = matrix_reader<ScalarType>::read_matrix_from_file("../matrix_market/mini_test.mtx", &A, &minfo);
	stat = matrix_reader<ScalarType>::read_matrix_from_file("../matrix_market/rand9x9.mtx", &A, &minfo);
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
	V = (ScalarType *)mkl_calloc(n * (m+1), sizeof(ScalarType), 64);if(V == NULL){return 1;}
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

	stat = gmres<ScalarType>::init_gmres(n, A, r, H, V, theta_vals, 2*s, m);
	
	// printf("\n============= H:\n");
	// for(size_t i = 0; i < m+1; ++i) {
		// for(size_t j = 0; j < m; ++j) {
			// printf("%2.2f ", H[i*m + j]);
		// }
		// printf("\n");
	// }
	
	std::cout << std:: endl << "theta_vals (FINAL):" << std:: endl;
	for (auto p: theta_vals) {
		std::cout << p.second << "   \toutlist: " << p.first << std:: endl;
	}
	
	ScalarType **diag_ptr = minfo.diag_ptr;
	lambda_old = 0;
	
	for (i = 2*s-1; i < 3*s-1 && i < m; ++i) {

		
		lambda = theta_vals.at((i - 2*s + 1)%s).second.real();
		lambda_imag = theta_vals.at((i - 2*s + 1)%s).second.imag();
		
		if (n < 15) {
			printf("\n============= V col major:\n");
			for(k = 0; k < n; ++k) {
				for(j = 0; j < m+1; ++j) {
					std::cout << V[j*n + k] << " ";
				}
				std::cout << std::endl;
			}
		}

		#pragma omp parallel for if(n > 100000000000000) default(none) \
		shared(diag_ptr, n, lambda, lambda_old) schedule(auto)
		for (j = 0; j < n; ++j) {
			// *diag_ptr[j] = *diag_ptr[j] - lambda + lambda_old;
		}
		
		for(k = 0; k < n; ++k) {
			// *minfo.values[k] = 1000000;
		}
		
		for(k = 0; k < n; ++k) {
			std::cout << V[i*n + k] << " VECCC!!!\n";
		}
		
		// std::cout << "LAMBDAS!!!:   " << minfo.values[8*n+8] << "  ****\n";
		
		gmres<ScalarType>::mv(A, &V[i*n], &V[(i+1)*n], s);
		
		if (lambda_imag < 0) {
			// cblas_daxpy(n, lambda_imag*lambda_imag, &V[(i-1)*n], 1, &V[(i+1)*n], 1);	
		}
		
		lambda_old = lambda;
	}
	
	for (j = 0; j < n; ++j) {
		// *diag_ptr[j] = *diag_ptr[j] + lambda_old;
	}
	

	if (n < 15) {
		printf("\n============= V col major after:\n");
		for(i = 0; i < n; ++i) {
			for(j = 0; j < m+1; ++j) {
				std::cout << V[j*n + i] << " ";
			}
			std::cout << std::endl;
		}
	}

	stat = mkl_sparse_destroy(A);
	mkl_free(V);
	mkl_free(r);
	mkl_free(H);
	mkl_free(tx);
	
	mkl_free_buffers();

	return 0;
}

#endif


/*
	transpose matrix in place
*/
// void mkl_dimatcopy (const char ordering, const char trans, size_t rows, size_t cols, const double alpha, double * AB, size_t lda, size_t ldb);
// mkl_dimatcopy('R', 'T', s+1, n, 1, V, n, s+1);

// printf("\n============= V row major:\n");
// for(size_t j = 0; j < n*(s+1); ++j) {
	// std::cout << V[j] << ", ";
// }
// std::cout << std::endl;