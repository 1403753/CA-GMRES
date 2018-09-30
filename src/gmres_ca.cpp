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
	ScalarType                    *V;           // matrix with vectors {A^0v, A^1v, ..., A^kv}
	ScalarType                    **V_col_ptr;  // necessary for MKL spmv!
	ScalarType                    *r;           // residual vector b - Ax
	ScalarType                    *tx;          // true solution
	ScalarType                    beta;         // L2-norm of r_0
	ScalarType                    *H;           //
	ScalarType                    *Q;           //
	// ScalarType                    rTol;         //
	// size_t                        itTol;        //
	ScalarType 										**diag_ptr;

	std::vector<pair_t, mkl_allocator<pair_t>>  theta_vals;     // ritz values in modified leja ordering	

	sparse_status_t               stat;
	sparse_matrix_t               A;            // n x n matrix A
	MatrixInfo<ScalarType>        minfo;        //
	size_t                        n;            // dim(A)
	const size_t                  t = 2;        // number of 'outer iterations' before restart
	const size_t                  m = s*t;      // restart length
	size_t                        indx;         // necessary for MKL spmv! Connected to V_col_ptr.
	size_t                        i;            // index used in for-loops

	
	if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
		exit(1);
	
	stat = matrix_reader<ScalarType>::read_matrix_from_file("../matrix_market/matlab_example.mtx", &A, &minfo);
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
	Q = (ScalarType *)mkl_calloc(n * (m+1), sizeof(ScalarType), 64);if(Q == NULL){return 1;}
	H = (ScalarType *)mkl_calloc((m+1) * m, sizeof(ScalarType), 64);if(H == NULL){return 1;}	
	
	V_col_ptr = (ScalarType **)mkl_malloc(n * sizeof(ScalarType *), 64);if(V_col_ptr == NULL){return 1;}
	V = (ScalarType *)mkl_calloc(n * (s+1), sizeof(ScalarType), 64);if(V == NULL){return 1;}
	
	
	/* set true solution */
	for(i = 0; i < n; ++i)
		tx[i] = 1;
	
	/********************************************/
	/*  compute RHS b                           */
	/*  compute residual vector (r = b - Ax_0)  */
	/*  therefore, b == r                       */
	/********************************************/
	gmres<ScalarType>::mv(A, tx, r, 1);

	stat = gmres<ScalarType>::init_gmres(n, A, r, H, Q, theta_vals, 2*s, m);
	
	
	printf("\n============= H:\n");
	for(size_t i = 0; i < m+1; ++i) {
		for(size_t j = 0; j < m; ++j) {
			printf("%2.2f ", H[i*m + j]);
		}
		printf("\n");
	}
	
	std::cout << std:: endl << "theta_vals (FINAL):" << std:: endl;
	for (auto p: theta_vals) {
		std::cout << p.second << "   \toutlist: " << p.first << std:: endl;
	}	
	
	
	indx = 0;
	for(i = 0; i < n * (s+1); i += n)
		V_col_ptr[indx++] = &V[i];

	diag_ptr = minfo.diag_ptr;

	for (i = 0; i < n; ++i)
		*diag_ptr[i] -= theta_vals.at(0).second.real();


/*
	compute L2-norm beta
*/
// double cblas_dnrm2(const MKL_INT n, const double *x, const MKL_INT incx);
	beta = cblas_dnrm2(n, r, 1);
/*
	copy and scale residual vector 'r' into V matrix
*/
	cblas_dcopy(n, r, 1, *V_col_ptr, 1);
// void	cblas_dscal(const int N, const double alpha, double *X, const int incX)
	cblas_dscal(n, 1 / beta, *V_col_ptr, 1);
	
	if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
			exit(1);
	
	for(i = 0; i < s; ++i) {

		gmres<ScalarType>::mv(A, V_col_ptr[i], V_col_ptr[i+1], s);
		
		// beta = cblas_dnrm2 (n, V_col_ptr[i + 1], 1);
		// cblas_dscal (n, 1 / beta, V_col_ptr[i + 1], 1);
		
		// printf("\n============= V:\n");
		// for(size_t j = 0; j < n; ++j) {
			// for(size_t k = 0; k < s+1; ++k) {
				// std::cout << V[k*n + j] << '\t';
			// }
			// std::cout << std::endl;
		// }
	}
	
	if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
			exit(1);
	
	PAPI_shutdown();	

	// printf("\n============= V col major:\n");
	// for(size_t j = 0; j < n*(s+1); ++j) {
		// std::cout << V[j] << ", ";
	// }
	// std::cout << std::endl;

	
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

	// if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
			// exit(1);

	// tsqr::qr(&V, n, s+1);
	
	// if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
		// exit(1);
	
	// PAPI_shutdown();
	
	// arnoldi_ca::givens_rotations();
	// arnoldi_ca::modified_leja_ordering();

	stat = mkl_sparse_destroy(A);
	mkl_free(V);
	mkl_free(r);
	mkl_free(H);
	mkl_free(Q);
	mkl_free(tx);
	mkl_free(V_col_ptr);
	printf("runtime SpMV: %f\n", rtime);
	mkl_free_buffers();

	return 0;
}

#endif