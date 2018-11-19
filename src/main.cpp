//============================================================================
// Name        : main.cpp
// Author      : Robert Ernstbrunner
// Version     :
// Copyright   : Your copyright notice
// Description : main in C++
//============================================================================

/**********/
/*  TODO  */
/**********/

/* check sparse_status_t in every class/function */
/* check fees for allocated arrays in blackbox functions*/
/* make template functions virtual, or leave 'em */
/* remove static functions and initialize objects */
/* change allocation of 'R' back to malloc */
/* change allocation of 'H_reduced' back to malloc */
/* implement standard GMRES */

/**************/
/*  END TODO  */
/**************/

/*********/
/* notes */
/*********/
/* givens works, but produces switched signs for different values compared to matlab. the zetas are also different so the solution y is still the same */
/* all matrices are stored in col-major */

#ifndef MAIN_HPP

#include "gmres_ca.hpp"
#include "matrix_reader.hpp"
#include <papi.h>
#include <omp.h>

#include <cstddef> /* NULL */
#include <metis.h>

#define MAIN_HPP

#define s 5
#define ScalarType double

int main() {
	
	const idx_t nV = 6;
	
	idx_t nVertices = 6;
	const idx_t nEdges    = 7;
	idx_t nWeights  = 1;
	idx_t nParts    = 2;

	idx_t objval;
	idx_t part[nV];


	// Indexes of starting points in adjacent array
	idx_t xadj[nV+1] = {0,2,5,7,9,12,14};

	// Adjacent vertices in consecutive index order
	idx_t adjncy[2 * nEdges] = {1,3,0,4,2,1,5,0,4,3,1,5,4,2};
    

	// int ret = METIS_PartGraphRecursive(&nVertices,& nWeights, xadj, adjncy,
	// 				       NULL, NULL, NULL, &nParts, NULL,
	// 				       NULL, NULL, &objval, part);

	int ret = METIS_PartGraphKway(&nVertices,& nWeights, xadj, adjncy,
				       NULL, NULL, NULL, &nParts, NULL,
				       NULL, NULL, &objval, part);

  std::cout << ret << std::endl;
    
  for(unsigned part_i = 0; part_i < nVertices; part_i++){
		std::cout << part_i << " " << part[part_i] << std::endl;
  }
	
	
	
	// std::vector<ScalarType> test;
	// test.reserve(2);
	
	std::cout.setf(std::ios_base::fixed);
	std::cout.precision(5);
	
	float                         rtime, ptime, mflops;
	long long                     flpops;
	ScalarType                    *V;                             // basis for Krylov subspace K(A, v) = span{v, (A^1)v, (A^2)v,...,(A^s-1)v}
	ScalarType                    *Q;                             // orthonormal basis for Krylov subspace K(A,v), Arnoldi output
	ScalarType                    *R;                             // upper triangular matrix containing construction info for H
                                                                   // (at each outer iteration only the last s columns of R are needed)
	ScalarType                    *H;                             // Hessenberg matrix, Arnoldi output
	ScalarType                    *H_reduced;                     // QR factorized H, upper triangular matrix
  ScalarType                    *R_k;                           // s x s submatrix needed for updating H
	ScalarType                    *r;                             // residual vector b - Ax
	ScalarType                    *x;                             // solution vector x
	ScalarType                    *tx;                            // true solution
	ScalarType 						        *zeta;                          // rotated beta*e1
	ScalarType                    imag;                           // imaginary part of Ritz value for Newton Basis	
	// ScalarType                    rTol;                           //
	// ScalarType                    aTol;                           //
	// ScalarType                    dTol;                           //
	// size_t                        maxit;                          //
	bool                          restart = true;                 // restart as long as value is true

	std::vector<ic_pair_t>  theta_vals; // ritz values in modified leja ordering	

	sparse_status_t               stat;
	sparse_matrix_t               A;                              // n x n matrix A
	sparse_matrix_t               M;                              // n x n preconditioned matrix M == ilu0(A)
	MatrixInfo<ScalarType>        minfo;                          //
	size_t                        n;                              // dim(A)
	const size_t                  t = 12;                         // number of 'outer iterations' before restart
	const size_t                  m = s*t;                        // restart length
	const size_t                  ritz_num = s;                   // number of Ritz-values
	size_t                        i, k;                           // indices used in for-loops	
	
	if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
		exit(1);

	stat = matrix_reader<ScalarType>::read_matrix_from_file("../matrix_market/sparse9x9complex.mtx", &A, &minfo);
	// stat = matrix_reader<ScalarType>::read_matrix_from_file("../matrix_market/matlab_example.mtx", &A, &minfo);
	// stat = matrix_reader<ScalarType>::read_matrix_from_file("../matrix_market/mini_test.mtx", &A, &minfo);
	// stat = matrix_reader<ScalarType>::read_matrix_from_file("../matrix_market/goodwin.mtx", &A, &minfo);
	// stat = matrix_reader<ScalarType>::read_matrix_from_file("../matrix_market/nasa4704.mtx", &A, &minfo);
	// stat = matrix_reader<ScalarType>::read_matrix_from_file("../matrix_market/dwb512.mtx", &A, &minfo);
	// stat = matrix_reader<ScalarType>::read_matrix_from_file("../matrix_market/bmw7st_1.mtx", &A, &minfo);

	if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
			exit(1);

	PAPI_shutdown();

	printf("runtime readfile: %f\n", rtime);

	if(stat != SPARSE_STATUS_SUCCESS) throw std::invalid_argument("MatCreate failed");

	n = minfo.n;

	V = (ScalarType *)mkl_malloc(n*s * sizeof(ScalarType), 64);if(V == NULL){return 1;}
	Q = (ScalarType *)mkl_calloc(n*(m + 1), sizeof(ScalarType), 64);if(Q == NULL){return 1;}
	R = (ScalarType *)mkl_calloc((m + s + 1)*s, sizeof(ScalarType), 64);if(R == NULL){return 1;}
	H = (ScalarType *)mkl_calloc((m + 1)*m, sizeof(ScalarType), 64);if(H == NULL){return 1;}
	H_reduced = (ScalarType *)mkl_calloc((m + 1)*m, sizeof(ScalarType), 64);if(H_reduced == NULL){return 1;}
	r = (ScalarType *)mkl_malloc(n * sizeof(ScalarType), 64);if(r == NULL){return 1;}
	x = (ScalarType *)mkl_malloc(n * sizeof(ScalarType), 64);if(x == NULL){return 1;}
	tx = (ScalarType *)mkl_malloc(n * sizeof(ScalarType), 64);if(tx == NULL){return 1;}
	zeta = (ScalarType *)mkl_calloc(m + 1, sizeof(ScalarType), 64);if(zeta == NULL){return 1;}
	R_k = (ScalarType *)mkl_calloc(s*s, sizeof(ScalarType), 64);if(R_k == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	R_k[0] = 1;

	std::vector<std::pair<ScalarType, ScalarType>> cs(m+1); // cs stores givens rotations for reduction of H

	/* set true solution */
	for (i = 0; i < n; ++i)
		tx[i] = 1;

	/********************************************/
	/*  compute RHS b                           */
	/*  compute residual vector (r = b - Ax_0)  */
	/*  therefore, b == r                       */
	/********************************************/

	struct matrix_descr descr;

	descr.type = SPARSE_MATRIX_TYPE_GENERAL;

	mkl_sparse_set_mv_hint(A, SPARSE_OPERATION_NON_TRANSPOSE, descr, 1);
	mkl_sparse_optimize(A);	
	
	kernels<ScalarType>::mv(A, tx, r, 1);

	ScalarType beta = cblas_dnrm2(n, r, 1);

	zeta[0] = beta;

	cblas_daxpy(n, 1 / beta, r, 1, Q, 1);

	stat = kernels<ScalarType>::gmres_init(n, A, H, H_reduced, Q, theta_vals, ritz_num, m);

	// after init Q contains s+1 orthonormal basis vectors for the Krylov subspace

	std::cout << std:: endl << "theta_vals (FINAL):" << std::endl;
	for (auto p: theta_vals) {
		std::cout << p.second << "   \toutlist: " << p.first << std:: endl;
	}

	// printf("\n\n============= H:\n");
	// for(size_t i = 0; i < ritz_num + 1; ++i) {
		// for(size_t j = 0; j < ritz_num; ++j) {
			// printf("%2.2f ", H[i*m + j]);
		// }
		// std::cout << std::endl;
	// }	
	// std::cout << std::endl;

	/**************/
  /*  reduce H  */
	/**************/

	stat = kernels<ScalarType>::reduce_H(H_reduced, s, m, 0, zeta, cs);	

	// after H_reduced is reduced, zeta contains s+1 values

	// std::cout << "\n\n============= zeta:\n";
	// for(size_t i = 0; i < m + 1; ++i)
		// std::cout << zeta[i] << std::endl;
	// std::cout << std::endl;

	// std::cout << "\n\n============= H reduced:\n";
	// for(size_t i = 0; i < m+1; ++i) {
		// for(size_t j = 0; j < m; ++j) {
			// printf("%2.2f ", H_reduced[(m+1)*j + i]);
		// }
		// std::cout << std::endl;
	// }
	// std::cout << std::endl;	

	do{

		// outer iterations before restart
		for (k = 1; k < t; ++k) {

			// compute A*q_(s+1) and store result in V[:,0]

			mkl_sparse_set_mv_hint(A, SPARSE_OPERATION_NON_TRANSPOSE, descr, s);
			mkl_sparse_optimize(A);				
			
			kernels<ScalarType>::mv(A, &Q[n*(s*k)], V, s);
			cblas_daxpy(n, -theta_vals.at(0).second.real(), &Q[n*(s*k)], 1, V, 1);		

			for (i = 1; i < s; ++i) {

				imag = theta_vals.at(i).second.imag();

				// leave in for debug, matrix may not be sparse enough
				// stat = mkl_sparse_d_create_csr (&A, SPARSE_INDEX_BASE_ZERO, n, n, minfo.rows_start, minfo.rows_end, minfo.col_indx, minfo.values);

				kernels<ScalarType>::mv(A, &V[n*(i - 1)], &V[n*i], s);

				cblas_daxpy(n, -theta_vals.at(i).second.real(), &V[n*(i - 1)], 1, &V[n*i], 1);

				if (imag < 0) {
					cblas_daxpy(n, imag*imag, &V[n*(i - 2)], 1, &V[n*i], 1);
				}
			}

			/**************************/
			/*  "BLOCK" GRAM SCHMIDT  */
			/**************************/

			// void cblas_dgemm (const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa,
									// const CBLAS_TRANSPOSE transb, const MKL_INT m, 
                  // const MKL_INT n, const MKL_INT k,
									// const double alpha, const double *a,
									// const MKL_INT lda, const double *b,
									// const MKL_INT ldb, const double beta,
									// double *c, const MKL_INT ldc);

			cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, (s*k) + 1, s, n, 1, Q, n, V, n, 0, R, m + s + 1);
			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, s, (s*k) + 1, -1, Q, n, R, m + s + 1, 1, V, n);

			/**********/
			/*  TSQR  */
			/**********/
			
			kernels<ScalarType>::tsqr(V, &Q[n*(s*k + 1)], &R[s*k + 1], n, s, m);

			// if (n < 15) {
				// printf("\n============= final Q (col major):\n");
				// for(size_t o = 0; o < n; ++o) {
					// for(size_t j = 0; j < m + 1; ++j) {
						// std::cout << Q[n*j + o] << " ";
					// }
					// std::cout << std::endl;
				// }
			// }		
			
			stat = kernels<ScalarType>::update_H(H, H_reduced, R, R_k, theta_vals, s, m, k);

			stat = kernels<ScalarType>::reduce_H(H_reduced, s, m, k, zeta, cs);

			std::cout << "\n============= rel. res.:\n";
			printf("%e\n", std::abs(zeta[s*k + 1]) / beta);
			
		} // end for "outer iteration"

		restart = false;

	} while(restart);

	/**************************************/
	/*  solve the system (standard GMRES) */
	/**************************************/

	// std::cout << "\n\n============= H final\n";
	// for (size_t i = 0; i < m + 1; ++i) {
		// for (size_t j = 0; j < m; ++j) {
			// printf("%e ", H[(m + 1)*j + i]);
		// }
		// std::cout << std::endl;
	// }
	// std::cout << std::endl;	

	// std::cout << "\n\n============= zeta:\n";
	// for(i = 0; i < m + 1; ++i)
		// std::cout << zeta[i] << std::endl;

// void cblas_dtrsv (const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const CBLAS_DIAG diag, 
                  // const MKL_INT n, const double *a, const MKL_INT lda, double *x, const MKL_INT incx);
	cblas_dtrsv (CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, s*k, H_reduced, m + 1, zeta, 1);

// void cblas_dgemv (const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE trans, const MKL_INT m, const MKL_INT n, const double alpha,
                  // const double *a, const MKL_INT lda, const double *x, const MKL_INT incx, const double beta, double *y, const MKL_INT incy);	
	cblas_dgemv (CblasColMajor, CblasNoTrans, n, s*k, 1, Q, n, zeta, 1, 0, x, 1);

	if (n < 15) {
		
		std::cout.precision(20);
		
		std::cout << "\n\n============= sol x:\n";
		for(i = 0; i < n; ++i)
			std::cout << x[i] << std::endl;

		std::cout << std::endl;
	} else {
		std::cout.precision(20);

		std::cout << "\n\n============= sol x:\n";
		for(i = n - 30; i < n; ++i)
			std::cout << x[i] << std::endl;

		std::cout << std::endl;
	}

	stat = mkl_sparse_destroy(A);

	mkl_free(V);
	mkl_free(Q);
	mkl_free(R);
	mkl_free(H);
	mkl_free(H_reduced);
	mkl_free(R_k);
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
// mkl_dimatcopy('R', 'T', s + 1, n, 1, Q, n, s + 1);

/*****************************************/
/*  cout 'V', 'Q' and result matrix 'R'  */
/*****************************************/

			// if (n < 15) {
				// printf("\n============= final Q (col major):\n");
				// for(size_t k = 0; k < n; ++k) {
					// for(size_t j = 0; j < m + 1; ++j) {
						// std::cout << Q[n*j + k] << " ";
					// }
					// std::cout << std::endl;
				// }
			// }

			// if (n < 15) {
				// printf("\n============= V (col major):\n");
				// for(size_t k = 0; k < n; ++k) {
					// for(size_t j = 0; j < s; ++j) {
						// std::cout << V[n*j + k] << " ";
					// }
					// std::cout << std::endl;
				// }
			// }

			// printf("\n============= R (col major):\n");
			// for(size_t o = 0; o < m + s + 1; ++o) {
				// for(size_t j = 0; j < s; ++j) {
					// std::cout << R[(m + s + 1)*j + o] << " ";
				// }
				// std::cout << std::endl;
			// }

			// printf("\n============= R_ (col major):\n");
			// for(size_t o = 0; o < s; ++o) {
				// for(size_t j = 0; j < s; ++j) {
					// std::cout << R[(m + s + 1)*j + (o + (s*k + 1))] << " ";
				// }
				// std::cout << std::endl;
			// }
		
			// printf("\n============= Q row major:\n");
			// for(size_t j = 0; j < n*(s + 1); ++j) {
				// std::cout << Q[j] << ", ";
			// }
			// std::cout << std::endl;