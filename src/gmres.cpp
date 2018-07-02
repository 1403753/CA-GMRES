//============================================================================
// Name        : ca-gmres.cpp
// Author      : Robert Ernstbrunner
// Version     :
// Copyright   : Your copyright notice
// Description : ca-gmres in C++
//============================================================================

#include "gmres.hpp"

#define s 500

int main() {
	std::cout.setf(std::ios_base::fixed);
	std::cout.precision(3);
	
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(rng, time(NULL));
	
	float							rtime, ptime, mflops;
	long long 				flpops;
	double 						*V, **V_col_ptr, *r, beta;
	sparse_status_t 	stat;	
	sparse_matrix_t 	A;
	MatrixInfo				minfo;
	size_t						indx;

	// stat = matrix_reader::read_matrix_from_file("../matrix_market/mini_test.mtx", &A, &minfo);
	stat = matrix_reader::read_matrix_from_file("../matrix_market/bmw7st_1.mtx", &A, &minfo);
	if(stat != SPARSE_STATUS_SUCCESS) throw std::invalid_argument("MatCreate failed");

	const size_t m = minfo.rows;

	r = (double *)mkl_calloc(m, sizeof(double), 64);if(r == NULL){return 1;}
	V_col_ptr = (double **)mkl_malloc(m * sizeof(double *), 64);if(V_col_ptr == NULL){return 1;}
	V = (double *)mkl_calloc(m * (s+1), sizeof(double), 64);if(V == NULL){return 1;}

	for(size_t i = 0; i < m; ++i)
		r[i] = 1;
	
	indx = 0;
	for(size_t i = 0; i < m * (s+1); i += m)
		V_col_ptr[indx++] = &V[i];

/*
	compute L2-norm beta
*/
// double cblas_dnrm2 (const MKL_INT n, const double *x, const MKL_INT incx);
	beta = cblas_dnrm2 (m, r, 1);
/*
	copy and scale residual vector into V matrix
*/
	cblas_dcopy (m, r, 1, V_col_ptr[0], 1);
// void	cblas_dscal (const int N, const double alpha, double *X, const int incX)
	cblas_dscal (m, 1 / beta, V_col_ptr[0], 1);
	
	if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
			exit(1);
	
	for(size_t i = 0; i < s; ++i) {

		spmv::mv (A, V_col_ptr[i], &V_col_ptr[i+1], s);
		
		// beta = cblas_dnrm2 (m, V_col_ptr[i + 1], 1);
		// cblas_dscal (m, 1 / beta, V_col_ptr[i + 1], 1);
		
		// printf("\n============= V:\n");
		// for(size_t j = 0; j < m; ++j) {
			// for(size_t k = 0; k < s+1; ++k) {
				// std::cout << V[k*m + j] << '\t';
			// }
			// std::cout << std::endl;
		// }
	}
	
	if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
			exit(1);
	
	PAPI_shutdown();	

	// printf("\n============= V col major:\n");
	// for(size_t j = 0; j < m*(s+1); ++j) {
		// std::cout << V[j] << ", ";
	// }
	// std::cout << std::endl;

	
	/*
		transpose matrix in place
	*/
	// void mkl_dimatcopy (const char ordering, const char trans, size_t rows, size_t cols, const double alpha, double * AB, size_t lda, size_t ldb);
	// mkl_dimatcopy('R', 'T', s+1, m, 1, V, m, s+1);

	// printf("\n============= V row major:\n");
	// for(size_t j = 0; j < m*(s+1); ++j) {
		// std::cout << V[j] << ", ";
	// }
	// std::cout << std::endl;

	// if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
			// exit(1);
		
	tsqr::qr(&V, m, s+1);
	
	// if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
		// exit(1);
	
	// PAPI_shutdown();	
	
	arnoldi_ca::givens_rotations();
	arnoldi_ca::modified_leja_ordering();

	stat = mkl_sparse_destroy(A);
	gsl_rng_free(rng);
	mkl_free(V);
	mkl_free(r);
	mkl_free(V_col_ptr);
	printf("runtime SpMV: %f\n", rtime);

	return 0;
}

/* in place matrix / vector creation

	size_t rows_ptr[] = {0, 1, 2, 3, 4, 5};
	size_t col_indx[] = {0, 0, 1, 0, 0};
	double values[] = {3, 2, 4, 9, 10};	
	double beta;
	double r[M] = {3,2,4,9,10};
	// double r[M] = {1,1,1,1,1};
  // 5 x 3 Matrix testcode:	
		// MKL_INT rows_ptr[] = {0, 3, 5, 8, 11, 13};
		// MKL_INT col_indx[] = {0, 1, 3, 0, 1, 2, 3, 4, 0, 2, 3, 1, 4};
		// double values[] = 	 {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	
	stat = mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO, M, N, rows_ptr, rows_ptr + 1, col_indx, values);
	if(stat != SPARSE_STATUS_SUCCESS)
		throw std::invalid_argument("MatCreate failed");
	
	V_col_ptr = (double **)mkl_malloc(M * sizeof(double *), 64);if(V_col_ptr == NULL){return 1;}
	V = (double *)mkl_calloc(M*N, sizeof(double), 64);if(V == NULL){return 1;}
	
	size_t idx = 0;
	for(size_t i = 0; i < M*N; i += M) {
		V_col_ptr[idx++] = &V[i];
	}
*/
