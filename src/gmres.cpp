//============================================================================
// Name        : ca-gmres.cpp
// Author      : Robert Ernstbrunner
// Version     :
// Copyright   : Your copyright notice
// Description : ca-gmres in C++
//============================================================================

#include "gmres.hpp"

#define M 5
#define N 3

int main() {
	std::cout.setf(std::ios_base::fixed);
	std::cout.precision(3);
	
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(rng, time(NULL));
	
	// double r[M] = {3,2,4,9,10};
	double r[M] = {1,1,1,1,1};
	double *V, **V_col_ptr;
	
	V_col_ptr = (double **)mkl_malloc(M * sizeof(double *), 64);if(V_col_ptr == NULL){return 1;}
	V = (double *)mkl_calloc(M*N, sizeof(double), 64);if(V == NULL){return 1;}
	
	size_t idx = 0;
	for(size_t i = 0; i < M*N; i += M) {
		V_col_ptr[idx++] = &V[i];
	}
	
	sparse_matrix_t A;
	
	size_t rows_ptr[] = {0, 1, 2, 3, 4, 5};
	size_t col_indx[] = {0, 0, 1, 0, 0};
	double values[] = 	 {3, 2, 4, 9, 10};	
/*  5 x 3 Matrix testcode:	
		MKL_INT rows_ptr[] = {0, 3, 5, 8, 11, 13};
		MKL_INT col_indx[] = {0, 1, 3, 0, 1, 2, 3, 4, 0, 2, 3, 1, 4};
		double values[] = 	 {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
*/
	
		
	sparse_status_t stat;	
	stat = mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO, M, N, rows_ptr, rows_ptr + 1, col_indx, values);
	if(stat != SPARSE_STATUS_SUCCESS)
		throw std::invalid_argument("MatCreate failed");
	
	printf("\n============= V:\n");
	for(size_t i = 0; i < M*N; ++i) {
		printf("%f, ", V[i]);
	}
	std::cout << std::endl;
	
	cblas_dcopy (M, r, 1, V_col_ptr[0], 1);
	
	for(size_t i = 0; i < N - 1; ++i) {
		
		spmv::mv(A, V_col_ptr[i], &V_col_ptr[i+1]);

		printf("\n============= V:\n");
		for(size_t j = 0; j < M; ++j) {
			for(size_t k = 0; k < N; ++k) {
				std::cout << V[k*M + j] << '\t';
			}
			std::cout << std::endl;
		}
	}
	
	/*
		transpose matrix in place
	*/
	// void mkl_dimatcopy (const char ordering, const char trans, size_t rows, size_t cols, const double alpha, double * AB, size_t lda, size_t ldb);
	// mkl_dimatcopy('R', 'T', N, M, 1, V, M, N);



//	mpk::mv(A, u, y);
	tsqr::qr(&V, M, N);
	arnoldi_ca::givens_rotations();
	arnoldi_ca::modified_leja_ordering();
	
	
	
	gsl_rng_free(rng);
	stat = mkl_sparse_destroy(A);
	mkl_free(V);
	mkl_free(V_col_ptr);
	
	return 0;
}
