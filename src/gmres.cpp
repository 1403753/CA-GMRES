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
	
	gsl_matrix *A = gsl_matrix_alloc(M, N);
	gsl_vector *u = gsl_vector_alloc(N);
	
	for (size_t i = 0; i < M; ++i)
		for (size_t j = 0; j < N; ++j)
      gsl_matrix_set(A, i, j, 1);//gsl_rng_uniform_int(rng, 100) + gsl_rng_uniform(rng));
	
	for (size_t j = 0; j < N; ++j)
      gsl_vector_set(u, j, 1);

	gsl_vector *y = gsl_vector_alloc(A->size1);
	
	
	// std::vector<double> A_(M*N);
	// std::vector<double> C_(M*N);

	double u_[N];
	double **C_ptr, *C_;
	C_ptr= (double **)mkl_malloc(M * sizeof(double *), 64); //check malloc!
	C_= (double *)mkl_calloc(M*N, sizeof(double), 64); //check malloc!
	
	size_t idx = 0;
	for(size_t i = 0; i < M*N; i += M) {
		C_ptr[idx++] = &C_[i];
	}
	
	sparse_matrix_t A_;
	
	MKL_INT rows_ptr[] = {0, 1, 2, 3, 4, 5};
	MKL_INT col_indx[] = {0, 0, 1, 0, 0};
	double values[] = 	 {3, 2, 4, 9, 10};
		
	sparse_status_t stat;	
	stat = mkl_sparse_d_create_csr(&A_, SPARSE_INDEX_BASE_ZERO, M, N, rows_ptr, rows_ptr + 1, col_indx, values);
	if(stat != SPARSE_STATUS_SUCCESS)
		throw std::invalid_argument("MatCreate failed");
		
	int k = 4;
	for(auto &v: u_)
		v = k--;
	
	printf("\n============= C_:\n");
	for(size_t i = 0; i < M*N; ++i) {
		printf("%f, ", C_[i]);
	}
	std::cout << std::endl;
	
	for(size_t i = 0; i < N; ++i) {
		std::cout << "U: ";
		for(auto u: u_)
			std::cout << u << " ";
		std::cout << std::endl;
		
		spmv::mv(A_, u_, &C_ptr[i]);
		/*
		for(size_t j = 0; j < M; ++j) {
			C_[i*M+j] = k++;
		}
		k+=100;
		*/
		printf("\n============= C_:\n");
		for(size_t i = 0; i < M*N; ++i) {
			printf("%f, ", C_[i]);
		}
		std::cout << std::endl;
		
		
		
		for(size_t j = 0; j < N; ++j)
			u_[j] = C_[j+i];
	}
	
	/*
		transpose matrix in place
	*/
	// void mkl_dimatcopy (const char ordering, const char trans, size_t rows, size_t cols, const double alpha, double * AB, size_t lda, size_t ldb);
	// mkl_dimatcopy('R', 'T', N, M, 1, C_, M, N);
	
	printf("\n============= C_:\n");
	for(size_t i = 0; i < M*N; ++i) {
		printf("%f, ", C_[i]);
	}
	std::cout << std::endl;

	stat = mkl_sparse_destroy(A_);
		



//	mpk::mv(A, u, y);
	tsqr::qr(&C_, M, N);
//	arnoldi_ca::print();
	
	
//	for (size_t i = 0; i < M; ++i)
//		std::cout << gsl_vector_get(y, i) << ' ';
	/*
	std::cout << std::endl << "C_ptr: ";
	
	for (size_t i = 0; i < M; ++i)
		std::cout << C_ptr[i] << ' ';
	
	std::cout << std::endl << std::endl;
	*/
	gsl_rng_free(rng);
	gsl_matrix_free(A);
	gsl_vector_free(u);
	gsl_vector_free(y);
	
	
	return 0;
}
