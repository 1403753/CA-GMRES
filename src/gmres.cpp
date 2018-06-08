//============================================================================
// Name        : ca-gmres.cpp
// Author      : Robert Ernstbrunner
// Version     :
// Copyright   : Your copyright notice
// Description : ca-gmres in C++
//============================================================================

#include "gmres.hpp"

#define M 5
#define N 5

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
	
	
	std::vector<double> A_(M*N);
	std::vector<double> u_(N);
	std::vector<double> y_;
	y_.reserve(M);
	
	for(auto &v: A_)
		v = 1;
	
	for(auto &v: u_)
		v = 1;
	
//	spmv::mv(A_, u_, y_);
//	mpk::mv(A, u, y);

	arnoldi_ca::print();
	
	
//	for (size_t i = 0; i < M; ++i)
//		std::cout << gsl_vector_get(y, i) << ' ';
	/*
	std::cout << std::endl << "y_: ";
	
	for (size_t i = 0; i < M; ++i)
		std::cout << y_[i] << ' ';
	
	std::cout << std::endl << std::endl;
	*/
	gsl_rng_free(rng);
	gsl_matrix_free(A);
	gsl_vector_free(u);
	gsl_vector_free(y);
	
	
	return 0;
}
