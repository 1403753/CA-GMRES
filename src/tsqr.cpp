/*
 * tsqr.cpp
 *
 *  Created on: 06.05.2018
 *      Author: Robert
 */

#include "tsqr.hpp"

tsqr::tsqr() {
	// TODO Auto-generated constructor stub

}


void tsqr::qr(double **A, size_t M, size_t N){
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set(rng, time(NULL));

	const size_t m = M;
	const size_t n = N;
	// const size_t m = 6;
	// const size_t n = 5;
	size_t lwork = -1;
	size_t lmwork = -1;
	size_t tsize = -1;
	size_t info = 0;
	double *Q, *work, *t, *tquery, *workquery, *mwork, *mworkquery;
	char side = 'L', trans = 'N';
	
	// A = (double *)mkl_calloc(n*m, sizeof(double), 64); //check malloc!
	Q = (double *)mkl_calloc(m*n, sizeof(double), 64); //check malloc!
	for(size_t i = 0; i < n; ++i) 
		Q[i*m + i] = 1;
		
/* row-major
	double A[m*n] = {
										23, 58, 94, 76, 58,
										95, 70, 59, 26, 5,
										33, 49, 23, 72, 15,
										79, 22, 74, 10, 76,
										80, 11, 32, 67, 99,
										96, 98, 72, 40, 70
									};
*/

/* col-major
	double A[m*n] = {
										23,95,33,79,80,
										96,58,70,49,22,
										11,98,94,59,23,
										74,32,72,76,26,
										72,10,67,40,58,
										5,15,76,99,70
									};
*/

	tquery = (double *)mkl_malloc(5 * sizeof(double), 64); //check malloc!
	workquery = (double *)mkl_malloc(sizeof(double), 64); //check malloc!

/*
	for(size_t i = 0; i < m; ++i) 
		for(size_t j = 0; j < n; ++j) {
			A[i*n + j] = gsl_rng_uniform_int(rng, 100);// + gsl_rng_uniform(rng);
			if(i + 1 < n)
				A[(i + 1)*n + j] = gsl_rng_uniform_int(rng, 100);// + gsl_rng_uniform(rng);
		}
*/	

	printf("\nC:\n");
	for(size_t i = 0; i < m; ++i) {
		for(size_t j = 0; j < n; ++j) {
			printf("%f, ", (*A)[j*m + i]);
		}
		printf("\n");
	}
		
	/* Workspace query */
//	dgeqr(m, n, A, lda, t, tsize, work, lwork, info);
	dgeqr(&m, &n, *A, &m, tquery, &tsize, workquery, &lwork, &info);
	
	tsize = tquery[0];
	t = (double *)mkl_malloc(tsize * sizeof(double), 64); //check malloc!

	lwork = workquery[0];
	work = (double *)mkl_malloc(lwork * sizeof(double), 64); //check malloc!

	
	dgeqr(&m, &n, *A, &m, t, &tsize, work, &lwork, &info);
	
	printf("\n============= R:\n");
	for(size_t i = 0; i < m; ++i) {
		for(size_t j = 0; j < n; ++j) {
			printf("%f, ", (*A)[j*m + i]);
		}
		printf("\n");
	}
	
	mworkquery = (double *)mkl_malloc(sizeof(double), 64); //check malloc!

	/* Workspace query */
	dgemqr(&side, &trans, &m, &n, &n, *A, &m, t, &tsize, Q, &m, mworkquery, &lmwork, &info);

	lmwork = mworkquery[0];
	mwork = (double *)mkl_malloc(lmwork * sizeof(double), 64); //check malloc!

	
	
	dgemqr(&side, &trans, &m, &n, &n, *A, &m, t, &tsize, Q, &m, mwork, &lmwork, &info);

	printf("\n============= R:\n");
	for(size_t i = 0; i < m; ++i) {
		for(size_t j = 0; j < n; ++j) {
			printf("%f, ", (*A)[j*m + i]);
		}
		printf("\n");
	}
	printf("\n============= Q:\n");
	for(size_t i = 0; i < m; ++i) {
		for(size_t j = 0; j < n; ++j) {
			printf("%f, ", Q[j*m + i]);
		}
		printf("\n");
	}

	std::cout << "\n\n";	
	for(size_t j = 0; j < n*m; ++j)
		printf("%f, ", (*A)[j]);
	std::cout << "\n";
	
	mkl_free(Q);
	mkl_free(work);
	mkl_free(mwork);
	mkl_free(t);
	mkl_free(mworkquery);
	mkl_free(workquery);
	mkl_free(tquery);
}
tsqr::~tsqr() {
	// TODO Auto-generated destructor stub
	
}