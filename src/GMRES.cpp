/*
 * GMRES.cpp
 *
 *  Created on: 17.10.2018
 *      Author: Robert
 */

#include "GMRES.hpp"

GMRES::GMRES() {

}

GMRES::GMRES(size_t m) : m{m} {

}

GMRES::~GMRES() {

}

sparse_status_t GMRES::solve(double *x_0, double *b) {
	this->SpMV = 0;
	this->MGS = 0;
	double                                     *Q;                         // orthonormal basis for Krylov subspace K(A,v), Arnoldi output
	double                                     *H;                         // Hessenberg matrix, Arnoldi output
	double                                     *r;                         // residual vector b - Ax
	double 						                         *zeta;                      // rotated beta*e1
	bool                                       restart = true;             // restart as long as value is true


	sparse_status_t                            stat;
	sparse_matrix_t                            *A_mkl = ksp->getA_mkl();   // n x n matrix A
	const std::shared_ptr<Mtx_CSR>             A_ptr = ksp->getA_ptr();    // n x n matrix A
	size_t                                     n = A_ptr->n;               // dim(A)
	std::vector<std::pair<double, double>>     cs(this->m+1);              // cs stores givens rotations for reduction of H	
	struct matrix_descr                        descr;
	double                                     *x;
	double                                     rRes;
	size_t                                     iter = 0;                   // number of total iterations
	double                                     beta, r_0nrm;               // r-norms
	double                                     *w;
	double                                     h_ij;
	double                                     h_jp1j;
	double                                     x1;
	double                                     x2;
	size_t                                     i, j;          // indices in for-loops
	std::vector<std::pair<size_t, double>>*    rHist;

	Q = (double *)mkl_calloc(n*(m + 1), sizeof(double), 64);if(Q == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	H = (double *)mkl_calloc((m + 1)*m, sizeof(double), 64);if(H == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	r = (double *)mkl_malloc(n * sizeof(double), 64);if(r == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	zeta = (double *)mkl_calloc(m + 1, sizeof(double), 64);if(zeta == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	x = (double *)mkl_calloc(n, sizeof(double), 64);if(x == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	w = (double *)mkl_malloc(n * sizeof(double), 64);if(w == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	
	IPCType *pc = ksp->getIPCType();

	descr.type = SPARSE_MATRIX_TYPE_GENERAL;

	stat = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, *A_mkl, descr, x_0, 0, r);

	for (size_t i = 0; i < n; ++i) {
		r[i] = b[i] - r[i];
	}

	pc->precondition(r);

	r_0nrm = cblas_dnrm2(n, r, 1);

	if (ksp->getStoreHist()) {
		rHist = ksp->getRHist();
		rHist->clear();
		rHist->reserve(ksp->getMaxit());
		if (ksp->getStoreHist()) {
			rHist->push_back(std::pair<size_t, double>(iter, r_0nrm / r_0nrm));
		}
	}
	
	do{

		mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, *A_mkl, descr, x_0, 0, r);	

		for (size_t i = 0; i < n; ++i) {
			r[i] = b[i] - r[i];
		}
	
		pc->precondition(r);

		beta = cblas_dnrm2(n, r, 1);
	
		zeta[0] = beta;
		beta = 1 / beta;
		cblas_daxpy(n, beta, r, 1, Q, 1);

		for(j = 0; j < m; ++j) {

			if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
				exit(1);			

			pc->mv(&Q[n*j], w, descr);

			if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
				exit(1);

			PAPI_shutdown();

			this->SpMV += rtime;

			////////////////////////////
			//  Modfied Gram-Schmidt  //
			////////////////////////////

			if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
				exit(1);

			for(i = 0; i < j + 1; ++i) {

				h_ij = cblas_ddot(n, w, 1, &Q[i*n], 1);
				H[(m + 1)*j + i] = h_ij;
				
				cblas_daxpy (n, -h_ij, &Q[i*n], 1, w, 1);
			}

			h_jp1j = cblas_dnrm2(n, w, 1);
			H[(m + 1)*j + (j + 1)] = h_jp1j;

			h_jp1j = 1/h_jp1j;

			cblas_daxpy(n, h_jp1j, w, 1, &Q[n*(j + 1)], 1);

			if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
				exit(1);

			PAPI_shutdown();

			this->MGS += rtime;

			for (size_t i = 0; i < j; ++i) {
		// void cblas_drot (const size_t n, double *x, const size_t incx, double *y, const size_t incy, const double c, const double s);
				cblas_drot(1, &H[(m + 1)*j + i], m + 1, &H[(m + 1)*j + (i + 1)], m + 1, cs.at(i).first, cs.at(i).second);
			}

			x1 = H[(m + 1)*j + j];
			x2 = H[(m + 1)*j + (j + 1)];

		// void cblas_drotg (double *x1, double *x2, double *c, double *s);
			cblas_drotg(&x1, &x2, &cs.at(j).first, &cs.at(j).second);

		// void cblas_drot (const size_t n, double *x, const size_t incx, double *y, const size_t incy, const double c, const double s);
			cblas_drot(1, &H[(m + 1)*j + j], m + 1, &H[(m + 1)*j + (j + 1)], m + 1, cs.at(j).first, cs.at(j).second);

		// void cblas_drot (const size_t n, double *x, const size_t incx, double *y, const size_t incy, const double c, const double s);
			cblas_drot(1, &zeta[j], 1, &zeta[j + 1], 1, cs.at(j).first, cs.at(j).second);	

			rRes = std::abs(zeta[j + 1]) / r_0nrm;

			// std::cout << "rel. res.: " << rRes << ", tol: " << ksp->getRTol() << std::endl;

			++iter;

			if (ksp->getStoreHist()) {
				rHist->push_back(std::pair<size_t, double>(iter, rRes));
			}					
			
			if (rRes < ksp->getRTol() || iter >= ksp->getMaxit() || rRes > ksp->getDTol()) {
				++j;
				restart = false;
				break;
			}
		}

		////////////////////////
		//  solve the system  //
		////////////////////////

		cblas_dtrsv (CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, j, H, m + 1, zeta, 1);
		cblas_dgemv (CblasColMajor, CblasNoTrans, n, m, 1, Q, n, zeta, 1, 0, x, 1);

		for (size_t i = 0; i < n; ++i) {
			x_0[i] = x_0[i] + x[i];
		}

		if(restart) {
			cs.clear();
			cs.resize(m+1);
			std::fill(zeta, zeta+(m+1), 0);
			std::fill(H, H+((m + 1)*m), 0);
	    std::fill(Q, Q+(n*(m + 1)), 0);
		}

	} while(restart);

	std::cout << "gmres(m) iter: " << iter << ", maxit: " << ksp->getMaxit() << ", r_0nrm: " << r_0nrm << ", r_knrm: " << std::abs(zeta[j]) << std::endl;

	mkl_free(w);
	mkl_free(x);
	mkl_free(Q);
	mkl_free(H);
	mkl_free(r);
	mkl_free(zeta);	
	
	return stat;
}