/*
 * GMRES_ca.cpp
 *
 *  Created on: 17.10.2018
 *      Author: Robert
 */

 
#include "GMRES_ca.hpp"
#include "solve.cc"
#include "tsqr.cc"

GMRES_ca::GMRES_ca() {
	
}

sparse_status_t GMRES_ca::mpk(double *x, double *dest, size_t s) {
	ILU0_ca *prec = (ILU0_ca*) ksp->getIPCType();
	prec->test();	
	

	
	return SPARSE_STATUS_SUCCESS;
}


GMRES_ca::~GMRES_ca() {

}

sparse_status_t GMRES_ca::update_H(double *H, double *H_reduced, double *R, double *R_k, std::vector<ic_pair_t>  theta_vals, size_t s, size_t m, size_t k) {
	sparse_status_t stat = SPARSE_STATUS_SUCCESS;
	double      imag;

	for (size_t i = 0; i < s - 1; ++i) {
		cblas_dcopy(s, &R[(m + s + 1)*i + s*k], 1, &R_k[s*(i + 1)], 1);		
	}
	
	// h_underscore_blackletter_(k - 1,k): -h_blackletter_0 * R_(k - 1, k)
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, s*k, s - 1, s*k, -1, H, m + 1, R, m + s + 1, 1, &H[(m + 1)*(s*k + 1)], m + 1);
	
	// h_underscore_blackletter_(k - 1,k): R_(k-1,k)*B_
	cblas_daxpy(s*k, 1, R, 1, &H[(m + 1)*(s*k)], 1);

	for (size_t i = 1; i < s; ++i) {
		cblas_daxpy(s*k, theta_vals.at(i).second.real(), &R[(m + s + 1)*(i - 1)], 1, &H[(m + 1)*(s*k + i)], 1);
		cblas_daxpy(s*k, 1, &R[(m + s + 1)*i], 1, &H[(m + 1)*(s*k + i)], 1);
		imag = theta_vals.at(i).second.imag();
		if (imag < 0)
			cblas_daxpy(s*k, -(imag*imag), &R[(m + s + 1)*(i - 2)], 1, &H[(m + 1)*(s*k + i)], 1);
	}
		
	// H_k: -h0*e1*e_sk^T*R_blackletter_(0,1)
	cblas_daxpy(s - 1, -H[(m + 1)*(s*k - 1) + s*k], &R[s*k - 1], m + s + 1, &H[(m + 1)*(s*k + 1) + s*k], m + 1);	
	
	// H_k: R_k * B_k
	for (size_t i = 0; i < s - 1; ++i) {
		cblas_daxpy(i + 1, theta_vals.at(i).second.real(), &R_k[s*i], 1, &H[(m + 1)*(s*k + i) + s*k], 1);
		cblas_daxpy(i + 2, 1, &R_k[s*(i + 1)], 1, &H[(m + 1)*(s*k + i) + s*k], 1);
		imag = theta_vals.at(i).second.imag();
		if (imag < 0)
			cblas_daxpy(s, -(imag*imag), &R_k[s*(i - 1)], 1, &H[(m + 1)*(s*k + i) + s*k], 1);
	}

	cblas_daxpy(s, theta_vals.at(s - 1).second.real(), &R_k[s*(s - 1)], 1, &H[(m + 1)*(s*k + (s - 1)) + s*k], 1);
	imag = theta_vals.at(s - 1).second.imag();
	if (imag < 0)
		cblas_daxpy(s, -(imag*imag), &R_k[s*(s - 2)], 1, &H[(m + 1)*(s*k + (s - 1)) + s*k], 1);

	// h_k: roh_k * b_k (here b_k is 1)
	H[(m + 1)*(s*(k + 1) - 1) + s*(k + 1)] = R[(m + s + 1)*(s - 1) + s*(k + 1)] * 1;
	
	int info;
	info = LAPACKE_dtrtri (LAPACK_COL_MAJOR, 'U', 'N', s, R_k, s);

	if (info != 0)
		stat = SPARSE_STATUS_EXECUTION_FAILED;

	// h_underscore_blackletter_(k - 1,k): times R^-1
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, s*k, s, s, 1, &H[(m + 1)*(s*k)], m + 1, R_k, s, 0, &H_reduced[(m + 1)*(s*k)], m + 1);

	// H_k: times R_k^-1
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, s, s, s, 1, &H[(m + 1)*(s*k) + s*k], m + 1, R_k, s, 0, &H_reduced[(m + 1)*(s*k) + s*k], m + 1);	

	// H_k: + roh~^-1*b_k*z*es^T
	cblas_daxpy(s, R_k[s*s - 1]*1, &R[(m + s + 1)*(s - 1) + s*k], 1, &H_reduced[(m + 1)*(s*(k + 1) - 1) + s*k], 1);

	// h_k: times roh~^-1
	H[(m + 1)*(s*(k + 1) - 1) + s*(k + 1)] *= R_k[s*s - 1];
	H_reduced[(m + 1)*(s*k + s - 1) + s*(k + 1)] = H[(m + 1)*(s*(k + 1) - 1) + s*(k + 1)];

	for (size_t i = 0; i < s; ++i) {
		cblas_dcopy(s*(k + 1), &H_reduced[(m + 1)*(s*k + i)], 1, &H[(m + 1)*(s*k + i)], 1);
	}

	return stat;
}


sparse_status_t GMRES_ca::reduce_H(double *H, size_t s, size_t m, size_t k, double *zeta, std::vector<std::pair<double, double>> &cs) {
	sparse_status_t stat = SPARSE_STATUS_SUCCESS;
	double x1;
	double x2;

	for (size_t i = 0; i < s*k; ++i) {
// void cblas_drot (const size_t n, double *x, const size_t incx, double *y, const size_t incy, const double c, const double s);
		cblas_drot(s, &H[(m + 1)*s*k + i], m + 1, &H[(m + 1)*s*k + (i + 1)], m + 1, cs.at(i).first, cs.at(i).second);
	}

	for (size_t i = s*k; i < s*(k+1); ++i) {

		x1 = H[(m + 1)*i + i];
		x2 = H[(m + 1)*i + (i + 1)];

// void cblas_drotg (double *x1, double *x2, double *c, double *s);
		cblas_drotg(&x1, &x2, &cs.at(i).first, &cs.at(i).second);

// void cblas_drot (const size_t n, double *x, const size_t incx, double *y, const size_t incy, const double c, const double s);
		cblas_drot(s*(k + 1) - i, &H[(m + 1)*i + i], m + 1, &H[(m + 1)*i + (i + 1)], m + 1, cs.at(i).first, cs.at(i).second);

// void cblas_drot (const size_t n, double *x, const size_t incx, double *y, const size_t incy, const double c, const double s);
		cblas_drot(1, &zeta[i], 1, &zeta[i + 1], 1, cs.at(i).first, cs.at(i).second);
	}

	return stat;
}


sparse_status_t GMRES_ca::gmres_init(size_t n,
																		 const sparse_matrix_t A,
																		 double *H,
																		 double *H_reduced,
																		 double *Q,
																		 std::vector<ic_pair_t> &theta_vals,
																		 size_t s,
																		 size_t m)
{
	sparse_status_t 			                  stat = SPARSE_STATUS_SUCCESS;
	double                                  *w;
	double                                  *H_s;          // is square with odd dimension to ensure at least one real value bc. we initially needed to pick 's' out of '2s' ritzvalues.
	double                                  h_ij;
	double                                  h_jp1j;
	struct matrix_descr                     descr;
	
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	
	/******************************/
	/*  LAPACKE_dhseqr variables  */
	/******************************/

	size_t                                  ilo = 1;       // if A hasn't been balanced by ?gebal: ilo = 1
	size_t                                  ihi = s;       // if A hasn't been balanced by ?gebal: ihi = s
	char                                    job = 'E';     // eigenvalues only are required
	const char                              compz = 'N';   // no Schur vectors are computed
	double                                  *wr, *wi;      // real and imag part of ritz values
	// double																   *scale;
	size_t                                  i, j;          // index in for-loops

	w = (double *)mkl_calloc(n, sizeof(double), 64);if(w == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	H_s = (double *)mkl_calloc((s + 1)*s, sizeof(double), 64);if(H_s == NULL){return SPARSE_STATUS_ALLOC_FAILED;}

	/**************************/
	/*  Modfied Gram-Schmidt  */
	/**************************/

	for(j = 0; j < s; ++j) {

		mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descr, &Q[n*j], 0, w);		
		
		for(i = 0; i < j + 1; ++i) {

			h_ij = cblas_ddot(n, w, 1, &Q[i*n], 1);
			H_s[s*i + j] = h_ij;
			H[(m+1)*j + i] = h_ij;
			H_reduced[(m + 1)*j + i] = h_ij;
			
			cblas_daxpy (n, -h_ij, &Q[i*n], 1, w, 1);
		}

		h_jp1j = cblas_dnrm2(n, w, 1);
		H_s[s*(j + 1) + j] = h_jp1j;
		H[(m + 1)*j + (j + 1)] = h_jp1j;
		H_reduced[(m + 1)*j + (j + 1)] = h_jp1j;

		cblas_daxpy(n, 1/h_jp1j, w, 1, &Q[n*(j + 1)], 1);
	}

	wr = (double *)mkl_malloc(s*sizeof(double), 64);if(wr == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	wi = (double *)mkl_malloc(s*sizeof(double), 64);if(wi == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	// scale = (double *)mkl_malloc((s + 1)*sizeof(double), 64);if(scale == NULL){return SPARSE_STATUS_ALLOC_FAILED;}

// balancing prob. not needed!
	// job = 'N'; // neither scaled nor pivoted
//LAPACKE_dgebal( int matrix_layout, char job, lapack_int n, double* a, lapack_int lda, lapack_int* ilo, lapack_int* ihi, double* scale );
	// LAPACKE_dgebal(LAPACK_ROW_MAJOR, job, s, H_s, s, &ilo, &ihi, scale);	

	// job = 'E';

//LAPACKE_dhseqr(int matrix_layout, char job, char compz, lapack_int n, lapack_int ilo, lapack_int ihi,
							// double *h, lapack_int ldh, double *wr, double *wi, double *z, lapack_int ldz);
	LAPACKE_dhseqr(LAPACK_ROW_MAJOR, job, compz, s, ilo, ihi, H_s, s, wr, wi, nullptr, s);	

	// std::cout << std::endl;
	// for(size_t o = 0; o < s; ++o)
		// printf("%.10f,\t%.10f\n", wr[o], wi[o]);	

	modified_leya_ordering(s, wr, wi, theta_vals);

	// mkl_free(scale);
	mkl_free(w);
	mkl_free(wr);
	mkl_free(wi);
	mkl_free(H_s);
	
	return stat;
}

sparse_status_t GMRES_ca::solve(double *x, double *b) {
	
	double                        *V;                         // basis for Krylov subspace K(A, v) = span{v, (A^1)v, (A^2)v,...,(A^s-1)v}
	double                        *Q;                         // orthonormal basis for Krylov subspace K(A,v), Arnoldi output
	double                        *R;                         // upper triangular matrix containing construction info for H
                                                            // (at each outer iteration only the last s columns of R are needed)
	double                        *H;                         // Hessenberg matrix, Arnoldi output
	double                        *H_reduced;                 // QR factorized H, upper triangular matrix
  double                        *R_k;                       // s x s submatrix needed for updating H
	double                        *r;                         // residual vector b - Ax
	// double                        *x;                         // solution vector x
	// double                        *tx;                        // true solution
	double 						            *zeta;                      // rotated beta*e1
	double                        imag;                       // imaginary part of Ritz value for Newton Basis	
	bool                          restart = true;             // restart as long as value is true

	std::vector<ic_pair_t>        theta_vals;                 // ritz values in modified leja ordering	

	// sparse_status_t               stat;
	sparse_matrix_t               *A_mkl = ksp->getA_mkl();       // n x n matrix A
	Mtx_CSR                       *A_mtx = ksp->getA_mtx();   // n x n matrix A
	// sparse_matrix_t               *M;                         // n x n preconditioned matrix M == ilu0(A)
	size_t                        n = A_mtx->n;               // dim(A)
	const size_t                  s = ksp->getS();            // stepsize number of 'inner iterations'
	const size_t                  t = ksp->getT();            // number of 'outer iterations' before restart
	const size_t                  m = s*t;                    // restart length
	std::vector<std::pair<double, double>> cs(m+1);           // cs stores givens rotations for reduction of H	
	struct matrix_descr           descr;

	V = (double *)mkl_malloc(n*s * sizeof(double), 64);if(V == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	Q = (double *)mkl_calloc(n*(m + 1), sizeof(double), 64);if(Q == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	R = (double *)mkl_calloc((m + s + 1)*s, sizeof(double), 64);if(R == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	H = (double *)mkl_calloc((m + 1)*m, sizeof(double), 64);if(H == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	H_reduced = (double *)mkl_calloc((m + 1)*m, sizeof(double), 64);if(H_reduced == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	r = (double *)mkl_malloc(n * sizeof(double), 64);if(r == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	zeta = (double *)mkl_calloc(m + 1, sizeof(double), 64);if(zeta == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	R_k = (double *)mkl_calloc(s*s, sizeof(double), 64);if(R_k == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	R_k[0] = 1;

	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	
	mkl_sparse_set_mv_hint(*A_mkl, SPARSE_OPERATION_NON_TRANSPOSE, descr, 1);
	mkl_sparse_optimize(*A_mkl);	

	// kernels<double>::mv(A_mkl, tx, r, 1); --> r == b

	// mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A_mkl, descr, b, 0, b);	
	
	double beta = cblas_dnrm2(n, b, 1);

	zeta[0] = beta;

	cblas_daxpy(n, 1 / beta, b, 1, Q, 1);

	gmres_init(n, *A_mkl, H, H_reduced, Q, theta_vals, s, m);

	// after init Q contains s+1 orthonormal basis vectors for the Krylov subspace

	/**************/
  /*  reduce H  */
	/**************/

	reduce_H(H_reduced, s, m, 0, zeta, cs);	

	// after H_reduced is reduced, zeta contains s+1 values

	do{

		// outer iterations before restart
		for (size_t k = 1; k < t; ++k) {

			// compute A_mkl*q_(s+1) and store result in V[:,0]

			mkl_sparse_set_mv_hint(*A_mkl, SPARSE_OPERATION_NON_TRANSPOSE, descr, s);
			mkl_sparse_optimize(*A_mkl);				
		  
			mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, *A_mkl, descr, &Q[n*(s*k)], 0, V);

			cblas_daxpy(n, -theta_vals.at(0).second.real(), &Q[n*(s*k)], 1, V, 1);

			for (size_t i = 1; i < s; ++i) {

				imag = theta_vals.at(i).second.imag();

				// leave in for debug, matrix may not be sparse enough
				// stat = mkl_sparse_d_create_csr (&A_mkl, SPARSE_INDEX_BASE_ZERO, n, n, minfo.rows_start, minfo.rows_end, minfo.col_indx, minfo.values);

				mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, *A_mkl, descr, &V[n*(i - 1)], 0, &V[n*i]);

				cblas_daxpy(n, -theta_vals.at(i).second.real(), &V[n*(i - 1)], 1, &V[n*i], 1);

				if (imag < 0) {
					cblas_daxpy(n, imag*imag, &V[n*(i - 2)], 1, &V[n*i], 1);
				}
			}

			/**************************/
			/*  "BLOCK" GRAM SCHMIDT  */
			/**************************/

			cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, (s*k) + 1, s, n, 1, Q, n, V, n, 0, R, m + s + 1);
			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, s, (s*k) + 1, -1, Q, n, R, m + s + 1, 1, V, n);

			/**********/
			/*  TSQR  */
			/**********/
			
			tsqr(V, &Q[n*(s*k + 1)], &R[s*k + 1], n, s, m);
			
			update_H(H, H_reduced, R, R_k, theta_vals, s, m, k);
			reduce_H(H_reduced, s, m, k, zeta, cs);

			std::cout << "\n============= rel. res.:\n";
			printf("%e\n", std::abs(zeta[s*k + 1]) / beta);
			
		} // end for "outer iteration"

		restart = false;

	} while(restart);

	
	/**************************************/
	/*  solve the system (standard GMRES) */
	/**************************************/
	
	cblas_dtrsv (CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, m, H_reduced, m + 1, zeta, 1);	
	cblas_dgemv (CblasColMajor, CblasNoTrans, n, m, 1, Q, n, zeta, 1, 0, x, 1);	

	mkl_free(V);
	mkl_free(Q);
	mkl_free(R);
	mkl_free(H);
	mkl_free(H_reduced);
	mkl_free(R_k);
	mkl_free(r);
	mkl_free(zeta);	
	
	return SPARSE_STATUS_SUCCESS;
}

