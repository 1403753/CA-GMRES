/*
 * GMRES_ca.cpp
 *
 *  Created on: 17.10.2018
 *      Author: Robert
 */

#include "GMRES_ca.hpp"

GMRES_ca::GMRES_ca() {
	
}

GMRES_ca::~GMRES_ca() {

}

sparse_status_t GMRES_ca::mpk(double *x, double *dest) {
	
	ILU0_ca *prec = (ILU0_ca*) ksp->getIPCType();
	std::vector<std::vector<Mtx_CSR>> *A_mtxArr = prec->getA_mtxArr();
	std::vector<std::vector<Mtx_CSR>> *L_mtxArr = prec->getL_mtxArr();
	std::vector<std::vector<Mtx_CSR>> *U_mtxArr = prec->getU_mtxArr();
	std::vector<std::vector<size_t>> *alphaArr = prec->getAlphaArr();
	size_t *Ap;
	size_t *Ai;
	double *Ax;
	size_t n, nz;
	size_t nParts = prec->getNParts();
	size_t globalN = ksp->getA_mtx()->n;
	// #pragma omp parallel for
	for (size_t i = 0; i < nParts; ++i) {
		for (size_t j = 0; j < 1; ++j) {
			n = A_mtxArr->at(i).at(j).n;
			size_t length = n / nParts;
			// size_t offset = length * omp_get_thread_num();
			nz = A_mtxArr->at(i).at(j).nz;
			Ap = A_mtxArr->at(i).at(j).row_ptr;
			Ai = A_mtxArr->at(i).at(j).col_indx;
			Ax = A_mtxArr->at(i).at(j).values;
			alphaArr.at(i);
			for (size_t i = 0; i < n; ++i) {
				double rowsum = 0.0;
				for (size_t j = Ap[i]; j < Ap[i+1]; ++j) {
					// if (x + Ai[j] + offset - x < n)
						// rowsum += Ax[j] * *(x + Ai[j] + offset);
						rowsum += Ax[j] * x[Ai[j]];
				}
				
				x = rowsum;
			}
			// mkl_dcsrtrsv ('L', 'N', 'N', &A_mtxArr->at(i).at(j).n, A_mtxArr->at(i).at(j).ilu0_values, A_mtxArr->at(i).at(j).row_ptr, A_mtxArr->at(i).at(j).col_indx, x, dest);
		}
	}
		
	for (size_t i = 0; i < globalN; ++i)
		std::cout << dest[i] << " .. " << x[i] << "\n";
	
	
	return SPARSE_STATUS_SUCCESS;
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
	
	////////////////////////////////
	//  LAPACKE_dhseqr variables  //
	////////////////////////////////

	size_t                                  ilo = 1;       // if A hasn't been balanced by ?gebal: ilo = 1
	size_t                                  ihi = s;       // if A hasn't been balanced by ?gebal: ihi = s
	char                                    job = 'E';     // eigenvalues only are required
	const char                              compz = 'N';   // no Schur vectors are computed
	double                                  *wr, *wi;      // real and imag part of ritz values
	// double																   *scale;
	size_t                                  i, j;          // index in for-loops

	w = (double *)mkl_calloc(n, sizeof(double), 64);if(w == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	H_s = (double *)mkl_calloc((s + 1)*s, sizeof(double), 64);if(H_s == NULL){return SPARSE_STATUS_ALLOC_FAILED;}

	////////////////////////////
	//  Modfied Gram-Schmidt  //
	////////////////////////////

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
	
	double beta = cblas_dnrm2(n, b, 1);

	zeta[0] = beta;

	cblas_daxpy(n, 1 / beta, b, 1, Q, 1);

	gmres_init(n, *A_mkl, H, H_reduced, Q, theta_vals, s, m);

	// after init Q contains s+1 orthonormal basis vectors for the Krylov subspace

	////////////////
  //  reduce H  //
	////////////////

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

			////////////////////////////
			//  "BLOCK" GRAM SCHMIDT  //
			////////////////////////////

			cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, (s*k) + 1, s, n, 1, Q, n, V, n, 0, R, m + s + 1);
			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, s, (s*k) + 1, -1, Q, n, R, m + s + 1, 1, V, n);

			////////////
			//  TSQR  //
			////////////
			
			tsqr(V, &Q[n*(s*k + 1)], &R[s*k + 1], n, s, m);
			
			update_H(H, H_reduced, R, R_k, theta_vals, s, m, k);
			reduce_H(H_reduced, s, m, k, zeta, cs);

			std::cout << "\n============= rel. res.:\n";
			printf("%e\n", std::abs(zeta[s*k + 1]) / beta);
			
		} // end for "outer iteration"

		restart = false;

	} while(restart);

	
	/////////////////////////////////////////
	//  solve the system (standard GMRES)  //
	/////////////////////////////////////////
	
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
	
	// ILU0_ca *prec = (ILU0_ca*) ksp->getIPCType();
	
	// prec->destroy_MtxArrs();

	return SPARSE_STATUS_SUCCESS;
}

bool GMRES_ca::is_conj_pair(complex_t a, complex_t b) {
	return (a.real() == b.real() && a.imag() == -b.imag() && a.imag() != 0);
}

sparse_status_t GMRES_ca::modified_leya_ordering(size_t s, double *wr, double *wi, std::vector<ic_pair_t> &theta_vals) {
	std::vector<ic_pair_t> ritz_vals;

//	ritz_vals.reserve(sizeof(wr) / sizeof(wr[0]));
	ritz_vals.reserve(s);
	
	for(size_t i = 0; i < s; ++i)
		ritz_vals.emplace_back(ic_pair_t(1, complex_t(wr[i], wi[i])));		
		
	std::stable_sort(ritz_vals.begin( ), ritz_vals.end( ), [ ]( const ic_pair_t &lhs, const ic_pair_t &rhs ) {
			return lhs.second.real() < rhs.second.real();
	});
	
	for(size_t i = 0; i < ritz_vals.size() - 1;) {
		if(ritz_vals.at(i).second.imag() > 0) {
			if(i + 3 < ritz_vals.size() && ritz_vals.at(i).second == ritz_vals.at(i + 2).second && ritz_vals.at(i + 1).second == ritz_vals.at(i + 3).second) {
				ritz_vals.at(i).first++;
				ritz_vals.erase(ritz_vals.begin() + i + 2);
				ritz_vals.at(i + 1).first++;
				ritz_vals.erase(ritz_vals.begin() + i + 2);
			} else {
				i += 2;
			}
		} else {
			if(ritz_vals.at(i).second == ritz_vals.at(i + 1).second) {
				ritz_vals.at(i).first++;
				ritz_vals.erase(ritz_vals.begin() + i + 1);
			} else {
				i++;
			}
		}
	}

	ritz_vals.shrink_to_fit();
	
	std::vector<size_t> k_index;

	theta_vals.reserve(ritz_vals.size());
	k_index.reserve(ritz_vals.size());
	
	for(size_t i = 0; i < ritz_vals.size(); ++i)
		k_index.push_back(i);

	/*
	*
	* MODIFIED LEJA ORDERING:
	*
	*	input: ritz_vals (unique shifts with multiplicity mu)
	*				 type: std::vector<size_t mu, complex<double> vals> 
	*
	*	output: theta_vals (in modified Leja order with outlist)
	*					type: std::vector<size_t outlist, complex<double> vals>
	*/
	
	complex_t Capacity_old, Capacity = 1;
	size_t L;


	if(std::abs(ritz_vals.back().second) < std::abs(ritz_vals.front().second)) {
		
		if (ritz_vals.front().second.imag() == 0) {

			theta_vals.push_back(ic_pair_t(std::move(k_index.front()), ritz_vals.front().second));
			k_index.erase(k_index.begin());
			L = 0;
		} else {
			
			theta_vals.push_back(ic_pair_t(std::move(k_index.front()), ritz_vals.front().second));
			theta_vals.push_back(ic_pair_t(std::move(k_index.begin()[1]), ritz_vals.begin()[1].second));

			k_index.erase(k_index.begin());
			k_index.erase(k_index.begin());
			L = 1;
		}
	} else {

		if (ritz_vals.back().second.imag() == 0) {

			theta_vals.push_back(ic_pair_t(std::move(k_index.back()), ritz_vals.back().second));
			k_index.pop_back();
			L = 0;
		} else {
			
			theta_vals.push_back(ic_pair_t(std::move(k_index.end()[-2]), ritz_vals.end()[-2].second));
			theta_vals.push_back(ic_pair_t(std::move(k_index.back()), ritz_vals.back().second));
			
			k_index.pop_back();
			k_index.pop_back();				
			L = 1;
		}
	}
	
	
	while(L < ritz_vals.size() - 1) {

		Capacity_old = Capacity;
		Capacity = 1;
		
		for(size_t j = 0; j < L; ++j) {
			Capacity *= std::pow(std::abs(theta_vals.at(L).second - theta_vals.at(j).second), (double) ritz_vals.at(theta_vals.at(j).first).first / (L+1));
		}
				
		for (auto &z: ritz_vals) {
			z.second /= (Capacity / Capacity_old);
		}
		
		for (auto &t: theta_vals)
			t.second /= (Capacity / Capacity_old);
		
		std::vector<complex_t> zprod;
		
		for (auto k: k_index) {
			complex_t prod = 1;
			for(size_t j = 0; j < L+1; ++j) {
				prod *= std::pow(std::abs(ritz_vals.at(k).second - ritz_vals.at(theta_vals.at(j).first).second) / Capacity, ritz_vals.at(theta_vals.at(j).first).first);
			}
			zprod.push_back(prod);
		}
		
		std::vector<complex_t>::iterator max_zprod;
		max_zprod = std::max_element(zprod.begin( ), zprod.end( ), [ ]( const complex_t& lhs, const complex_t& rhs ) {
			return std::abs(lhs) < std::abs(rhs);
		});
		
		if (*max_zprod == (complex_t)0 ) {
			throw std::invalid_argument("Product to maximize is zero; either there are multiple shifts, or the product underflowed");
		} else if (std::isinf((*max_zprod).real()) ||	std::isnan((*max_zprod).real()) )
			throw std::invalid_argument("Product to maximize is Inf; must have overflowed");
		
		size_t idx = max_zprod - zprod.begin();
		

		
		if(ritz_vals.at(k_index.at(idx)).second.imag() != 0) {
			
			if(ritz_vals.at(k_index.at(idx)).second.imag() < 0) {
				if (!is_conj_pair(ritz_vals.at(k_index.at(idx) - 1).second, ritz_vals.at(k_index.at(idx)).second)) 
					throw std::invalid_argument( "Input out of order");
				
				theta_vals.push_back(ic_pair_t(std::move(k_index.at(idx - 1)), ritz_vals.at(k_index.at(idx) - 1).second));
				theta_vals.push_back(ic_pair_t(std::move(k_index.at(idx)), ritz_vals.at(k_index.at(idx)).second));			

				k_index.erase(k_index.begin() + idx - 1);
				k_index.erase(k_index.begin() + idx - 1);
				
				L += 2;	
			
			
			} else {
				if (!is_conj_pair(ritz_vals.at(k_index.at(idx)).second, ritz_vals.at(k_index.at(idx + 1)).second))
					throw std::invalid_argument( "Input out of order");
				
				theta_vals.push_back(ic_pair_t(std::move(k_index.at(idx)), ritz_vals.at(k_index.at(idx)).second));
				theta_vals.push_back(ic_pair_t(std::move(k_index.at(idx + 1)), ritz_vals.at(k_index.at(idx + 1)).second));			
				
				k_index.erase(k_index.begin() + idx);
				k_index.erase(k_index.begin() + idx);

				L += 2;
			}
		} else {
		
			theta_vals.push_back(ic_pair_t(k_index.at(idx), ritz_vals.at(k_index.at(idx)).second));
			k_index.erase(k_index.begin() + idx);
		
			L++;
		}
	} //end while
	
	// std::cout << std:: endl << "ritz_vals :" << std:: endl;
	for (auto &p: ritz_vals) {
		p.second *= Capacity;
		// std::cout << p.second << "   \tmult: " << p.first << std::endl;
	}
	
	// std::cout << std:: endl << "theta_vals (FINAL):" << std:: endl;
	for (auto &p: theta_vals) {
		p.second *= Capacity;
		// std::cout << p.second << "   \toutlist: " << p.first << std:: endl;
	}
	
	return SPARSE_STATUS_SUCCESS;
}

sparse_status_t GMRES_ca::tsqr(double *V, double *Q, double *R_, const size_t m, const size_t n, const size_t st) {

	size_t lwork = -1;
	size_t lmwork = -1;
	size_t tsize = -1;
	size_t info = 0;
	size_t i, j;
	double *work, *t, *tquery, *workquery, *mwork, *mworkquery;
	char side = 'L', trans = 'N';
	sparse_status_t stat = SPARSE_STATUS_SUCCESS;
		
	for(i = 0; i < n; ++i) 
		Q[i*m + i] = 1;

	tquery = (double *)mkl_malloc(5 * sizeof(double), 64);if (tquery == NULL)throw std::invalid_argument("malloc error: allocating memory for tquery failed.");
	workquery = (double *)mkl_malloc(sizeof(double), 64);if (workquery == NULL)throw std::invalid_argument("malloc error: allocating memory for workquery failed.");
	
	/* Workspace query */
	dgeqr(&m, &n, V, &m, tquery, &tsize, workquery, &lwork, &info);
	
	tsize = tquery[0];
	mkl_free(tquery);
	
	t = (double *)mkl_malloc(tsize * sizeof(double), 64);if (t == NULL)throw std::invalid_argument("malloc error: allocating memory for t failed.");
	
	lwork = workquery[0];
	mkl_free(workquery);	
	
	work = (double *)mkl_malloc(lwork * sizeof(double), 64);if (work == NULL)throw std::invalid_argument("malloc error: allocating memory for work failed.");
	
	dgeqr(&m, &n, V, &m, t, &tsize, work, &lwork, &info);
	mkl_free(work);
	
	mworkquery = (double *)mkl_malloc(sizeof(double), 64);if (mworkquery == NULL)throw std::invalid_argument("malloc error: allocating memory for mworkquery failed.");
	
	/************************/
	/*  store Q explicitly  */
	/************************/
	
	/* Workspace query */
	dgemqr(&side, &trans, &m, &n, &n, V, &m, t, &tsize, Q, &m, mworkquery, &lmwork, &info);

	lmwork = mworkquery[0];
	mkl_free(mworkquery);
		
	mwork = (double *)mkl_malloc(lmwork * sizeof(double), 64);if (mwork == NULL)throw std::invalid_argument("malloc error: allocating memory for mwork failed.");
	
	dgemqr(&side, &trans, &m, &n, &n, V, &m, t, &tsize, Q, &m, mwork, &lmwork, &info);
	mkl_free(mwork);

	mkl_free(t);
	
	for(i = 0; i < n; ++i) {
		for(j = i; j < n; ++j) {
			R_[(st + n + 1)*j + i] = V[m*j + i];
		}
	}
	
	if (info < 0) stat = SPARSE_STATUS_EXECUTION_FAILED;

	return stat;
}