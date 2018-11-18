template <typename ScalarType>
sparse_status_t kernels<ScalarType>::update_H(ScalarType *H, ScalarType *H_reduced, ScalarType *R, ScalarType *R_k, std::vector<ic_pair_t>  theta_vals, size_t s, size_t m, size_t k) {
	sparse_status_t stat = SPARSE_STATUS_SUCCESS;
	ScalarType      imag;

	for (size_t i = 0; i < s - 1; ++i) {
		cblas_dcopy(s, &R[(m + s + 1)*i + s*k], 1, &R_k[s*(i + 1)], 1);		
	}
	
	// printf("\n============= R_k (col major):\n");
	// for (size_t o = 0; o < s; ++o) {
		// for (size_t j = 0; j < s; ++j) {
			// std::cout << R_k[s*j + o] << " ";
		// }
		// std::cout << std::endl;
	// }
	
	// std::cout.precision(3);

	// printf("\n============= R (col major):\n");
	// for (size_t o = 0; o < m + s + 1; ++o) {
		// for (size_t j = 0; j < s; ++j) {
			// std::cout << R[(m + s + 1)*j + o] << " ";
		// }
		// std::cout << std::endl;
	// }

	// std::cout << "\n\n============= H update: (before)\n";
	// for (size_t i = 0; i < m+1; ++i) {
		// for (size_t j = 0; j < m; ++j) {
			// printf("%2.5f ", H[(m+1)*j + i]);
		// }
		// std::cout << std::endl;
	// }
	// std::cout << std::endl;	

// void cblas_dgemm (const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa,
									// const CBLAS_TRANSPOSE transb, const MKL_INT m, 
                  // const MKL_INT n, const MKL_INT k,
									// const double alpha, const double *a,
									// const MKL_INT lda, const double *b,
									// const MKL_INT ldb, const double beta,
									// double *c, const MKL_INT ldc);

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

	// printf("\n============= R_k (col major):\n");
	// for (size_t o = 0; o < s; ++o) {
		// for (size_t j = 0; j < s; ++j) {
			// std::cout << R_k[s*j + o] << " ";
		// }
		// std::cout << std::endl;
	// }

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

template <typename ScalarType>
sparse_status_t kernels<ScalarType>::reduce_H(ScalarType *H, size_t s, size_t m, size_t k, ScalarType *zeta, std::vector<std::pair<ScalarType, ScalarType>> &cs) {
	sparse_status_t stat = SPARSE_STATUS_SUCCESS;
	ScalarType x1;
	ScalarType x2;

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

	// std::cout << "\n\n============= zeta:\n";
	// for(size_t o = 0; o < m + 1; ++o)
		// std::cout << zeta[o] << std::endl;

	return stat;
}

template <typename ScalarType>
sparse_status_t kernels<ScalarType>::gmres_init(size_t n,
                                              const sparse_matrix_t A,
																							ScalarType *H,
																							ScalarType *H_reduced,
																							ScalarType *Q,
																							std::vector<ic_pair_t> &theta_vals,
																							size_t s,
																							size_t m) {

	sparse_status_t 			                      stat = SPARSE_STATUS_SUCCESS;
	ScalarType                                  *w;
	ScalarType                                  *H_s;          // is square with odd dimension to ensure at least one real value bc. we initially needed to pick 's' out of '2s' ritzvalues.
	ScalarType                                  h_ij;
	ScalarType                                  h_jp1j;

	/******************************/
	/*  LAPACKE_dhseqr variables  */
	/******************************/

	size_t                                      ilo = 1;        // if A hasn't been balanced by ?gebal: ilo = 1
	size_t                                      ihi = s;        // if A hasn't been balanced by ?gebal: ihi = s
	char                                        job = 'E';      // eigenvalues only are required
	const char                                  compz = 'N';    // no Schur vectors are computed
	ScalarType                                  *wr, *wi;       // real and imag part of ritz values
	// ScalarType																	*scale;
	size_t                                      i, j;           // index in for-loops

	w = (ScalarType *)mkl_calloc(n, sizeof(ScalarType), 64);if(w == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	H_s = (ScalarType *)mkl_calloc((s + 1)*s, sizeof(ScalarType), 64);if(H_s == NULL){return SPARSE_STATUS_ALLOC_FAILED;}

	/**************************/
	/*  Modfied Gram-Schmidt  */
	/**************************/

	for(j = 0; j < s; ++j) {
		mv(A, &Q[n*j], w, 1);
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
	mkl_free_buffers();
	
	return stat;
}

template <typename ScalarType>
sparse_status_t kernels<ScalarType>::mv(const sparse_matrix_t A, const ScalarType *x, ScalarType *y, size_t s) {
	sparse_status_t 			stat = SPARSE_STATUS_SUCCESS;
	struct matrix_descr 	descr;

	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	
//sparse_status_t mkl_sparse_d_mv (sparse_operation_t operation, double alpha, const sparse_matrix_t A, struct matrix_descr descr, const double *x, double beta, double *y);
	stat = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descr, x, 0, y);
	if(stat != SPARSE_STATUS_SUCCESS) throw std::invalid_argument("MatMult failed");

	mkl_free_buffers();

	return stat;
}


/****************/
/* print output */
/****************/

	// printf("\n============= H_s:\n");
	// for(size_t i = 0; i < s+1; ++i) {
		// for(size_t j = 0; j < s; ++j) {
			// printf("%2.2f ", H_s[s*i + j]);
		// }
		// printf("\n");
	// }

	// std::cout << "\n\n============= H_col_maj:\n";
	// for(i = 0; i < m+1; ++i) {
		// for(j = 0; j < m; ++j) {
			// printf("%2.2f ", H[(m+1)*j + i]);
		// }
		// std::cout << std::endl;
	// }
	// std::cout << std::endl;

		// if (n < 15) {
		// for(i = 0; i < n; ++i) {
			// for(size_t k = 0; k < m+1; ++k) {
				// std::cout << Q[k*n + i] << " ";
			// }
			// std::cout << std::endl;
		// }
		// std::cout << std::endl;
	// }

	// std::cout << std::endl;
	// for(size_t o = 0; o < n; ++o)
		// std::cout << w[o] << "\n";

	// std::cout << "\n\n============= H update: checkpoint\n";
	// for(size_t i = 0; i < m + 1; ++i) {
		// for(size_t j = 0; j < m; ++j) {
			// printf("%e ", H[(m + 1)*j + i]);
		// }
		// std::cout << std::endl;
	// }
	// std::cout << std::endl;	

	// std::cout << "\n\n============= H update: checkpoint reduced)\n";
	// for(size_t i = 0; i < m + 1; ++i) {
		// for(size_t j = 0; j < m; ++j) {
			// printf("%e ", H_reduced[(m + 1)*j + i]);
		// }
		// std::cout << std::endl;
	// }
	// std::cout << std::endl;