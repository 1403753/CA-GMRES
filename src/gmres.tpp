template <typename ScalarType>
sparse_status_t gmres<ScalarType>::reduce_H(ScalarType *H, size_t s, size_t m, size_t k, ScalarType *zeta) {
	sparse_status_t stat = SPARSE_STATUS_EXECUTION_FAILED;
	ScalarType c_giv = 0;
	ScalarType s_giv = 0;
	ScalarType x1 = 0;
	ScalarType x2 = 0;

	for (size_t i = 0; i < s; ++i) {
		
		x1 = H[(i + s*k)*m + (i + s*k)];
		x2 = H[(i + s*k + 1)*m + (i + s*k)];

// void cblas_drotg (double *x1, double *x2, double *c, double *s);
		cblas_drotg(&x1, &x2, &c_giv, &s_giv);

// void cblas_drot (const size_t n, double *x, const size_t incx, double *y, const size_t incy, const double c, const double s);
		cblas_drot(m - (i + s*k), &H[(i + s*k)*m + (i + s*k)], 1, &H[(i + s*k + 1)*m + (i + s*k)], 1, c_giv, s_giv);

// void cblas_drot (const size_t n, double *x, const size_t incx, double *y, const size_t incy, const double c, const double s);
		cblas_drot(1, &zeta[i + s*k], 1, &zeta[i + s*k + 1], 1, c_giv, s_giv);

	}

	stat = SPARSE_STATUS_SUCCESS;
	return stat;
}

template <typename ScalarType>
sparse_status_t gmres<ScalarType>::gmres_init(size_t n,
                                              const sparse_matrix_t A,
																							ScalarType *H,
																							ScalarType *Q,
																							std::vector<pair_t, mkl_allocator<pair_t>> &theta_vals,
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
	
	/**************************/
	/*  Modfied Gram-Schmidt  */
	/**************************/
	
	for(j = 0; j < s; ++j) {
		mv(A, &Q[j*n], w, 1);
		for(i = 0; i < j+1; ++i) {
			
			h_ij = cblas_ddot(n, w, 1, &Q[i*n], 1);
			H_s[s*i + j] = h_ij;
			H[m*i + j] = h_ij;
			
			cblas_daxpy (n, h_ij*(-1), &Q[i*n], 1, w, 1);
		}
		
		h_jp1j = cblas_dnrm2(n, w, 1);
		H_s[s*(j+1) + j] = h_jp1j;
		H[m*(j+1) + j] = h_jp1j;
		
// void	cblas_dscal(const int N, const double alpha, double *X, const int incX)
		// cblas_dscal(n, 1 / h_jp1j, w, 1);		
		// cblas_dcopy(n, w, 1, &Q[(j+1)*n], 1);
		cblas_daxpy(n, 1 / h_jp1j, w, 1, &Q[(j+1)*n], 1);	
	}

	printf("\n============= H_s:\n");
	for(size_t i = 0; i < s+1; ++i) {
		for(size_t j = 0; j < s; ++j) {
			printf("%2.2f ", H_s[s*i + j]);
		}
		printf("\n");
	}
	
	
	wr = (double *)mkl_malloc(s*sizeof(double), 64);if(wr == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	wi = (double *)mkl_malloc(s*sizeof(double), 64);if(wi == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	// scale = (double *)mkl_malloc((s + 1)*sizeof(double), 64);if(scale == NULL){return SPARSE_STATUS_ALLOC_FAILED;}



// balancing prob. not needed!
	// job = 'N'; // neither scaled nor pivoted
//LAPACKE_dgebal( int matrix_layout, char job, lapack_int n, double* a, lapack_int lda, lapack_int* ilo, lapack_int* ihi, double* scale );
	// LAPACKE_dgebal(LAPACK_ROW_MAJOR, job, s, H_s, s, &ilo, &ihi, scale);	
	
	job = 'E';
	
//LAPACKE_dhseqr(int matrix_layout, char job, char compz, lapack_int n, lapack_int ilo, lapack_int ihi, double *h, lapack_int ldh, double *wr, double *wi, double *z, lapack_int ldz);
	LAPACKE_dhseqr(LAPACK_ROW_MAJOR, job, compz, s, ilo, ihi, H_s, s, wr, wi, nullptr, s);	
	
	std::cout << std::endl;
	for(size_t o = 0; o < s; ++o)
		printf("%.10f,\t%.10f\n", wr[o], wi[o]);	

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
sparse_status_t gmres<ScalarType>::mv(const sparse_matrix_t A, const ScalarType *x, ScalarType *y, size_t s) {
	sparse_status_t 			stat = SPARSE_STATUS_SUCCESS;
	struct matrix_descr 	descr;
	
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	
	mkl_sparse_set_mv_hint(A, SPARSE_OPERATION_NON_TRANSPOSE, descr, s);
	mkl_sparse_optimize(A);
		
//sparse_status_t mkl_sparse_d_mv (sparse_operation_t operation, double alpha, const sparse_matrix_t A, struct matrix_descr descr, const double *x, double beta, double *y);
	stat = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descr, x, 0, y);
	if(stat != SPARSE_STATUS_SUCCESS) throw std::invalid_argument("MatMult failed");
	
	mkl_free_buffers();
	return stat;
}
