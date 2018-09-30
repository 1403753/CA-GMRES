template <typename ScalarType>
bool gmres<ScalarType>::is_conj_pair(complex_t a, complex_t b) {
	return (a.real() == b.real() && a.imag() == -b.imag() && a.imag() != 0);
}

template <typename ScalarType>
sparse_status_t gmres<ScalarType>::init_gmres(size_t n, const sparse_matrix_t A, const ScalarType *v, ScalarType *H, ScalarType *Q, std::vector<pair_t, mkl_allocator<pair_t>> &theta_vals, size_t s, size_t m) {
	sparse_status_t 			                      stat = SPARSE_STATUS_SUCCESS;
	ScalarType 						                      beta;
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
	H_s = (ScalarType *)mkl_calloc((s + 1)*(s + 1), sizeof(ScalarType), 64);if(H_s == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	
	beta = cblas_dnrm2 (n, v, 1);
	
	// Q is in col-major!
	for(i = 0; i < n; ++i)
		Q[i] = v[i] / beta;
	
	if (n < 15) {
		for(i = 0; i < n; ++i) {
			for(size_t k = 0; k < m+1; ++k) {
				std::cout << Q[k*n + i] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	
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
			H_s[(s+1)*i + j] = h_ij;
			H[m*i + j] = h_ij;
			
			cblas_daxpy (n, h_ij*(-1), &Q[i*n], 1, w, 1);
		}
		h_jp1j = cblas_dnrm2(n, w, 1);
		H_s[(s+1)*(j+1) + j] = h_jp1j;
		H[m*(j+1) + j] = h_jp1j;
		
		cblas_daxpy(n, 1 / h_jp1j, w, 1, &Q[(j+1)*n], 1);	
	}
	
  /********************************************************************/
	/* compute one more column of H_s to get odd number of Ritz values  */
	/********************************************************************/
	
	mv(A, &Q[j*n], w, 1);
	for(i = 0; i < j+1; ++i) {
		h_ij = cblas_ddot(n, w, 1, &Q[i*n], 1);
		H_s[(s+1)*i + j] = h_ij;			
		cblas_daxpy (n, h_ij*(-1), &Q[i*n], 1, w, 1);
	}
	h_jp1j = cblas_dnrm2(n, w, 1);

	printf("\n============= H_s:\n");
	for(size_t i = 0; i < s+2; ++i) {
		for(size_t j = 0; j < s+1; ++j) {
			printf("%2.2f ", H_s[i*(s+1) + j]);
		}
		printf("\n");
	}
	
	
	wr = (double *)mkl_malloc((s + 1)*sizeof(double), 64);if(wr == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	wi = (double *)mkl_malloc((s + 1)*sizeof(double), 64);if(wi == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	// scale = (double *)mkl_malloc((s + 1)*sizeof(double), 64);if(scale == NULL){return SPARSE_STATUS_ALLOC_FAILED;}



// balancing prob. not needed!
	// job = 'N'; // neither scaled nor pivoted
//LAPACKE_dgebal( int matrix_layout, char job, lapack_int n, double* a, lapack_int lda, lapack_int* ilo, lapack_int* ihi, double* scale );
	// LAPACKE_dgebal(LAPACK_ROW_MAJOR, job, s+1, H_s, s+1, &ilo, &ihi, scale);	
	
	job = 'E';
	
//LAPACKE_dhseqr(int matrix_layout, char job, char compz, lapack_int n, lapack_int ilo, lapack_int ihi, double *h, lapack_int ldh, double *wr, double *wi, double *z, lapack_int ldz);
	LAPACKE_dhseqr(LAPACK_ROW_MAJOR, job, compz, s+1, ilo, ihi, H_s, s+1, wr, wi, nullptr, s+1);	
	
	std::cout << std::endl;
	for(size_t o = 0; o < s; ++o)
		printf("%.10f,\t%.10f\n", wr[o], wi[o]);	

	modified_leya_ordering(s/2+1, wr, wi, theta_vals);
	
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
