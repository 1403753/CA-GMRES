template <typename ScalarType>
sparse_status_t gmres<ScalarType>::init_gmres(size_t n, const sparse_matrix_t A, const ScalarType *v, ScalarType **H_s, ScalarType **Q, size_t s) {
	sparse_status_t 			stat = SPARSE_STATUS_SUCCESS;
	ScalarType 						beta;
	ScalarType            *w;
	// ScalarType            theta_vals;
	size_t                i, j;
	
	w = (ScalarType *)mkl_calloc(n, sizeof(ScalarType), 64);if(w == NULL){return SPARSE_STATUS_ALLOC_FAILED;}

	
	beta = cblas_dnrm2 (n, v, 1);
	
	for(i = 0; i < n; ++i)
		(*Q)[i] = v[i] / beta;
		
	for(j = 0; j < s; ++j) {
		mv(A, Q[0], &w, 1);
	}
	
	mkl_free(w);
	mkl_free_buffers();
	
	return stat;
}




template <typename ScalarType>
sparse_status_t gmres<ScalarType>::mv(const sparse_matrix_t A, const ScalarType *x, ScalarType **y, size_t s) {
	sparse_status_t 			stat;
	struct matrix_descr 	descr;
	
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	
	mkl_sparse_set_mv_hint(A, SPARSE_OPERATION_NON_TRANSPOSE, descr, s);
	mkl_sparse_optimize(A);
		
//sparse_status_t mkl_sparse_d_mv (sparse_operation_t operation, double alpha, const sparse_matrix_t A, struct matrix_descr descr, const double *x, double beta, double *y);
	stat = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descr, x, 0, *y);
	if(stat != SPARSE_STATUS_SUCCESS) throw std::invalid_argument("MatMult failed");
	
	mkl_free_buffers();
	return stat;
}
