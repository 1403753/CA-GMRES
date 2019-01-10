template <typename ScalarType>
sparse_status_t matrix_reader<ScalarType>::read_matrix_from_file(std::string fname, sparse_matrix_t *A, MatrixInfo<ScalarType> *minfo) {
	
	std::cout.precision(3);
	
	size_t            *row_indx, *col_indx, i, n, cols, nnz;
	ScalarType        *values;
	std::ifstream     file (fname);
	sparse_status_t   stat;
	std::string       line;
	
	if(!file.is_open())throw std::invalid_argument("Unable to open file");
	if(file.peek() != '%')throw std::invalid_argument("Matrix Market Converter : you must verify the format of entry file.");
	
	getline (file, line);
	std::cout << "line: " << line << std::endl;
	
	while (file.peek() == '%') getline(file, line);

	file >> n >> cols >> nnz;

	if(n != cols)throw std::invalid_argument("Matrix Market Converter : input matrix is not square.");	
	if(nnz < 1)throw std::invalid_argument("Matrix Market Converter : you must verify the format of entry file.");
	
	std::cout << "n: " << n << std::endl;
	
	row_indx = (size_t *) mkl_malloc(nnz * sizeof(size_t), 64);if(row_indx == NULL){throw std::invalid_argument("Matrix Market Converter : malloc on 'row_indx' failed.");}
	col_indx = (size_t *) mkl_malloc(nnz * sizeof(size_t), 64);if(col_indx == NULL){throw std::invalid_argument("Matrix Market Converter : malloc on 'col_indx' failed.");}
	values = (ScalarType *) mkl_malloc(nnz * sizeof(ScalarType), 64);if(values == NULL){throw std::invalid_argument("Matrix Market Converter : malloc on 'values' failed.");}	
	
	for(i = 0; i < nnz; ++i) {
		file >> row_indx[i] >> col_indx[i] >> values[i];
		--row_indx[i];
		--col_indx[i];
	}
	
	file.close();
	
	minfo->n = n;
	minfo->nnz = nnz;
	
	sparse_matrix_t B;
	
	stat = mkl_sparse_d_create_coo(&B, SPARSE_INDEX_BASE_ZERO, n, n, nnz, row_indx, col_indx, values);
	if (stat != SPARSE_STATUS_SUCCESS) throw std::invalid_argument("Matrix Market Converter : matrix creation failed.");
	
	stat = mkl_sparse_convert_csr (B, SPARSE_OPERATION_NON_TRANSPOSE, A);
	if (stat != SPARSE_STATUS_SUCCESS) throw std::invalid_argument("Matrix Market Converter : matrix conversion to csr failed.");

	/***************************/
	/* BUILT IN PRECONDITIONER */
	/***************************/
	// void dcsrilut (const MKL_INT *n, const double *a, const MKL_INT *ia, const MKL_INT *ja, double *bilut, MKL_INT *ibilut, MKL_INT *jbilut, const double *tol, const MKL_INT *maxfil, const MKL_INT *ipar, const double *dpar, MKL_INT *ierr);
	
	sparse_index_base_t indexing;
	
	size_t *rows_start;
	size_t *rows_end;

	mkl_sparse_order(*A);	

	// mkl_sparse_d_export_csr (*A, &indexing, &n, &cols, &rows_start, &rows_end, &col_indx, &values);

	
	// minfo->col_indx = col_indx;
	// minfo->rows_start = rows_start;
	// minfo->rows_end = rows_end;
	// minfo->values = values;

	mkl_free(row_indx);
	mkl_free(col_indx);
	mkl_free(values);
	mkl_sparse_destroy(B);
	mkl_free_buffers();
	
	return stat;
}