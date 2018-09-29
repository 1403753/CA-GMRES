template <typename ScalarType>
void diagonal_pointer_csr (size_t n, size_t *rows_start, size_t *rows_end, size_t *col_indx, ScalarType *values, ScalarType **diag_ptr);

template <typename ScalarType>
sparse_status_t matrix_reader<ScalarType>::read_matrix_from_file(std::string fname, sparse_matrix_t *A, MatrixInfo<ScalarType> *minfo) {
	
	std::cout.precision(3);
	
	size_t            *row_indx, *col_indx, indx, n, cols, nnz;
	ScalarType        *values, **diag_ptr;
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
	
	row_indx = (size_t *) mkl_malloc(nnz * sizeof(size_t), 64);if(row_indx == NULL){throw std::invalid_argument("Matrix Market Converter : malloc on 'row_indx' failed.");}
	col_indx = (size_t *) mkl_malloc(nnz * sizeof(size_t), 64);if(col_indx == NULL){throw std::invalid_argument("Matrix Market Converter : malloc on 'col_indx' failed.");}
	diag_ptr = (ScalarType **) mkl_malloc(n * sizeof(ScalarType*), 64);if(col_indx == NULL){throw std::invalid_argument("Matrix Market Converter : malloc on 'diag_ptr' failed.");}
	values = (ScalarType *) mkl_malloc(nnz * sizeof(ScalarType), 64);if(values == NULL){throw std::invalid_argument("Matrix Market Converter : malloc on 'values' failed.");}	
	
	// size_t next = 0;	
	
	for(indx = 0; indx < nnz; ++indx) {
		file >> row_indx[indx] >> col_indx[indx] >> values[indx];
		--row_indx[indx];
		--col_indx[indx];

		// if (row_indx[indx] == col_indx[indx])
			// diag_ptr[next++] = &values[indx];
	}
	
	// if(next != n)throw std::invalid_argument("Matrix Market Converter : matrix has missing diagonal entries.");
	
	file.close();
	
	minfo->n = n;
	minfo->nnz = nnz;
	
	stat = mkl_sparse_d_create_coo (A, SPARSE_INDEX_BASE_ZERO, n, n, nnz, row_indx, col_indx, values);
	if (stat != SPARSE_STATUS_SUCCESS) throw std::invalid_argument("Matrix Market Converter : matrix creation failed.");
	std::cout << indx << std::endl;

	stat = mkl_sparse_convert_csr (*A, SPARSE_OPERATION_NON_TRANSPOSE, A);
	std::cout << indx << std::endl;
	
	/***************************/
	/* BUILT IN PRECONDITIONER */
	/***************************/
	// void dcsrilut (const MKL_INT *n, const double *a, const MKL_INT *ia, const MKL_INT *ja, double *bilut, MKL_INT *ibilut, MKL_INT *jbilut, const double *tol, const MKL_INT *maxfil, const MKL_INT *ipar, const double *dpar, MKL_INT *ierr);

	mkl_free(row_indx);
	mkl_free(col_indx);
	mkl_free(values);
	
	sparse_index_base_t indexing;
	
	size_t *rows_start;
	size_t *rows_end;
	
	mkl_sparse_d_export_csr (*A, &indexing, &n, &cols, &rows_start, &rows_end, &col_indx, &values);

	diagonal_pointer_csr(n, rows_start, rows_end, col_indx, values, diag_ptr);
	
	minfo->col_indx = col_indx;
	minfo->rows_start = rows_start;
	minfo->rows_end = rows_end;
	minfo->values = values;
	minfo->diag_ptr = diag_ptr;

	mkl_free_buffers();
	
	return stat;
}

template <typename ScalarType>
void diagonal_pointer_csr (size_t n, size_t *rows_start, size_t *rows_end, size_t *col_indx, ScalarType *values, ScalarType **diag_ptr) {
	size_t j1, j2;
	for (size_t i = 0; i < n; ++i) {
    j1 = rows_start[i];
    j2 = rows_end[i];

    for (size_t j = j1; j < j2; ++j) {
			if (col_indx[j] == i) {
				diag_ptr[i] = &values[j];
				break;
			}
    }
	}
}