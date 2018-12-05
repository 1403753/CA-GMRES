/*
 * MmtReader.cpp
 *
 *  Created on: 06.05.2018
 *      Author: Robert
 */
#include "MmtReader.hpp"

MmtReader::MmtReader() {
	
}

sparse_status_t MmtReader::read_matrix_from_file(std::string fname, sparse_matrix_t *A_mkl) {
	
	std::cout.precision(3);
	
	size_t            *row_indx, *col_indx, n, cols, nz;
	double        		*values;
	std::ifstream     file (fname);
	sparse_status_t   stat;
	std::string       line;
	
	if(!file.is_open())throw std::invalid_argument("Unable to open file");
	if(file.peek() != '%')throw std::invalid_argument("Matrix Market Converter : you must verify the format of entry file.");
	
	getline (file, line);
	
	while (file.peek() == '%') getline(file, line);

	file >> n >> cols >> nz;

	if(n != cols)throw std::invalid_argument("Matrix Market Converter : input matrix is not square.");	
	if(nz < 1)throw std::invalid_argument("Matrix Market Converter : you must verify the format of entry file.");
		
	row_indx = (size_t *) mkl_malloc(nz * sizeof(size_t), 64);if(row_indx == NULL){throw std::invalid_argument("Matrix Market Converter : malloc on 'row_indx' failed.");}
	col_indx = (size_t *) mkl_malloc(nz * sizeof(size_t), 64);if(col_indx == NULL){throw std::invalid_argument("Matrix Market Converter : malloc on 'col_indx' failed.");}
	values = (double *) mkl_malloc(nz * sizeof(double), 64);if(values == NULL){throw std::invalid_argument("Matrix Market Converter : malloc on 'values' failed.");}	
	
	for(size_t i = 0; i < nz; ++i) {
		file >> row_indx[i] >> col_indx[i] >> values[i];
		--row_indx[i];
		--col_indx[i];
	}
	
	file.close();
	
	stat = mkl_sparse_d_create_coo(&B, SPARSE_INDEX_BASE_ZERO, n, n, nz, row_indx, col_indx, values);
	
	stat = mkl_sparse_convert_csr(B, SPARSE_OPERATION_NON_TRANSPOSE, A_mkl);

	mkl_sparse_destroy(B);
	mkl_free(row_indx);
	mkl_free(col_indx);
	mkl_free(values);

	mkl_sparse_order(*A_mkl);
	
	// sparse_index_base_t indexing;

	// size_t *rows_start;
	// size_t *rows_end;

	
	// mkl_sparse_d_export_csr(A_mkl, &indexing, &n, &cols, &rows_start, &rows_end, &col_indx, &values);

	// A_mtx->n = n;
	// A_mtx->nz = nz;	
	
	// A_mtx->row_ptr = (size_t *) mkl_malloc((n + 1) * sizeof(size_t), 64);if(A_mtx->row_ptr == NULL){throw std::invalid_argument("Matrix Market Converter : malloc on 'row_indx' failed.");}
	// A_mtx->col_indx = (size_t *) mkl_malloc(nz * sizeof(size_t), 64);if(A_mtx->col_indx == NULL){throw std::invalid_argument("Matrix Market Converter : malloc on 'col_indx' failed.");}
	// A_mtx->values = (double *) mkl_malloc(nz * sizeof(double), 64);if(A_mtx->values == NULL){throw std::invalid_argument("Matrix Market Converter : malloc on 'values' failed.");}	
	
	// for (size_t i = 0; i < n + 1; ++i)
		// A_mtx->row_ptr[i] = rows_start[i];
	
	// for (size_t i = 0; i < nz; ++i) {
		// A_mtx->col_indx[i] = col_indx[i];
		// A_mtx->values[i] = values[i];
	// }

	// A_mtx->values = values;
	// A_mtx->row_ptr = rows_start;
	// A_mtx->col_indx = col_indx;
	
	mkl_free_buffers();

	
	return stat;
}


MmtReader::~MmtReader() {
	mkl_free_buffers();
}

