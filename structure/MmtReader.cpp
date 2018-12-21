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
	
	size_t            n, cols, nz;
	std::vector<size_t> row_indx;
	std::vector<size_t> col_indx;
	std::vector<double> values;
	// double        		*values, *row_indx, *col_indx;
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
		
	// row_indx = (size_t *) mkl_malloc(nz * sizeof(size_t), 64);if(row_indx == NULL){throw std::invalid_argument("Matrix Market Converter : malloc on 'row_indx' failed.");}
	// col_indx = (size_t *) mkl_malloc(nz * sizeof(size_t), 64);if(col_indx == NULL){throw std::invalid_argument("Matrix Market Converter : malloc on 'col_indx' failed.");}
	// values = (double *) mkl_malloc(nz * sizeof(double), 64);if(values == NULL){throw std::invalid_argument("Matrix Market Converter : malloc on 'values' failed.");}	

	row_indx.reserve(nz + 10);
	col_indx.reserve(nz + 10);
	values.reserve(nz + 10);

	size_t row_t, col_t;
	double val_t;
	bool diag = false;
	size_t col_ctr = 1;
	
	while (file >> row_t >> col_t >> val_t) {
		--row_t;
		--col_t;

		// while (col_t > col_ctr) { // as long as there are empty columns add 0 to the diagonal
			// row_indx.push_back(col_ctr - 1);
			// col_indx.push_back(col_ctr - 1);
			// values.push_back(0);
			// std::cout << "empty col!: " << col_ctr << ", "<< col_t << "\n";
			// ++col_ctr;
			// ++nz;
		// }			
			
			
		if (col_t == col_ctr) {
			if (!diag) {
				row_indx.push_back(col_t - 1);
				col_indx.push_back(col_t - 1);
				values.push_back(0);
				++nz;
			}
			// std::cout << col_t << " == " << col_ctr << std::endl;
			diag = false;
			++col_ctr;		
		}
	
		row_indx.push_back(row_t);
		col_indx.push_back(col_t);
		values.push_back(val_t);
		
		if (col_t == row_t) {
			diag = true;
		}

	}
		
	file.close();
	
	
	stat = mkl_sparse_d_create_coo(&B, SPARSE_INDEX_BASE_ZERO, n, n, nz, &row_indx[0], &col_indx[0], &values[0]);

	mkl_sparse_order(B);

	stat = mkl_sparse_convert_csr(B, SPARSE_OPERATION_NON_TRANSPOSE, A_mkl);

	mkl_sparse_destroy(B);
	// mkl_free(row_indx);
	// mkl_free(col_indx);
	// mkl_free(values);

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
	
	return stat;
}


MmtReader::~MmtReader() {
}

