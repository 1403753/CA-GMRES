/*
 * matrix_reader.cpp
 *
 *  Created on: 06.05.2018
 *      Author: Robert
 */
#include "matrix_reader.hpp"
matrix_reader::matrix_reader() {
	
}

sparse_status_t matrix_reader::read_matrix_from_file(std::string fname, sparse_matrix_t *A, MatrixInfo *minfo) {
	
	std::cout.precision(3);
	
	size_t						*row_indx, *col_indx, indx, rows, cols, nnz;
	double						*values;
	std::ifstream			file (fname);
	sparse_status_t 	stat;
	std::string 			line;
	sparse_matrix_t		Acoo;
	
	if (!file.is_open()) {throw std::invalid_argument("Unable to open file");}

	getline (file, line);
	std::cout << "line: " << line << '\n';
	
	while(line.front() == '%') {getline(file, line);}
	
	std::istringstream iss (line);

	iss >> rows >> cols >> nnz;
	if(nnz < 1)throw std::invalid_argument("Matrix Market Converter : you must verify the format of entry file.");
	
	row_indx = (size_t *) mkl_malloc(nnz * sizeof(size_t), 64);if(row_indx == NULL){throw std::invalid_argument("Matrix Market Converter : malloc on 'row_indx' failed.");}
	col_indx = (size_t *) mkl_malloc(nnz * sizeof(size_t), 64);if(col_indx == NULL){throw std::invalid_argument("Matrix Market Converter : malloc on 'col_indx' failed.");}
	values = (double *) mkl_malloc(nnz * sizeof(double), 64);if(values == NULL){throw std::invalid_argument("Matrix Market Converter : malloc on 'values' failed.");}	

	indx = 0;
	while (getline(file, line)) {
		iss.clear();
		iss.str(line);
		iss >> row_indx[indx] >> col_indx[indx] >> values[indx];		
		--row_indx[indx];
		--col_indx[indx];
		++indx;
	}
	
	file.close();
	
	minfo->rows = rows;
	minfo->cols = cols;
	minfo->nnz = nnz;
	
	std::cout << "rows: " << rows << ", cols: " << cols << ", nnz: " << nnz << "\n";

	// for(size_t i = 0; i < nnz; ++i)
		// std::cout << "row: " << row_indx[i] << ", col: " << col_indx[i] << ", val: " << values[i] << "\n";	
	
	stat = mkl_sparse_d_create_coo (&Acoo, SPARSE_INDEX_BASE_ZERO, rows, cols, nnz, row_indx, col_indx, values);
	if (stat != SPARSE_STATUS_SUCCESS) throw std::invalid_argument("Matrix Market Converter : matrix creation failed.");

	stat = mkl_sparse_convert_csr (Acoo, SPARSE_OPERATION_NON_TRANSPOSE, A);
	
	mkl_free(row_indx);
	mkl_free(col_indx);
	mkl_free(values);
	mkl_free_buffers();
	
	return stat;
}


matrix_reader::~matrix_reader() {
	
}

