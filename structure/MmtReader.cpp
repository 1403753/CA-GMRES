/*
 * MmtReader.cpp
 *
 *  Created on: 06.05.2018
 *      Author: Robert
 */
#include "MmtReader.hpp"

MmtReader::MmtReader() {
	
}
MmtReader::~MmtReader() {

}

sparse_status_t MmtReader::read_matrix_from_file(std::string fname, sparse_matrix_t *A_mkl) {
		
	size_t                  n, cols, nz;
	std::vector<size_t>     row_indx;
	std::vector<size_t>     col_indx;
	std::vector<double>     values;
	std::ifstream           file (fname);
	sparse_status_t         stat;
	std::string             line;
	size_t                  row_t, col_t;
	double                  val_t;
	bool                    diag = false;
	size_t                  col_ctr = 1;

	if(!file.is_open())throw std::invalid_argument("Unable to open file");
	if(file.peek() != '%')throw std::invalid_argument("Matrix Market Converter : you must verify the format of entry file.");

	getline (file, line);

	while (file.peek() == '%') getline(file, line);

	file >> n >> cols >> nz;

	if(n != cols)throw std::invalid_argument("Matrix Market Converter : input matrix is not square.");	
	if(nz < 1)throw std::invalid_argument("Matrix Market Converter : you must verify the format of entry file.");

	row_indx.reserve(nz + 10);
	col_indx.reserve(nz + 10);
	values.reserve(nz + 10);

	while (file >> row_t >> col_t >> val_t) {
		--row_t;
		--col_t;
		
		if (col_t == col_ctr) {
			if (!diag) {
				row_indx.push_back(col_t - 1);
				col_indx.push_back(col_t - 1);
				values.push_back(0);
				++nz;
			}
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
	mkl_sparse_order(*A_mkl);
		
	return stat;
}