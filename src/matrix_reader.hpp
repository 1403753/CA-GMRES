/*
 * matrix_reader.hpp
 *
 *  Created on: 06.05.2018
 *      Author: Robert
 */

#ifndef MATRIX_READER_HPP_

#include "includes.hpp"
#include <fstream>
#include <string>

#define MATRIX_READER_HPP_

template <typename ScalarType>
struct MatrixInfo{
	size_t n;
	size_t nnz;
	size_t *rows_start;
	size_t *rows_end;
	size_t *col_indx;
	ScalarType *values;
	ScalarType **diag_ptr;
};

template <typename ScalarType>
class matrix_reader {
public:
	matrix_reader();
	MatrixInfo<ScalarType> minfo;
	static sparse_status_t read_matrix_from_file(std::string fname, sparse_matrix_t *A, MatrixInfo<ScalarType> *minfo);
	virtual ~matrix_reader();
};

#include "matrix_reader.tpp"

#endif /* MATRIX_READER_HPP_ */
