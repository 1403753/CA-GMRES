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

typedef struct _MatrixInfo{
	size_t rows;
	size_t cols;
	size_t nnz;
} MatrixInfo;

template <typename ScalarType>
class matrix_reader {
public:
	matrix_reader();
	static sparse_status_t read_matrix_from_file(std::string fname, sparse_matrix_t *A, MatrixInfo *minfo);
	virtual ~matrix_reader();
};

#include "matrix_reader.tpp"

#endif /* MATRIX_READER_HPP_ */
