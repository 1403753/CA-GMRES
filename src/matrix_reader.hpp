/*
 * matrix_reader.hpp
 *
 *  Created on: 06.05.2018
 *      Author: Robert
 */

#ifndef MATRIX_READER_HPP_

#include "includes.hpp"

#define MATRIX_READER_HPP_

class matrix_reader {
public:
	matrix_reader();
	static void read_matrix_from_file();
	virtual ~matrix_reader();
};

#endif /* MATRIX_READER_HPP_ */
