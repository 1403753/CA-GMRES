/*
 * MmtReader.hpp
 *
 *  Created on: 06.05.2018
 *      Author: Robert
 */

#ifndef MMT_READER_HPP

#include "KSP.hpp"
#include <fstream>
#include <string>

#define MMT_READER_HPP

class MmtReader {
	sparse_matrix_t B;
	sparse_matrix_t A;
public:
	MmtReader();
	virtual ~MmtReader();
	sparse_status_t read_matrix_from_file(std::string fname, sparse_matrix_t *A_mkl);
};

#endif // MMT_READER_HPP