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
public:
	MmtReader();
	static sparse_status_t read_matrix_from_file(std::string fname, Mtx_CSR *A_mtx);
	virtual ~MmtReader();
};

#endif // MMT_READER_HPP