/*
 * spmv.hpp
 *
 *  Created on: 06.05.2018
 *      Author: Robert
 */

#ifndef SPMV_HPP_

#include "includes.hpp"

#define SPMV_HPP_

class spmv {
public:
	spmv();
	// static void mv(const std::vector<double>, const std::vector<double>, std::vector<double>&);
	static sparse_status_t mv(const sparse_matrix_t A, const double *x, double **y);
	virtual ~spmv();
};

#endif /* SPMV_HPP_ */
