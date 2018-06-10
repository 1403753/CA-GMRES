/*
 * tsqr.hpp
 *
 *  Created on: 06.05.2018
 *      Author: Robert
 */

#ifndef TSQR_HPP_
#include "mkl.h"
#include "mkl_spblas.h"
#define TSQR_HPP_

class tsqr {
public:
	tsqr();
	static void qr(double **A, size_t M, size_t N);
	virtual ~tsqr();
};

#endif /* TSQR_HPP_ */
