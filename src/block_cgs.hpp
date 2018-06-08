/*
 * block_cgs.hpp
 *
 *  Created on: 07.05.2018
 *      Author: Robert
 */

#ifndef BLOCK_CGS_HPP_
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#define BLOCK_CGS_HPP_

class block_cgs {
public:
	block_cgs();
	virtual ~block_cgs();
};

#endif /* BLOCK_CGS_HPP_ */
