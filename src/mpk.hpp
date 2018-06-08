/*
 * mpk.hpp
 *
 *  Created on: 06.05.2018
 *      Author: Robert
 */

#ifndef MPK_HPP_
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
//#include <gsl/gsl_blas.h>
#define MPK_HPP_

class mpk {
public:
	mpk();
	static void mv(const gsl_matrix*, const gsl_vector*, gsl_vector*);
	virtual ~mpk();
};

#endif /* MPK_HPP_ */
