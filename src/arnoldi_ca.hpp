/*
 * arnoldi_ca.hpp
 *
 *  Created on: 11.05.2018
 *      Author: Robert
 */

#ifndef ARNOLDI_CA_HPP_
#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <papi.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define MKL_INT size_t
#define MKL_Complex16 std::complex<double>
#include "mkl.h"
#define ARNOLDI_CA_HPP_

class arnoldi_ca {
public:
	arnoldi_ca();
	static void givens_rotations();	
	static void modified_leja_ordering();	
	virtual ~arnoldi_ca();
};

#endif /* ARNOLDI_CA_HPP_ */