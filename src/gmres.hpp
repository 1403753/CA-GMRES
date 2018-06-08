/*
 * gmres.hpp
 *
 *  Created on: 06.05.2018
 *      Author: Robert
 */

#ifndef GMRES_HPP_
#include <iostream>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include <iomanip>
#include "mpk.hpp"
#include "arnoldi_ca.hpp"
#include "spmv.hpp"
#define GMRES_HPP_

class gmres {
public:
	gmres();
	virtual ~gmres();
};

#endif /* GMRES_HPP_ */
