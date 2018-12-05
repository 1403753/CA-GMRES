/*
 * GMRES_ca.cpp
 *
 *  Created on: 17.10.2018
 *      Author: Robert
 */

 
#include "GMRES_ca.hpp"
#include "solve.cc"
#include "tsqr.cc"

GMRES_ca::GMRES_ca() { 

}

// void GMRES_ca::setA_mtx(Mtx_CSR *A_mtx) {
	// this->A_mtx = A_mtx;
// };

GMRES_ca::~GMRES_ca() {

}

sparse_status_t GMRES_ca::solve(double *x, double *b) {
	
	return SPARSE_STATUS_SUCCESS;
}

