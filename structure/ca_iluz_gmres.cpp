/*
 * ca_iluz_gmres.cpp
 *
 *  Created on: 17.10.2018
 *      Author: Robert
 */

 
#include "ca_iluz_gmres.hpp"
#include "solve.cc"
#include "tsqr.cc"

ca_iluz_gmres::ca_iluz_gmres(Mtx_CSR *A_mtx) : A_mtx{A_mtx} { 

}

ca_iluz_gmres::~ca_iluz_gmres() {

}

sparse_status_t ca_iluz_gmres::setUp() {
	
	return SPARSE_STATUS_SUCCESS;
}

sparse_status_t ca_iluz_gmres::solve(double *x, double *b) {
	
	return SPARSE_STATUS_SUCCESS;
}

