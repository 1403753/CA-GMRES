/*
 * mpk.cpp
 *
 *  Created on: 06.05.2018
 *      Author: Robert
 */
#include "mpk.hpp"
mpk::mpk() {
	
}

void mpk::mv(const gsl_matrix *A, const gsl_vector *u, gsl_vector *y) {
	/*
	* TODO:
	*/
	//mkl_d_gemv(CblasNoTrans, 1, A, u, 0, y);
}


mpk::~mpk() {
	
}

