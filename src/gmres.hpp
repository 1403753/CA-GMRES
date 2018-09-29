/*
 * gmres.hpp
 *
 *  Created on: 06.05.2018
 *      Author: Robert
 */

#ifndef GMRES_HPP_

#include "includes.hpp"

#define GMRES_HPP_

template <typename ScalarType>
class gmres {
public:
	gmres();
	static sparse_status_t mv(const sparse_matrix_t A, const ScalarType *x, ScalarType **y, size_t s);
	static sparse_status_t init_gmres(size_t n, const sparse_matrix_t A, const ScalarType *v, ScalarType **H_s, ScalarType **Q, size_t s);
	virtual ~gmres();
};

#include "gmres.tpp"

#endif /* GMRES_HPP_ */
