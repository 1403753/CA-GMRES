/*
 * gmres.hpp
 *
 *  Created on: 06.05.2018
 *      Author: Robert
 */

#ifndef GMRES_HPP_

#include "includes.hpp"
#include <vector>
#include <algorithm>

#define GMRES_HPP_

typedef std::complex<double> complex_t; 
typedef std::pair<size_t, complex_t> ic_pair_t;

template <typename ScalarType>
class gmres {
public:
	gmres();
	static bool is_conj_pair(complex_t a, complex_t b);
	static sparse_status_t mv(const sparse_matrix_t A, const ScalarType *x, ScalarType *y, size_t s);
	static sparse_status_t gmres_init(size_t n,
	                                  const sparse_matrix_t A,
																		ScalarType *H,
																		ScalarType *H_reduced,
																		ScalarType *Q,
																		std::vector<ic_pair_t, mkl_allocator<ic_pair_t>> &theta_vals,
																		size_t s,
																		size_t m);
	static sparse_status_t modified_leya_ordering(size_t s, ScalarType *wr, ScalarType *wi, std::vector<ic_pair_t, mkl_allocator<ic_pair_t>> &theta_vals);
	static sparse_status_t reduce_H(ScalarType *H, size_t s, size_t m, size_t k, ScalarType *zeta, std::vector<std::pair<ScalarType, ScalarType>, mkl_allocator<std::pair<ScalarType, ScalarType>>> &sc);
	static sparse_status_t update_H(ScalarType *H, ScalarType *H_reduced, ScalarType *R, ScalarType *R_k, std::vector<ic_pair_t, mkl_allocator<ic_pair_t>>  theta_vals, size_t s, size_t m, size_t k);
	static sparse_status_t tsqr(ScalarType *V, ScalarType *Q, ScalarType *R_, const size_t m, const size_t n, const size_t st);
	virtual ~gmres();
};

#include "gmres.tpp"
#include "gmres_modified_leja_ordering.tpp"
#include "gmres_tsqr.tpp"

#endif /* GMRES_HPP_ */
