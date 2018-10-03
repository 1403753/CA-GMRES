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
typedef std::pair<size_t, complex_t> pair_t;

template <typename ScalarType>
class gmres {
public:
	gmres();
	static bool is_conj_pair(complex_t a, complex_t b);
	static sparse_status_t mv(const sparse_matrix_t A, const ScalarType *x, ScalarType *y, size_t s);
	static sparse_status_t gmres_init(size_t n,
	                                  const sparse_matrix_t A,
																		ScalarType *H,
																		ScalarType *Q,
																		std::vector<pair_t, mkl_allocator<pair_t>> &theta_vals,
																		size_t s,
																		size_t m);

	static sparse_status_t modified_leya_ordering(size_t s, ScalarType *wr, ScalarType *wi, std::vector<pair_t, mkl_allocator<pair_t>> &theta_vals);
	// TODO -> static reduce_H(size_t ritz_num, ScalarType *H, size_t m, ScalarType *zeta);
	virtual ~gmres();
};

#include "gmres.tpp"
#include "gmres_modified_leja_ordering.tpp"

#endif /* GMRES_HPP_ */
