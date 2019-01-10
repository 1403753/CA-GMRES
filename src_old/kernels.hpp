/*
 * kernels.hpp
 *
 *  Created on: 06.05.2018
 *      Author: Robert
 */

#ifndef KERNELS_HPP_

#include "includes.hpp"
#include <vector>
#include <algorithm>

#define KERNELS_HPP_

template <typename ScalarType>
class kernels {
public:
	kernels();
	static bool is_conj_pair(complex_t a, complex_t b);
	static sparse_status_t gmres_init(size_t n,
	                                  const sparse_matrix_t A,
																		ScalarType *H,
																		ScalarType *H_reduced,
																		ScalarType *Q,
																		std::vector<ic_pair_t> &theta_vals,
																		size_t s,
																		size_t m);
	static sparse_status_t reduce_H(ScalarType *H,
	                                size_t s,
																	size_t m,
																	size_t k,
																	ScalarType *zeta,
																	std::vector<std::pair<ScalarType, ScalarType>> &sc);
	static sparse_status_t update_H(ScalarType *H,
																	ScalarType *H_reduced,
																	ScalarType *R,
																	ScalarType *R_k,
																	std::vector<ic_pair_t>  theta_vals,
																	size_t s,
																	size_t m,
																	size_t k);
	static sparse_status_t mv(const sparse_matrix_t A, const ScalarType *x, ScalarType *y, size_t s);																	
	static sparse_status_t modified_leya_ordering(size_t s, ScalarType *wr, ScalarType *wi, std::vector<ic_pair_t> &theta_vals);																	
	static sparse_status_t tsqr(ScalarType *V, ScalarType *Q, ScalarType *R_, const size_t m, const size_t n, const size_t st);
	virtual ~kernels();
};

#include "kernels.tpp"
#include "kernels_modified_leja_ordering.tpp"
#include "kernels_tsqr.tpp"

#endif /* KERNELS_HPP_ */
