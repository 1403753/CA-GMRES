/*
 * ca_iluz_gmres.hpp
 *
 *  Created on: 17.10.2018
 *      Author: Robert
 */

#ifndef CA_ILUZ_GMRES_HPP_

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
// #include <gsl/gsl_randist.h> 	// <- not needed everywhere!
// #include <gsl/gsl_rng.h> 			// <- not needed everywhere!
#include <complex>
#include <exception>

#define MKL_MAX_PATH_LEN 4096
#define MKL_Complex16 std::complex<double>
#define MKL_INT size_t

#include "mkl.h"

#define CA_ILUZ_GMRES_HPP_

struct Mtx_CSR{
	size_t n;
	size_t nz;
	size_t *row_ptr;
	size_t *col_indx;
	double *values;
	Mtx_CSR(size_t n, size_t nz, size_t *row_ptr, size_t *col_indx, double *values)
		: n{n}, nz{nz}, row_ptr{row_ptr}, col_indx{col_indx}, values{values}
	{ }
};

typedef std::complex<double> complex_t; 
typedef std::pair<size_t, complex_t> ic_pair_t;

enum Gr_Part {GRAPH_LOWER, GRAPH_UPPER, GRAPH_COMPLETE};


class ca_iluz_gmres {
	size_t                    s = 5;
	size_t                    t = 12;
	double                    rTol = 1e-10;
	double                    aTol = 1e-50;
	double                    dTol = 1e+4;
	size_t                    maxit = 10000;
	const Mtx_CSR             *A_mtx;
	double                    *iluz_vals;
public:
	ca_iluz_gmres(Mtx_CSR *A_mtx);
	sparse_status_t setUp();
	sparse_status_t solve(double *x, double *b);
	virtual ~ca_iluz_gmres();
private:
	bool is_conj_pair(complex_t a, complex_t b);
	// sparse_status_t gmres_init(size_t n,
	                                  // const sparse_matrix_t A,
																		// double *H,
																		// double *H_reduced,
																		// double *Q,
																		// std::vector<ic_pair_t> &theta_vals,
																		// size_t s,
																		// size_t m);
	// sparse_status_t reduce_H(double *H,
	                                // size_t s,
																	// size_t m,
																	// size_t k,
																	// double *zeta,
																	// std::vector<std::pair<double, double>> &sc);
	// sparse_status_t update_H(double *H,
																	// double *H_reduced,
																	// double *R,
																	// double *R_k,
																	// std::vector<ic_pair_t>  theta_vals,
																	// size_t s,
																	// size_t m,
																	// size_t k);
	// sparse_status_t mv(const sparse_matrix_t A, const double *x, double *y, size_t s);																	
	sparse_status_t modified_leya_ordering(size_t s, double *wr, double *wi, std::vector<ic_pair_t> &theta_vals);																	
	sparse_status_t tsqr(double *V, double *Q, double *R_, const size_t m, const size_t n, const size_t st);
	
	};

#endif /* CA_ILUZ_GMRES_HPP_ */