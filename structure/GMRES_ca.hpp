#ifndef GMRES_CA_HPP
#define GMRES_CA_HPP

#include "IKSPType.hpp"

class GMRES_ca : public IKSPType {
	const Mtx_CSR *A_mtx;
public:
	GMRES_ca();
	// void setA_mtx(Mtx_CSR *A_mtx);
	sparse_status_t solve(double *x, double *b);
	virtual ~GMRES_ca();
// private:
	// bool is_conj_pair(complex_t a, complex_t b);
	// sparse_status_t modified_leya_ordering(size_t s, double *wr, double *wi, std::vector<ic_pair_t> &theta_vals);																	
	// sparse_status_t tsqr(double *V, double *Q, double *R_, const size_t m, const size_t n, const size_t st);

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
	
};

#endif