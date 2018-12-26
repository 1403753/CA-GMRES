#ifndef GMRES_HPP
#define GMRES_HPP

#include "IKSPType.hpp"
#include "PCILU0_ca.hpp"

class GMRES : public IKSPType {
	size_t										       m = 60;
public:
	GMRES();
	GMRES(size_t m);
	virtual ~GMRES();
	sparse_status_t solve(double *x, double *b);
// private:
	// sparse_status_t gmres_init(double *H, double *H_reduced, double *Q, std::vector<ic_pair_t> &theta_vals, size_t s, size_t m);
	// sparse_status_t reduce_H(double *H, size_t s, size_t m, size_t k, double *zeta, std::vector<std::pair<double, double>> &sc);
};

#endif