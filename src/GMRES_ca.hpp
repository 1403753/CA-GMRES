#ifndef GMRES_CA_HPP
#define GMRES_CA_HPP

#include "IKSPType.hpp"
#include "PCILU0_ca.hpp"
#include <papi.h>

enum Basis {MONOMIAL, NEWTON};

class GMRES_ca : public IKSPType {
	size_t										       s = 5;
	size_t										       t = 12;
	Basis                            basis = NEWTON;
	float                            SDO = 0;
	float                            SpMV = 0;
	float                            BCGS = 0;
	float                            TSQR = 0;
	float                            INIT = 0;
	float                            rtime, ptime, mflops;
	long long                        flpops;
	double                           rcond_max = 0;
	double                           rcond_min = 1;
	
public:
	GMRES_ca();
	GMRES_ca(size_t s, size_t t, Basis basis);
	virtual ~GMRES_ca();
	sparse_status_t solve(double *x, double *b);
	sparse_status_t setS(size_t s) {this->s = s; return SPARSE_STATUS_SUCCESS;};
	sparse_status_t setT(size_t t) {this->t = t; return SPARSE_STATUS_SUCCESS;};
	sparse_status_t setBasis(Basis basis) {this->basis = basis; return SPARSE_STATUS_SUCCESS;};
	float getSDO() {return this->SDO;};
	float getSpMV() {return this->SpMV;};
	float getBCGS() {return this->BCGS;};
	float getTSQR() {return this->SDO;};
	float getINIT() {return this->INIT;};
	double getRcondMin() {return this->rcond_min;};
	double getRcondMax() {return this->rcond_max;};
private:
	bool is_conj_pair(complex_t a, complex_t b);
	sparse_status_t modified_leya_ordering(size_t s, double *wr, double *wi, std::vector<ic_pair_t> &theta_vals);																	
	sparse_status_t tsqr(double *V, double *Q, double *R_, const size_t m, const size_t n, const size_t st);
	sparse_status_t gmres_init(double *H, double *H_reduced, double *Q, std::vector<ic_pair_t> &theta_vals, size_t s, size_t m);
	sparse_status_t reduce_H(double *H, size_t s, size_t m, size_t k, double *zeta, std::vector<std::pair<double, double>> &sc);
	sparse_status_t update_H(double *H, double *H_reduced, double *R, double *R_k, std::vector<ic_pair_t>  theta_vals, std::vector<double> *scales, size_t s, size_t m, size_t k);																	
};

#endif