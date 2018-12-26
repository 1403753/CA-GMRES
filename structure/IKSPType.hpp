#ifndef IKSPTYPE_HPP
#define IKSPTYPE_HPP

#include "KSP_ca.hpp"

class KSP_ca;

struct Mtx_CSR;

class IKSPType {
public:
	KSP_ca *ksp;
	virtual sparse_status_t solve(double *x, double *b) = 0;
};


#endif