#ifndef IKSPTYPE_HPP
#define IKSPTYPE_HPP

#include "KSM.hpp"

class KSM;

struct Mtx_CSR;

class IKSPType {
public:
	KSM *ksp;
	virtual sparse_status_t solve(double *x, double *b) = 0;
};


#endif