#ifndef IKSPTYPE_HPP
#define IKSPTYPE_HPP

#include "KSP.hpp"

class KSP;

struct Mtx_CSR;

class IKSPType {
public:
	KSP *ksp;
	// virtual void setA_mtx(Mtx_CSR *A_mtx);
	virtual sparse_status_t solve(double *x, double *b) = 0;
};


#endif