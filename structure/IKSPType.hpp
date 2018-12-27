#ifndef IKSPTYPE_HPP
#define IKSPTYPE_HPP

#include "KSP_.hpp"

class KSP_;

struct Mtx_CSR;

class IKSPType {
public:
	KSP_ *ksp;
	virtual sparse_status_t solve(double *x, double *b) = 0;
};


#endif