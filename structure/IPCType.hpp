#ifndef IPCTYPE_HPP
#define IPCTYPE_HPP

#include "KSP.hpp"

struct Mtx_CSR;
class KSP;

class IPCType {
public:
	KSP *ksp;
	virtual sparse_status_t mpk(double *x, double *y, size_t s) = 0;
	virtual sparse_status_t setUp() = 0;
	virtual sparse_status_t precondition(double *x) = 0; // remove this
};


#endif