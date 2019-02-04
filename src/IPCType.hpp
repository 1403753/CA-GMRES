#ifndef IPCTYPE_HPP
#define IPCTYPE_HPP

#include "KSM.hpp"

struct Mtx_CSR;
class KSM;

class IPCType {
public:
	KSM *ksp;
	virtual sparse_status_t mv(double *x, double *y, struct matrix_descr descr) = 0;
	virtual sparse_status_t setUp() = 0;
	virtual sparse_status_t precondition(double *x) = 0;
	virtual ~IPCType(){}
};

#endif