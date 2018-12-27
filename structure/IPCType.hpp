#ifndef IPCTYPE_HPP
#define IPCTYPE_HPP

#include "KSP_.hpp"

struct Mtx_CSR;
class KSP_;

class IPCType {
public:
	KSP_ *ksp;
	virtual sparse_status_t mv(double *x, double *y, struct matrix_descr descr) = 0;
	virtual sparse_status_t setUp() = 0;
	virtual sparse_status_t precondition(double *x) = 0;
	virtual ~IPCType(){}
};

#endif