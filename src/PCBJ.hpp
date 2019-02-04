#ifndef PCBJ_HPP
#define PCBJ_HPP

#include "IPCType.hpp"

class PCBJ : public IPCType{
	bool setup = false;
	sparse_matrix_t *M_mkl;
public:
	PCBJ();
	virtual ~PCBJ();
	sparse_status_t setUp();
	sparse_status_t precondition(double *x);
	sparse_status_t mv(double *x, double *y, struct matrix_descr descr);
};

#endif