#ifndef PCILU0_CA_HPP
#define PCILU0_CA_HPP

#include "IPCType.hpp"

class PCILU0_ca : public IPCType{
	bool setup = false;
	sparse_matrix_t *M_mkl;
public:
	PCILU0_ca();
	virtual ~PCILU0_ca();
	sparse_status_t setUp();
	sparse_status_t precondition(double *x);
	sparse_status_t mv(double *x, double *y, struct matrix_descr descr);
};

#endif