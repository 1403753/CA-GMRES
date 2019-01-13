#ifndef PCNONE_HPP
#define PCNONE_HPP

#include "IPCType.hpp"

class PCNone : public IPCType{
public:
	PCNone();
	virtual ~PCNone();
	sparse_status_t setUp();
	sparse_status_t precondition(double *x);
	sparse_status_t mv(double *x, double *y, struct matrix_descr descr);
};

#endif