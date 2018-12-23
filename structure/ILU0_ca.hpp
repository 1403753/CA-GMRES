#ifndef ILU0_CA_HPP
#define ILU0_CA_HPP

#include "IPCType.hpp"

class ILU0_ca : public IPCType{
public:
	ILU0_ca();
	virtual ~ILU0_ca();
	sparse_status_t setUp();
	sparse_status_t precondition(double *x);
	sparse_status_t mv(double *x, double *y, struct matrix_descr descr);
};

#endif