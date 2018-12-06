#ifndef ILU0_CA_HPP
#define ILU0_CA_HPP

#include "IPCType.hpp"

class ILU0_ca : public IPCType{
public:
	void test();
	ILU0_ca();
	sparse_status_t setUp();
};

#endif