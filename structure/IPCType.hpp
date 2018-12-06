#ifndef IPCTYPE_HPP
#define IPCTYPE_HPP

#include "KSP.hpp"

struct Mtx_CSR;
class KSP;

class IPCType {
public:
	KSP *ksp;
	virtual sparse_status_t setUp() = 0;
};


#endif