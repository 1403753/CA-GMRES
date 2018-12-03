#ifndef ILU0_CA_HPP
#define ILU0_CA_HPP

#include "IPCType.hpp"

enum Gr_Part {GRAPH_LOWER, GRAPH_UPPER, GRAPH_COMPLETE};

class ILU0_ca : public IPCType{
public:
	ILU0_ca();
	sparse_status_t setUp();
};

#endif