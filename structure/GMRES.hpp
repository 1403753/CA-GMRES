#ifndef GMRES_HPP
#define GMRES_HPP

#include "IKSPType.hpp"
#include "PCILU0_ca.hpp"
#include <papi.h>

class GMRES : public IKSPType {
	size_t										       m = 60;
	float                            SpMV = 0;
	float                            MGS = 0;
	float                            rtime, ptime, mflops;
	long long                        flpops;	
public:
	GMRES();
	GMRES(size_t m);
	virtual ~GMRES();
	sparse_status_t solve(double *x, double *b);
	float getSpMV() {return this->SpMV;};
	float getMGS() {return this->MGS;};
};

#endif