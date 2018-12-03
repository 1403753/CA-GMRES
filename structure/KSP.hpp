#ifndef KSP_HPP
#define KSP_HPP

#define MKL_MAX_PATH_LEN 4096
#define MKL_Complex16 std::complex<double>
#define MKL_INT size_t

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <complex>
#include <metis.h>
#include <omp.h>
// #include <exception>

#include "mkl.h"
#include "IKSPType.hpp"
#include "IPCType.hpp"

class IKSPType;
class IPCType;

struct Mtx_CSR{
	size_t n;
	size_t nz;
	size_t *row_ptr;
	size_t *col_indx;
	double *values;
	Mtx_CSR() {}
	Mtx_CSR(size_t n, size_t nz, size_t *row_ptr, size_t *col_indx, double *values)
		: n{n}, nz{nz}, row_ptr{row_ptr}, col_indx{col_indx}, values{values}
	{ }
};

typedef std::complex<double> complex_t; 
typedef std::pair<size_t, complex_t> ic_pair_t;

class KSP {
	double                    rTol = 1e-10;
	double                    aTol = 1e-50;
	double                    dTol = 1e+4;
	size_t                    maxit = 10000;
	size_t										s = 5;
	size_t										t = 12;
	Mtx_CSR										*A_mtx;
	Mtx_CSR										*M_mtx;
	IKSPType 								  *kspType;
	IPCType										*pcType;
public:
	KSP();
	void setOptions(size_t s, size_t t, double rTol, double aTol, double dTol, size_t maxit);
	void setOperators(Mtx_CSR *A_mtx, Mtx_CSR *M_mtx);
	void setPCType(IPCType *pcType);
	void setKSPType(IKSPType *kspType);
	void solve(double *b, double *x);
	void setUp();
	virtual ~KSP();
	size_t getMaxit() {return this->maxit;};
	size_t getS() {return this->s;};
	size_t getT() {return this->t;};
	Mtx_CSR* getA_mtx() {return this->A_mtx;};
	Mtx_CSR* getM_mtx() {return this->M_mtx;};
};

#endif