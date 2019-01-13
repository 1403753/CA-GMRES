#ifndef KSP__HPP
#define KSP__HPP

#define MKL_MAX_PATH_LEN 4096
#define MKL_Complex16 std::complex<double>
#define MKL_INT size_t

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <complex>

#include "mkl.h"
#include "IKSPType.hpp"
#include "IPCType.hpp"

#include <memory>


class IKSPType;
class IPCType;

struct Mtx_CSR{
	size_t n;
	size_t nz;
	size_t *row_ptr;
	size_t *col_indx;
	double *values;
	Mtx_CSR() : n{0}, nz{0}, row_ptr{nullptr}, col_indx{nullptr}, values{nullptr}{}
	virtual ~Mtx_CSR(){}
};

typedef std::complex<double>          complex_t; 
typedef std::pair<size_t, complex_t>  ic_pair_t;

class KSP_ {
	double                                  rTol = 1e-10;
	double                                  aTol = 1e-50;
	double                                  dTol = 1e+4;
	size_t                                  maxit = 1000;
	sparse_matrix_t						              *A_mkl;
	sparse_matrix_t						              M_mkl;
	const std::shared_ptr<Mtx_CSR>          A_ptr;
	const std::shared_ptr<Mtx_CSR>          M_ptr;
	IKSPType 								                *kspType;
	IPCType										              *pcType;
	bool                                    storeHist;
	std::vector<std::pair<size_t, double>>  rHist;	
public:
	KSP_();
	virtual ~KSP_();

	sparse_status_t setOptions(double rTol, double aTol, double dTol, size_t maxit, bool storeHist);
	sparse_status_t setOperator(sparse_matrix_t *A_mkl);
	sparse_status_t setPC(sparse_matrix_t *M_mkl);
	sparse_status_t setPCType(IPCType *pcType);
	sparse_status_t setKSPType(IKSPType *kspType);
	sparse_status_t solve(double *b, double *x);
	sparse_status_t setUp();

	size_t getMaxit() {return this->maxit;};
	double getRTol() {return this->rTol;};
	double getDTol() {return this->dTol;};
	bool getStoreHist() {return this->storeHist;};
	std::vector<std::pair<size_t, double>>* getRHist() {return &this->rHist;};
	sparse_matrix_t* getA_mkl() {return this->A_mkl;};
	sparse_matrix_t* getM_mkl() {return &this->M_mkl;};
	const std::shared_ptr<Mtx_CSR> getA_ptr() {return this->A_ptr;};
	const std::shared_ptr<Mtx_CSR> getM_ptr() {return this->M_ptr;};
	IPCType* getIPCType() {return this->pcType;};
	void print(Mtx_CSR *P);
	void createMtx(Mtx_CSR *Mtx, size_t n, size_t nz);
	void destroyMtx(Mtx_CSR *Mtx);
};

#endif