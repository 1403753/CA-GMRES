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
	Mtx_CSR() : n{0}, nz{0}, row_ptr{nullptr}, col_indx{nullptr}, values{nullptr} {
		
	}
	Mtx_CSR(size_t n, size_t nz) : n{n}, nz{nz} {
		row_ptr = (size_t *) mkl_malloc((n + 1) * sizeof(size_t), 64);if(row_ptr == NULL){return;}	
		col_indx = (size_t *) mkl_malloc(nz * sizeof(size_t), 64);if(col_indx == NULL){return;}	
		values = (double *) mkl_malloc(nz * sizeof(double), 64);if(values == NULL){return;}	
	}
	virtual ~Mtx_CSR() {
		if (row_ptr != nullptr) {
			mkl_free(row_ptr);
			row_ptr = nullptr;
		}
		if (col_indx != nullptr) {
			mkl_free(col_indx);
			col_indx = nullptr;
		}
		if (values != nullptr) {
			mkl_free(values);
			values = nullptr;
		}
	};
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
	sparse_matrix_t						*A_mkl;
	sparse_matrix_t						*M_mkl;
	Mtx_CSR										*A_mtx;
	Mtx_CSR										*M_mtx;
	IKSPType 								  *kspType;
	IPCType										*pcType;
public:
	KSP();
	virtual ~KSP();

	void setOptions(size_t s, size_t t, double rTol, double aTol, double dTol, size_t maxit);
	void setOperator(sparse_matrix_t *A_mkl);
	void setPC(sparse_matrix_t *M_mkl);
	void setPCType(IPCType *pcType);
	void setKSPType(IKSPType *kspType);
	void solve(double *b, double *x);
	void setUp();

	size_t getMaxit() {return this->maxit;};
	size_t getS() {return this->s;};
	size_t getT() {return this->t;};
	sparse_matrix_t* getA_mkl() {return this->A_mkl;};
	sparse_matrix_t* getM_mkl() {return this->M_mkl;};
	Mtx_CSR* getA_mtx() {return this->A_mtx;};
	Mtx_CSR* getM_mtx() {return this->M_mtx;};
	void mkl_order(sparse_matrix_t *O);
};

#endif