#include "KSP.hpp"

KSP::KSP() : kspType{nullptr}, pcType{nullptr} {}

void KSP::setOperators(Mtx_CSR *A_mtx, Mtx_CSR *M_mtx) {
	this->A_mtx = A_mtx;
	this->M_mtx = M_mtx;
}

void KSP::setOptions(size_t s, size_t t, double rTol, double aTol, double dTol, size_t maxit) {
	this->s = s;
	this->t = t;
	this->rTol = rTol;
	this->aTol = aTol;
	this->dTol = dTol;   
	this->maxit	= maxit;  
}

void KSP::setKSPType(IKSPType *kspType) {
	this->kspType = kspType;
	this->kspType->ctx = this;
}

void KSP::solve(double *b, double *x) {
	if (kspType)
		this->kspType->solve(b, x);
}

void KSP::setPCType(IPCType *pcType) {
	this->pcType = pcType;
	this->pcType->ctx = this;	
}

void KSP::setUp() {
	if (pcType)
		this->pcType->setUp();
}

KSP::~KSP() {
	mkl_free(A_mtx->values);
	mkl_free(A_mtx->row_ptr);
	mkl_free(A_mtx->col_indx);
	A_mtx->values = nullptr;
	A_mtx->row_ptr = nullptr;
	A_mtx->col_indx = nullptr;
	if (M_mtx->values == nullptr) {
		mkl_free(M_mtx->values);
		M_mtx->values = nullptr;
	}
	if (M_mtx->row_ptr == nullptr) {
		mkl_free(M_mtx->row_ptr);
		M_mtx->row_ptr = nullptr;
	}
	if (M_mtx->col_indx == nullptr) {
		mkl_free(M_mtx->col_indx);
		M_mtx->col_indx = nullptr;
	}
}