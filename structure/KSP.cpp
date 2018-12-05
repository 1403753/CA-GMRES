#include "KSP.hpp"

KSP::KSP() : A_mkl{nullptr}, M_mkl{nullptr}, A_mtx{new Mtx_CSR()}, M_mtx{new Mtx_CSR()}, kspType{nullptr}, pcType{nullptr} {}

void KSP::setOperator(sparse_matrix_t *A_mkl) {
	sparse_index_base_t indexing;
	size_t n, m;
	size_t *rows_start;
	size_t *rows_end;
	size_t *col_indx;
	double *values;
	
	this->A_mkl = A_mkl;
	mkl_sparse_d_export_csr(*this->A_mkl, &indexing, &n, &m, &rows_start, &rows_end, &col_indx, &values);
	
	this->A_mtx->n = n;
	this->A_mtx->row_ptr = rows_start;
	this->A_mtx->nz = rows_start[n];	
	this->A_mtx->col_indx = col_indx;
	this->A_mtx->values = values;
}

void KSP::setPC(sparse_matrix_t *M_mkl) {
	sparse_index_base_t indexing;
	size_t n, m;
	size_t *rows_start;
	size_t *rows_end;
	size_t *col_indx;
	double *values;
	
	this->M_mkl = M_mkl;
	mkl_sparse_d_export_csr(*this->M_mkl, &indexing, &n, &m, &rows_start, &rows_end, &col_indx, &values);
	
	this->M_mtx->n = n;
	this->M_mtx->row_ptr = rows_start;
	this->M_mtx->nz = rows_start[n];	
	this->M_mtx->col_indx = col_indx;
	this->M_mtx->values = values;
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

void KSP::mkl_order(sparse_matrix_t *O) {
	mkl_sparse_order(*O);	
}

KSP::~KSP() {
	this->A_mtx->row_ptr = nullptr;
	this->A_mtx->col_indx = nullptr;
	this->A_mtx->values = nullptr;	
	delete this->A_mtx;
	
	this->M_mtx->row_ptr = nullptr;
	this->M_mtx->col_indx = nullptr;
	this->M_mtx->values = nullptr;
	delete this->M_mtx;


}