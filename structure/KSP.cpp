#include "KSP.hpp"

KSP::KSP() : A_mkl{nullptr}, M_mkl{nullptr}, A_ptr{std::make_shared<Mtx_CSR>()}, M_ptr{std::make_shared<Mtx_CSR>()}, kspType{nullptr}, pcType{nullptr} { }

void KSP::print(Mtx_CSR *P) {
size_t kp = 0;
size_t nz = P->nz;
	for (size_t i = 0; i < nz; ++i) {
		if (i == P->row_ptr[kp]) {
			std::cout << P->row_ptr[kp++] << "\t: ";
		} else {
			std::cout << "\t: ";
		}
		while (i == P->row_ptr[kp]) {
			std::cout << std::endl << P->row_ptr[kp++] << "\t: ";
		}
		std::cout << P->col_indx[i] << "; " << P->values[i] << " ><";
		std::cout << std::endl;
	}
	std::cout << P->row_ptr[kp] << "  =>  ";	
	std::cout << P->nz << " nz" << std::endl;
}

void KSP::createMtx(Mtx_CSR *Mtx, size_t n, size_t nz) {
	Mtx->n = n;
	Mtx->nz = nz;
	Mtx->row_ptr = (size_t *) mkl_malloc((n + 1) * sizeof(size_t), 64);if(Mtx->row_ptr == NULL){return;}
	Mtx->col_indx = (size_t *) mkl_malloc(nz * sizeof(size_t), 64);if(Mtx->col_indx == NULL){return;}
	Mtx->values = (double *) mkl_malloc(nz * sizeof(double), 64);if(Mtx->values == NULL){return;}
}

void KSP::destroyMtx(Mtx_CSR *Mtx) {
	if(Mtx) {
		if (Mtx->row_ptr != nullptr) {
			mkl_free(Mtx->row_ptr);
			Mtx->row_ptr = nullptr;
		}
		if (Mtx->col_indx != nullptr) {
			mkl_free(Mtx->col_indx);
			Mtx->col_indx = nullptr;
		}
		if (Mtx->values != nullptr) {
			mkl_free(Mtx->values);
			Mtx->values = nullptr;
		}
	}
}

sparse_status_t KSP::setOperator(sparse_matrix_t *A_mkl) {
	sparse_status_t stat;
	sparse_index_base_t indexing;
	size_t n, m;
	size_t *rows_start;
	size_t *rows_end;
	size_t *col_indx;
	double *values;
	
	this->A_mkl = A_mkl;
	stat = mkl_sparse_d_export_csr(*this->A_mkl, &indexing, &n, &m, &rows_start, &rows_end, &col_indx, &values);

	this->A_ptr->n = n;
	this->A_ptr->row_ptr = rows_start;
	this->A_ptr->nz = rows_start[n];
	this->A_ptr->col_indx = col_indx;
	this->A_ptr->values = values;
	
	return stat;
}

sparse_status_t KSP::setPC(sparse_matrix_t *M_mkl) {
	this->M_mkl = *M_mkl;
	return SPARSE_STATUS_SUCCESS;
}

sparse_status_t KSP::setOptions(double rTol, double aTol, double dTol, size_t maxit) {
	this->rTol = rTol;
	this->aTol = aTol;
	this->dTol = dTol;   
	this->maxit	= maxit;  
	return SPARSE_STATUS_SUCCESS;
}

sparse_status_t KSP::setKSPType(IKSPType *kspType) {
	this->kspType = kspType;
	this->kspType->ksp = this;
	return SPARSE_STATUS_SUCCESS;
}

sparse_status_t KSP::solve(double *b, double *x) {
	sparse_status_t stat = SPARSE_STATUS_EXECUTION_FAILED;
	if (kspType)
		stat = this->kspType->solve(b, x);
	return stat;
}

sparse_status_t KSP::setPCType(IPCType *pcType) {
	this->pcType = pcType;
	this->pcType->ksp = this;
	return SPARSE_STATUS_SUCCESS;
}

sparse_status_t KSP::setUp() {
	sparse_status_t stat = SPARSE_STATUS_EXECUTION_FAILED;
	if (pcType)
		stat = this->pcType->setUp();
	return stat;
}

KSP::~KSP() {
	mkl_sparse_destroy(M_mkl);
	destroyMtx(this->M_ptr.get());
}