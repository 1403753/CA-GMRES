#include "PCILU0_ca.hpp"

PCILU0_ca::PCILU0_ca() {

}

PCILU0_ca::~PCILU0_ca() {
	if (setup) {
		const	std::shared_ptr<Mtx_CSR> M_ptr = ksp->getM_ptr();
		ksp->destroyMtx(M_ptr.get());
		mkl_sparse_destroy(*this->M_mkl);
	}
}

sparse_status_t PCILU0_ca::mv(double *x, double *y, struct matrix_descr descr) {

	sparse_matrix_t *A_mkl = ksp->getA_mkl();

	mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, *A_mkl, descr, x, 0, y);
	precondition(y);

	return SPARSE_STATUS_SUCCESS;
}

sparse_status_t PCILU0_ca::precondition(double *x) {

	struct matrix_descr descr;

	descr.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
	descr.mode = SPARSE_FILL_MODE_LOWER;
	descr.diag = SPARSE_DIAG_UNIT;
	
	sparse_matrix_t *M_mkl = ksp->getM_mkl();
	sparse_status_t stat;
	stat = mkl_sparse_d_trsv (SPARSE_OPERATION_NON_TRANSPOSE, 1, *M_mkl, descr, x, x);
	
	descr.mode = SPARSE_FILL_MODE_UPPER;
	descr.diag = SPARSE_DIAG_NON_UNIT;
	
	stat = mkl_sparse_d_trsv (SPARSE_OPERATION_NON_TRANSPOSE, 1, *M_mkl, descr, x, x);		
	
	return stat;
}

sparse_status_t PCILU0_ca::setUp() {
	this->setup = true;
	const std::shared_ptr<Mtx_CSR> A_ptr = ksp->getA_ptr();

	size_t n = A_ptr->n;
	size_t nz = A_ptr->nz;

	const	std::shared_ptr<Mtx_CSR> M_ptr = ksp->getM_ptr();
	ksp->createMtx(M_ptr.get(), n, nz);
	
	M_mkl = ksp->getM_mkl();
	
	for (size_t i = 0; i < nz; ++i) {
		M_ptr->col_indx[i] = A_ptr->col_indx[i];
	}

	for (size_t i = 0; i < n + 1; ++i) {
		M_ptr->row_ptr[i] = A_ptr->row_ptr[i];
	}

	///////////////
	// MKL ILU0  //
	///////////////

	size_t ipar[128];
	ipar[1] = 6; // 6 == display all error msgs on screen
	ipar[30] = 1; // use dpar[31] instead of 0 diagonal; set ipar[30] = 0 to abort computation instead

	double dpar[128];
	dpar[30] = 1.0e-16; // compare value to diagonal
	dpar[31] = 1.0e-10;	// the value that will be set on the diagonal

	size_t *ia, *ja, ierr = 0;
	ia = (size_t*) mkl_malloc((n + 1) * sizeof(size_t), 64); if(ia == NULL) {return SPARSE_STATUS_ALLOC_FAILED;}
	ja = (size_t*) mkl_malloc(nz * sizeof(size_t), 64); if(ja == NULL) {return SPARSE_STATUS_ALLOC_FAILED;}

	for(size_t i = 0; i < n + 1; ++i) {
		ia[i] = A_ptr->row_ptr[i] + 1;
	}
	
	for(size_t i = 0; i < nz; ++i) {
		ja[i] = A_ptr->col_indx[i] + 1;
	}
	
	dcsrilu0(&n, A_ptr->values, ia, ja, M_ptr->values, ipar, dpar, &ierr);
	if (ierr != 0)
		throw std::invalid_argument("couldn't compute preconditioner");
	std::cout << "factorization finished, status: " << ierr << std::endl;	

	mkl_sparse_d_create_csr(M_mkl, SPARSE_INDEX_BASE_ZERO, M_ptr->n, M_ptr->n, M_ptr->row_ptr, M_ptr->row_ptr + 1, M_ptr->col_indx, M_ptr->values);
	
	ksp->setPC(M_mkl);

	mkl_free(ia);
	mkl_free(ja);

	return SPARSE_STATUS_SUCCESS;
}