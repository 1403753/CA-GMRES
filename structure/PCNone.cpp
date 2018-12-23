#include "PCNone.hpp"

PCNone::PCNone() {

}

PCNone::~PCNone() {

}

sparse_status_t PCNone::mv(double *x, double *y, struct matrix_descr descr) {

	sparse_matrix_t *A_mkl = ksp->getA_mkl();

	mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, *A_mkl, descr, x, 0, y);

	return SPARSE_STATUS_SUCCESS;
}

sparse_status_t PCNone::precondition(double *x) {
	// do nothing
	return SPARSE_STATUS_SUCCESS;
}

sparse_status_t PCNone::setUp() {
	// do nothing
	return SPARSE_STATUS_SUCCESS;
}