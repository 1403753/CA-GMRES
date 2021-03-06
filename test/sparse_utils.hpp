#include "includes.hpp"

sparse_status_t extract(const Mtx_CSR *A,  // CSR input-Matrix. Column indices have to be in increasing order!
												Mtx_CSR *dest, // extracted CSR output-Matrix
												std::vector<size_t> &rows); // input rows

sparse_status_t extract(const Mtx_CSR *A,  // CSR input-Matrix
												Mtx_CSR *C, // extracted CSR output-Matrix
												std::vector<size_t> &rows, // input set
												std::vector<size_t> &cols);

sparse_status_t neighborhood(const Mtx_CSR *A,  // CSR MATRIX
						 std::vector<size_t> &v_In,  // input vector containing pointers to CSR rows
						 std::vector<size_t> &v_Res, // result vector for some level of order containing pointers to CSR rows
						 size_t order, // the reachabilty order
						 Gr_Part gPart
						);

sparse_status_t permute_Mtx(const Mtx_CSR *A, Mtx_CSR *dest, const size_t *pinv, const size_t *q);

sparse_status_t at_plus_a(
		const Mtx_CSR *A,
	  size_t *bnz,         // out - on exit, returns the actual number of nonzeros in matrix A'*A.
	  size_t **b_row_ptr,   // out - size n+1
	  size_t **b_col_indx    // out - size *bnz
);