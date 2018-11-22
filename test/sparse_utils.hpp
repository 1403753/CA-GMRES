#include "includes.hpp"

sparse_status_t my_permute(const RobMat *A, RobMat *dest, const size_t n, size_t *pinv, size_t *q);
// double *a, size_t *rowptr_, size_t *colind_, 

sparse_status_t at_plus_a(
	  const size_t n,      /* number of columns in matrix A. */
	  const size_t nz,     /* number of nonzeros in matrix A */
	  size_t *colptr,      /* column pointer of size n+1 for matrix A. */
	  size_t *rowind,      /* row indices of size nz for matrix A. */
	  size_t *bnz,         /* out - on exit, returns the actual number of
                               nonzeros in matrix A'*A. */
	  size_t **b_colptr,   /* out - size n+1 */
	  size_t **b_rowind    /* out - size *bnz */
);