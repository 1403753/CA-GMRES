#include "includes.hpp"

int at_plus_a(
	  const size_t n,      /* number of columns in matrix A. */
	  const size_t nz,     /* number of nonzeros in matrix A */
	  size_t *colptr,      /* column pointer of size n+1 for matrix A. */
	  size_t *rowind,      /* row indices of size nz for matrix A. */
	  size_t *bnz,         /* out - on exit, returns the actual number of
                               nonzeros in matrix A'*A. */
	  size_t **b_colptr,   /* out - size n+1 */
	  size_t **b_rowind    /* out - size *bnz */
);