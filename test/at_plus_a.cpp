#include "at_plus_a.hpp"

int at_plus_a(
	  const size_t n,      /* number of columns in matrix A. */
	  const size_t nz,     /* number of nonzeros in matrix A */
	  size_t *colptr,      /* column pointer of size n+1 for matrix A. */
	  size_t *rowind,      /* row indices of size nz for matrix A. */
	  size_t *bnz,         /* out - on exit, returns the actual number of
                               nonzeros in matrix A'*A. */
	  size_t **b_colptr,   /* out - size n+1 */
	  size_t **b_rowind    /* out - size *bnz */
	  )
{
/*
 * Purpose
 * =======
 *
 * Form the structure of A'+A. A is an n-by-n matrix in column oriented
 * format represented by (colptr, rowind). The output A'+A is in column
 * oriented format (symmetrically, also row oriented), represented by
 * (b_colptr, b_rowind).
 *
 */
	size_t i, j, k, col, num_nz;
	size_t *t_colptr, *t_rowind; /* a column oriented form of T = A' */
	size_t *marker;

	marker = (size_t*) mkl_malloc(n * sizeof(size_t), 64); if(marker == NULL) {return 1;}
	t_colptr = (size_t*) mkl_malloc((n+1) * sizeof(size_t), 64); if(t_colptr == NULL) {return 1;}
	t_rowind = (size_t*) mkl_malloc(nz * sizeof(size_t), 64); if(t_rowind == NULL) {return 1;}
	
	/* Get counts of each column of T, and set up column pointers */
	for (i = 0; i < n; ++i) marker[i] = 0;
	for (j = 0; j < n; ++j) {
		for (i = colptr[j]; i < colptr[j+1]; ++i)
			++marker[rowind[i]];
	}
	t_colptr[0] = 0;
	for (i = 0; i < n; ++i) {
		t_colptr[i+1] = t_colptr[i] + marker[i];
		marker[i] = t_colptr[i];
	}

	/* Transpose the matrix from A to T */
	for (j = 0; j < n; ++j)
		for (i = colptr[j]; i < colptr[j+1]; ++i) {
			col = rowind[i];
			t_rowind[marker[col]] = j;
			++marker[col];
		}


	/* ----------------------------------------------------------------
		 compute B = A + T, where column j of B is:

		 Struct (B_*j) = Struct (A_*k) UNION Struct (T_*k)

		 do not include the diagonal entry
		 ---------------------------------------------------------------- */

	/* Zero the diagonal flag */
	for (i = 0; i < n; ++i) marker[i] = -1;

	/* First pass determines number of nonzeros in B */
	num_nz = 0;
	for (j = 0; j < n; ++j) {
		/* Flag the diagonal so it's not included in the B matrix */
		marker[j] = j;

		/* Add pattern of column A_*k to B_*j */
		for (i = colptr[j]; i < colptr[j+1]; ++i) {
			k = rowind[i];
			if ( marker[k] != j ) {
				marker[k] = j;
				++num_nz;
			}
		}

		/* Add pattern of column T_*k to B_*j */
		for (i = t_colptr[j]; i < t_colptr[j+1]; ++i) {
			k = t_rowind[i];
			if ( marker[k] != j ) {
				marker[k] = j;
				++num_nz;
			}
		}
	}
		
	*bnz = num_nz;
  
	/* Allocate storage for A+A' */
	*b_colptr = (size_t*) mkl_malloc((n+1) * sizeof(size_t), 64); if(b_colptr == NULL) {return 1;}
	if ( *bnz) {
		*b_rowind = (size_t*) mkl_malloc(*bnz * sizeof(size_t), 64); if(b_rowind == NULL) {return 1;}
	}
    
	/* Zero the diagonal flag */
	for (i = 0; i < n; ++i) marker[i] = -1;
    
	/* Compute each column of B, one at a time */
	num_nz = 0;
	for (j = 0; j < n; ++j) {
		(*b_colptr)[j] = num_nz;
	
		/* Flag the diagonal so it's not included in the B matrix */
		marker[j] = j;

		/* Add pattern of column A_*k to B_*j */
		for (i = colptr[j]; i < colptr[j+1]; ++i) {
				k = rowind[i];
				if ( marker[k] != j ) {
			marker[k] = j;
			(*b_rowind)[num_nz++] = k;
				}
		}

		/* Add pattern of column T_*k to B_*j */
		for (i = t_colptr[j]; i < t_colptr[j+1]; ++i) {
			k = t_rowind[i];
			if ( marker[k] != j ) {
				marker[k] = j;
				(*b_rowind)[num_nz++] = k;
			}
		}
	}
	(*b_colptr)[n] = num_nz;
		 
	mkl_free(marker);
	mkl_free(t_colptr);
	mkl_free(t_rowind);	
	
	return 0;
}