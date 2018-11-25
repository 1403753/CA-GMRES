//============================================================================
// Name        : metisgraph.cpp
// Author      : Robert Ernstbrunner
// Version     :
// Copyright   : Your copyright notice
// Description : main in C++
//============================================================================

/**********/
/*  TODO  */
/**********/

/**************/
/*  END TODO  */
/**************/

/*********/
/* notes */
/*********/

#ifndef METIS_HPP

#include "sparse_utils.hpp"

#include <vector>
#include <algorithm>
#include <omp.h>

extern "C"
{
	#include <igraph.h>
}

// #include <cstddef> /* NULL */
#include <metis.h>

#define METIS_HPP

#define s 5

#define ScalarType double



int main(int argc, char *argv[]) {
		
	// A = [
  //    10     0     0     0    -2     0
  //     3     9     0     0     0     3
  //     0     7     8     7     0     0
  //     3     0     8    80     5     0
  //     0     8     0     9     9    13
  //     0     4     0     0     0    -1
	//     ]
	
	const size_t n_ = 6;
	const size_t nz_ = 18;

	////////////////
	// csr format //
	////////////////	
	
	ScalarType a[nz_] = {10,-2,3,9,3,7,8,7,3,8,80,5,8,9,9,13,4,-1};
	size_t rowptr_[n_+1] = {0,2,5,8,12,16,18};
	size_t colind_[nz_] = {0,4,0,1,5,1,2,3,0,2,3,4,1,3,4,5,1,5};
	
	////////////////
	// csc format //
	////////////////
	
	// ScalarType a[nz_] = {10,3,3,9,7,8,4,8,8,7,80,9,-2,5,9,2,3,13,-1};
	// size_t rowptr_[n_+1] = {0,3,7,9,12,16,19};
	// size_t colind_[nz_] = {0,1,3,1,2,4,5,2,3,2,3,4,0,3,4,5,1,4,5};

	
	size_t pinv[n_] = {0,1,2,3,4,5};
	
	Mtx_CSR R = {
		n_,
		rowptr_,
		colind_,
		a
	};
	
	size_t *Cp = new size_t[n_+1];
	size_t *Ci = new size_t[nz_];
	double *Cx = new double[nz_];
	
	Mtx_CSR C = {
		n_,
		Cp,
		Ci,
		Cx
	};
	
	my_permute(&R, &C, n_, pinv, pinv);
	
	// for (size_t i = 0; i < nz_; ++i)
		// std::cout << Ci[i] << "\n";
	// std::cout << std::endl;
	
	// std::cout << "permuted: \n"; 
	// for (size_t i = 0; i < nz_; ++i)
		// std::cout << Cx[i] << "\n";
	// std::cout << std::endl;

	// std::cout << "rowptr: \n"; 
	// for (size_t i = 0; i < n_+1; ++i)
		// std::cout << Cp[i] << "\n";
	// std::cout << std::endl;
	

	
	std::vector<size_t *> v_Res;
	std::vector<size_t *> v_In;
	size_t order = 1;
	
	v_In.push_back(Cp);
	v_In.push_back(Cp+1);
	// v_In.push_back(Cp+2);
	// v_In.push_back(Cp+3);
	// v_In.push_back(Cp+4);
	v_In.push_back(Cp+5);
	
	Gr_Part gPart = GRAPH_COMPLETE;
	
	neighborhood(&C,  // CSR MATRIX
							 v_In, // input set
							 v_Res, // result vertices for some level of order
							 order, // the reachabilty order
							 gPart
	);
	
	for(auto k:v_Res)
		std::cout << "RES:  " << std::distance(Cp, k) << "\n"; 

	delete[] Cp;
	delete[] Ci;
	delete[] Cx;	
	
	size_t bnz_;
	size_t *b_rowptr_;
	size_t *b_colind_;

	at_plus_a(
	  n_,      /* number of rows in matrix A. */
	  nz_,     /* number of nonzeros in matrix A */
	  rowptr_,      /* row pointer of size n+1 for matrix A. */
	  colind_,      /* column indices of size nz for matrix A. */
	  &bnz_,         /* out - on exit, returns the actual number of nonzeros in matrix A'*A. */
	  &b_rowptr_,   /* out - size n+1 */
	  &b_colind_    /* out - size *bnz */
	);
	
	size_t nWeights = 2;
	size_t nParts = 3;
	size_t options[METIS_NOPTIONS];
	size_t objval;
	size_t part[n_];
	struct matrix_descr descr;


	METIS_SetDefaultOptions((idx_t*) options);
	options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
	// options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB;
	options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL; // slower, but could effectively reduce overall communication volume
	// options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT; // faster, approximates minimum communication volume
	// options[METIS_OPTION_NUMBERING] = 1; // fortran
	
	// int METIS_PartGraphKway(
	// idx_t *nvtxs,   -> The number of vertices in the graph
	// idx_t *ncon,    -> The number of balancing constraints. It should be at least 1
	// idx_t *xadj,    -> The adjacency structure of the graph as described in Section 5.5
	// idx_t *adjncy,  -> The adjacency structure of the graph as described in Section 5.5
	// idx t *vwgt,    -> The weights of the vertices as described in Section 5.5
	// idx_t *vsize,   -> The size of the vertices for computing the total communication volume as described in Section 5.7
	// idx_t *adjwgt,  -> The weights of the edges as described in Section 5.5
	// idx_t *nparts,  -> The number of parts to partition the graph
	// real_t *tpwgts, -> This is an array of size nparts×ncon that specifies the desired weight for each partition and constraint.
	                   // The target partition weight for the ith partition and jth constraint is specified at tpwgts[i*ncon+j]
                     // (the numbering for both partitions and constraints starts from 0). For each constraint, the sum of the
                     // tpwgts[]entries must be 1.0 (i.e., \sum_i tpwgts[i*ncon+j] = 1.0)
										 // A NULL value can be passed to indicate that the graph should be equally divided among the partitions.
	// real_t ubvec,   -> This is an array of size ncon that speciﬁes the allowed load imbalance tolerance for each constraint.
                     // For the ith partition and jth constraint the allowed weight is the ubvec[j]*tpwgts[i*ncon+j] fraction
                     // of the jth’s constraint total weight. The load imbalances must be greater than 1.0.
                     // A NULLvalue can be passed indicating that the load imbalance tolerance for each constraint should
                     // be 1.001 (for ncon=1) or 1.01 (for ncon>1).
	// idx_t *options, -> This is the array of options as described in Section 5.4.
                     // The following options are valid for METIS PartGraphKway:
                     // METIS_OPTION_OBJTYPE, METIS_OPTION_CTYPE, METIS_OPTION_IPTYPE,
                     // METIS_OPTION_RTYPE, METIS_OPTION_NO2HOP, METIS_OPTION_NCUTS,
                     // METIS_OPTION_NITER, METIS_OPTION_UFACTOR, METIS_OPTION_MINCONN,
                     // METIS_OPTION_CONTIG, METIS_OPTION_SEED, METIS_OPTION_NUMBERING,
                     // METIS_OPTION_DBGLVL
	// idx_t *objval,  -> Upon successful completion, this variable stores the edge-cut or the total communication volume of
                     // the partitioning solution. The value returned depends on the partitioning’s objective function.
	// idx_t *part)    -> This is a vector of size nvtxs that upon successful completion stores the partition vector of the graph.
                     // The numbering of this vector starts from either 0 or 1, depending on the value of
                     // options[METIS_OPTION_NUMBERING].
	
	int ret = METIS_PartGraphKway((idx_t*) &n_, (idx_t*) &nWeights, (idx_t*) b_rowptr_, (idx_t*) b_colind_,
																NULL, NULL, NULL, (idx_t*) &nParts, NULL,
																NULL, (idx_t*) options, (idx_t*) &objval, (idx_t*) part);
	
	std::cout << ret << std::endl;
    
  for(size_t i = 0; i < n_; i++) {
		std::cout << i << ": " << part[i] << std::endl;
  }
	
	std::cout << "colindx:\n"; 
	for(size_t i = 0; i < bnz_; i++) {
		std::cout << b_colind_[i] << ", ";
  }
	std::cout << std::endl;

	std::cout << "rowptr:\n"; 
	for(size_t i = 0; i < n_+1; i++) {
		std::cout << b_rowptr_[i] << ", ";
  }
	std::cout << std::endl;
	
	mkl_free(b_rowptr_);
	if (bnz_) {
		mkl_free(b_colind_);
	}
	
	size_t ipar[128];		
	ipar[1] = 6; // 6 == display all error msgs on screen
	ipar[30] = 1; // use dpar[31] instead of 0 diagonal; set ipar[30] = 0 to abort computation instead

	double dpar[128];
	dpar[30] = 1.0e-16; // compare value to diagonal
	dpar[31] = 1.0e-10;	// the value that will be set on the diagonal
	
	///////////////
	// MKL ILU0  //
	///////////////
	
	// void dcsrilu0 (const MKL_INT *n,
								 // const double *a,			-> Array containing the set of elements of the matrix A
								 // const MKL_INT *ia,
								 // const MKL_INT *ja,
								 // double *bilu0,
								 // const MKL_INT *ipar,
								 // const double *dpar,
								 // MKL_INT *ierr );
	
	size_t *ia, *ja, ierr;
	ia = (size_t*) mkl_malloc((n_+1) * sizeof(size_t), 64); if(ia == NULL) {return 1;}
	ja = (size_t*) mkl_malloc(nz_ * sizeof(size_t), 64); if(ja == NULL) {return 1;}
	double *bilu0;
	bilu0 = (ScalarType*) mkl_malloc(nz_ * sizeof(ScalarType), 64); if(bilu0 == NULL) {return 1;}
	
	for(size_t i = 0; i < n_+1; ++i) {
		ia[i] = rowptr_[i] + 1;
	}
	
	for(size_t i = 0; i < nz_; ++i) {
		ja[i] = colind_[i] + 1;
	}
	
	dcsrilu0(&n_, a, ia, ja, bilu0, ipar, dpar, &ierr);

	for (size_t i = 0; i < nz_; ++i) {
		std::cout << bilu0[i] << ", ";
	}
	std::cout << std::endl;
	
	sparse_status_t   stat;
	sparse_matrix_t   A;
	sparse_matrix_t   B;
	
	// sparse_status_t mkl_sparse_d_create_csr (sparse_matrix_t *A,
																					 // sparse_index_base_t indexing,
																					 // MKL_INT rows,
																					 // MKL_INT cols,
																					 // MKL_INT *rows_start,
																					 // MKL_INT *rows_end,
																					 // MKL_INT *col_indx,
																					 // double *values);
	
	stat = mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO, n_, n_, rowptr_, rowptr_ + 1, colind_, a);
	stat = mkl_sparse_d_create_csr(&B, SPARSE_INDEX_BASE_ZERO, n_, n_, rowptr_, rowptr_ + 1, colind_, bilu0);

	ScalarType *x, *y;
	
	x = (ScalarType *) mkl_malloc(n_ * sizeof(ScalarType), 64);if(x == NULL){throw std::invalid_argument("Matrix Market Converter : malloc on 'x' failed.");}	
	y = (ScalarType *) mkl_malloc(n_ * sizeof(ScalarType), 64);if(y == NULL){throw std::invalid_argument("Matrix Market Converter : malloc on 'y' failed.");}	

	for (size_t i = 0; i < n_; ++i) {
		x[i] = 1;
	}
	
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;

	stat = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A, descr, x, 0, y);

	// sparse_status_t mkl_sparse_d_trsv (sparse_operation_t operation,
																		 // double alpha,
																		 // const sparse_matrix_t A,
																		 // struct matrix_descr descr,
																		 // const double *x,
																		 // double *y);
	descr.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
	descr.mode = SPARSE_FILL_MODE_LOWER;
	descr.diag = SPARSE_DIAG_UNIT;
	
	stat = mkl_sparse_d_trsv (SPARSE_OPERATION_NON_TRANSPOSE, 1, B, descr, y, y);

	descr.mode = SPARSE_FILL_MODE_UPPER;
	descr.diag = SPARSE_DIAG_NON_UNIT;

	stat = mkl_sparse_d_trsv (SPARSE_OPERATION_NON_TRANSPOSE, 1, B, descr, y, x);

	
	for (size_t i = 0; i < n_; ++i) {
		std::cout << " :" << x[i] << std::endl;
	}
	
	mkl_free(ia);
	mkl_free(ja);
	mkl_free(x);
	mkl_free(y);
	mkl_free(bilu0);
	mkl_sparse_destroy(A);
	mkl_sparse_destroy(B);
	mkl_free_buffers();
	
	return stat;
}


#endif