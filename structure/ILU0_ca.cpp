#include "ILU0_ca.hpp"
#include "SparseUtils.hpp"
#include <metis.h>

#include <papi.h>

ILU0_ca::ILU0_ca() {

}

void ILU0_ca::test(){
	std::cout << "SUCCESS!" << std::endl;
}

sparse_status_t ILU0_ca::setUp() {

	size_t bnz, *b_row_ptr, *b_col_indx;
	Mtx_CSR *A_mtx = ksp->getA_mtx();
	Mtx_CSR *M_mtx = ksp->getM_mtx();
	size_t n = A_mtx->n;
	size_t nz = A_mtx->nz;
	Mtx_CSR destM, destA;
	ksp->createMtx(&destM, n, nz);
	ksp->createMtx(&destA, n, nz);

	sparse_matrix_t *A_mkl = ksp->getA_mkl();
	sparse_matrix_t *M_mkl = ksp->getM_mkl();
	size_t nWeights = 2; // <- set this to 1 eventually
	size_t nParts = n < (size_t) mkl_get_max_threads() ? n : (size_t) mkl_get_max_threads();
	size_t options[METIS_NOPTIONS];
	size_t objval;
	size_t *part;

	///////////////
	// MKL ILU0  //
	///////////////

	size_t ipar[128];		
	ipar[1] = 6; // 6 == display all error msgs on screen
	ipar[30] = 1; // use dpar[31] instead of 0 diagonal; set ipar[30] = 0 to abort computation instead

	double dpar[128];
	dpar[30] = 1.0e-16; // compare value to diagonal
	dpar[31] = 1.0e-10;	// the value that will be set on the diagonal

	// void dcsrilu0 (const MKL_INT *n,
								 // const double *a,			-> Array containing the set of elements of the matrix A
								 // const MKL_INT *ia,
								 // const MKL_INT *ja,
								 // double *bilu0,
								 // const MKL_INT *ipar,
								 // const double *dpar,
								 // MKL_INT *ierr );

	size_t *ia, *ja, ierr;
	ia = (size_t*) mkl_malloc((n + 1) * sizeof(size_t), 64); if(ia == NULL) {return SPARSE_STATUS_ALLOC_FAILED;}
	ja = (size_t*) mkl_malloc(nz * sizeof(size_t), 64); if(ja == NULL) {return SPARSE_STATUS_ALLOC_FAILED;}
	double *bilu0;
	bilu0 = (double*) mkl_malloc(nz * sizeof(double), 64); if(bilu0 == NULL) {return SPARSE_STATUS_ALLOC_FAILED;}

	for(size_t i = 0; i < n + 1; ++i) {
		ia[i] = A_mtx->row_ptr[i] + 1;
		// std::cout << A_mtx->row_ptr[i] << ", ";
	}
	// std::cout << std::endl;

	for(size_t i = 0; i < nz; ++i) {
		ja[i] = A_mtx->col_indx[i] + 1;
		// std::cout << A_mtx->col_indx[i] << ", ";		
	}
	// std::cout << std::endl;
	// std::cout << std::endl;	

	// for(size_t i = 0; i < nz; ++i) {
		// std::cout << A_mtx->values[i] << ", ";		
	// }
	// std::cout << std::endl;
	std::cout << "allocations finished\n";

	dcsrilu0(&n, A_mtx->values, ia, ja, bilu0, ipar, dpar, &ierr);		

	std::cout << "factorization finished\n";

	/////////////
	//  METIS  //
	/////////////

	part = (size_t *) mkl_malloc(n * sizeof(double), 64);if(part == NULL){return SPARSE_STATUS_ALLOC_FAILED;}	

	at_plus_a(A_mtx, &bnz, &b_row_ptr, &b_col_indx);

	METIS_SetDefaultOptions((idx_t*) options);
	options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;     // part k-way
	// options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB;       // recursive bisection
	options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;  // faster, approximates minimum communication volume
	// options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;  // slower, but could effectively reduce overall communication volume

	METIS_PartGraphKway((idx_t*) &n, (idx_t*) &nWeights, (idx_t*) b_row_ptr, (idx_t*) b_col_indx, NULL, NULL,
											NULL, (idx_t*) &nParts, NULL, NULL, (idx_t*) options, (idx_t*) &objval, (idx_t*) part);	

	std::cout << "metis finished\n";

	size_t incp = 0;

	for (size_t i = n; i-- > 0;) {
		part[i] = incp++;
	}

	M_mtx->n = n;
	M_mtx->nz = nz;
	M_mtx->values = bilu0;
	M_mtx->row_ptr = A_mtx->row_ptr;
	M_mtx->col_indx = A_mtx->col_indx;	

	permute_Mtx(M_mtx, &destM, part, part);


	mkl_sparse_d_create_csr(M_mkl, SPARSE_INDEX_BASE_ZERO, destM.n, destM.n, destM.row_ptr, destM.row_ptr + 1, destM.col_indx, destM.values);
	mkl_sparse_order(*M_mkl);
	ksp->setPC(M_mkl);


	permute_Mtx(A_mtx, &destA, part, part);	

	std::cout << "permutation finished\n";

	mkl_sparse_destroy(*A_mkl);

	mkl_sparse_d_create_csr(A_mkl, SPARSE_INDEX_BASE_ZERO, destA.n, destA.n, destA.row_ptr, destA.row_ptr + 1, destA.col_indx, destA.values);
	mkl_sparse_order(*A_mkl);
	ksp->setOperator(A_mkl);

	// compute structure at_plus_a(A_mtx) ->output: b_rowptr and partition A with metis 
	// sort A_mtx with amml(s)
	// permute A_mtx -> col_indx out of order
	// create A_mkl and sort col_indx (maybe create own algorithm)
	// export A_mtx from A_mkl again
	// factorize A_mtx to get ILU(0)-values
	// find subsets alpha_p p == #threads depending on matrix size
	// openmp: find s-step dependencies for each alpha_p -> beta_p -> gamma_p -> delta_p 

	mkl_free(ia);
	mkl_free(ja);
	mkl_free(bilu0);
	mkl_free(part);
	mkl_free(b_row_ptr);
	mkl_free(b_col_indx);

	return SPARSE_STATUS_SUCCESS;
}
	////////////////////
	//  METIS OUTPUT  //
	////////////////////

	// size_t k = -1;
	// std::cout << "\nMetis K-Way:" << std::endl;
	// while (++k != (size_t) mkl_get_max_threads()) {
		// for(size_t i = 0; i < n; i++) {
			// if (part[i] == k) {
				// std::cout << i << ": " << part[i] << std::endl;
				// break;
			// }
		// }
	// }


	/////////////////////
	//  PAPI TEMPLATE  //
	/////////////////////

	// float                         rtime, ptime, mflops;
	// long long                     flpops;		
	// if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
		// exit(1);
	// // -> function here <-
	// if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
		// exit(1);
	// PAPI_shutdown();
	// std::cout << "runtime: " << rtime << " seconds." << std::endl;