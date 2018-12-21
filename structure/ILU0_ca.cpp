#include "ILU0_ca.hpp"
#include "SparseUtils.hpp"
#include <metis.h>

#include <papi.h>

ILU0_ca::ILU0_ca() {

}

ILU0_ca::~ILU0_ca() {

}

sparse_status_t ILU0_ca::mpk(double *x, double *y, size_t s) {
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;

	return SPARSE_STATUS_SUCCESS;
};

sparse_status_t ILU0_ca::precondition(double *x) {

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
	
	return SPARSE_STATUS_SUCCESS;
};

sparse_status_t ILU0_ca::setUp() {

	size_t bnz, *b_row_ptr, *b_col_indx;
	Mtx_CSR *A_mtx = ksp->getA_mtx();
	size_t n = A_mtx->n;
	size_t nz = A_mtx->nz;
	size_t s = ksp->getS();	
	Mtx_CSR *M_mtx = ksp->getM_mtx();
	ksp->createMtx(M_mtx, n, nz);

	const std::shared_ptr<Mtx_CSR> M_ptr = ksp->getM_ptr();
	ksp->createMtx(M_ptr.get(), n, nz);
	
	sparse_matrix_t *A_mkl = ksp->getA_mkl();
	sparse_matrix_t *M_mkl = ksp->getM_mkl();
	// size_t nWeights = 2; // <- set this to 1 eventually
	// this->nParts = n < (size_t) mkl_get_max_threads() ? 2 : (size_t) mkl_get_max_threads();
	size_t options[METIS_NOPTIONS];
	// size_t objval;
	// size_t *part;
	// double *P_indx;

	/////////////
	//  METIS  //
	/////////////

	// part = (size_t *) mkl_malloc(n * sizeof(double), 64);if(part == NULL){return SPARSE_STATUS_ALLOC_FAILED;}	

	at_plus_a(A_mtx, &bnz, &b_row_ptr, &b_col_indx);

	METIS_SetDefaultOptions((idx_t*) options);
	// options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;     // part k-way
	// options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB;       // recursive bisection
	// options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;  // faster, approximates minimum communication volume
	// options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;  // slower, but could effectively reduce overall communication volume

	// METIS_PartGraphKway((idx_t*) &n, (idx_t*) &nWeights, (idx_t*) b_row_ptr, (idx_t*) b_col_indx, NULL, NULL,
											// NULL, (idx_t*) &this->nParts, NULL, NULL, (idx_t*) options, (idx_t*) &objval, (idx_t*) part);	


	
	size_t *perm;
	size_t *iperm;
	perm = (size_t *) mkl_malloc(n * sizeof(double), 64);if(perm == NULL){return SPARSE_STATUS_ALLOC_FAILED;}		
	iperm = (size_t *) mkl_malloc(n * sizeof(double), 64);if(iperm == NULL){return SPARSE_STATUS_ALLOC_FAILED;}		
	
	// METIS NodeND(idx t*nvtxs, idx t*xadj, idx t*adjncy, idx t*vwgt, idx t*options, idx t*perm, idx t*iperm)	
	METIS_NodeND((idx_t*) &n, (idx_t*) b_row_ptr, (idx_t*) b_col_indx, NULL, NULL, (idx_t*) perm, (idx_t*) iperm);
	
	for(size_t i = 0; i < n; ++i)
		perm[i] = iperm[i] = i;
		
	permute_Mtx(A_mtx, M_ptr.get(), perm, iperm);	
	
	std::cout << "metis finished\n";

	// mkl_sparse_destroy(*A_mkl);

		///////////////
	// MKL ILU0  //
	///////////////

	size_t ipar[128];		
	ipar[1] = 6; // 6 == display all error msgs on screen
	ipar[30] = 1; // use dpar[31] instead of 0 diagonal; set ipar[30] = 0 to abort computation instead

	double dpar[128];
	dpar[30] = 1.0e-16; // compare value to diagonal
	dpar[31] = 1.0e-10;	// the value that will be set on the diagonal



	size_t *ia, *ja, ierr;
	ia = (size_t*) mkl_malloc((n + 1) * sizeof(size_t), 64); if(ia == NULL) {return SPARSE_STATUS_ALLOC_FAILED;}
	ja = (size_t*) mkl_malloc(nz * sizeof(size_t), 64); if(ja == NULL) {return SPARSE_STATUS_ALLOC_FAILED;}

	//////////////////////
	//  test vec start  //
	//////////////////////
	
	// double *vec, *vec2;	
	// vec = (double*) mkl_malloc(n * sizeof(double), 64); if(vec == NULL) {return SPARSE_STATUS_ALLOC_FAILED;}
	// vec2 = (double*) mkl_malloc(n * sizeof(double), 64); if(vec2 == NULL) {return SPARSE_STATUS_ALLOC_FAILED;}

	// for(size_t i = 0; i < n; ++i) {
		// vec[i] = i;
	// }
	
	// permute_Vec(vec, perm, n, vec2);
	// permute_Vec(vec, perm, n, NULL);
	
	// for(size_t i = 0; i < n; ++i) {
		// std::cout << vec2[i] << " <<\n";
	// }	
	// std::cout << "--------------\n";
	// for(size_t i = 0; i < n; ++i) {
		// std::cout << vec2[i] << " <<\n";
	// }	

	// mkl_free(vec);
	// mkl_free(vec2);

	////////////////////
	//  test vec end  //
	////////////////////	

	for(size_t i = 0; i < n + 1; ++i) {
		ia[i] = A_mtx->row_ptr[i] + 1;
	}
	
	for(size_t i = 0; i < nz; ++i) {
		ja[i] = A_mtx->col_indx[i] + 1;
	}
	
	// void dcsrilu0 (const MKL_INT *n,
								 // const double *a,			-> Array containing the set of elements of the matrix A
								 // const MKL_INT *ia,
								 // const MKL_INT *ja,
								 // double *bilu0,
								 // const MKL_INT *ipar,
								 // const double *dpar,
								 // MKL_INT *ierr );
	
	dcsrilu0(&n, A_mtx->values, ia, ja, M_ptr->values, ipar, dpar, &ierr);
	
	std::cout << "factorization finished\n";
	
	//////////////////////////////////////

	
	mkl_sparse_d_create_csr(M_mkl, SPARSE_INDEX_BASE_ZERO, M_ptr->n, M_mtx->n, M_ptr->row_ptr, M_ptr->row_ptr + 1, M_ptr->col_indx, M_ptr->values);
	
	mkl_sparse_order(*M_mkl);

	ksp->setPC(M_mkl);

	// sparse_index_base_t indexingA;
	// size_t nA, mA;
	// size_t *rows_startA;
	// size_t *rows_endA;
	// size_t *col_indxA;
	// double *valuesA;
	
	// mkl_sparse_d_export_csr(*A_mkl, &indexingA, &nA, &mA, &rows_startA, &rows_endA, &col_indxA, &valuesA);
	
	// for (size_t i = 0; i < nz; ++i) {
		// size_t idx = (size_t)valuesA[i];
		// A_mtx->ilu0_values[i] = destA.ilu0_values[idx];
		// valuesA[i] = destA.values[idx];
	// }
	
	// mkl_free(destA.values);
	// mkl_free(destA.ilu0_values);
	
	// ksp->setOperator(A_mkl);

	// this->alphaArr.resize(nParts);
	// std::vector<std::vector<std::vector<size_t>>> betaArr(nParts, std::vector<std::vector<size_t>>(ksp->getS()));
	// std::vector<std::vector<std::vector<size_t>>> gammaArr(nParts, std::vector<std::vector<size_t>>(ksp->getS()));
	// std::vector<std::vector<std::vector<size_t>>> deltaArr(nParts, std::vector<std::vector<size_t>>(ksp->getS()));
	
	// size_t alphaPartSize = (size_t)((n / nParts) * 1.03 + 5); // add 103% + 5 buffer
	
	// for (size_t i = 0; i < nParts; ++i) {
		// alphaArr.at(i).reserve(alphaPartSize); 
	// }
	
	// for (size_t i = 0; i < n; ++i) {
		// alphaArr.at(part[i]).push_back(i);
	// }
	
	// for (size_t i = 0; i < nParts; ++i) {
		// alphaArr.at(i).shrink_to_fit();
	// }	
	
	// size_t sum = 0;
	
	// for (auto a:alphaArr) {
		// sum += a.size();
		// std::cout << a.size() << std::endl;
		// for(size_t i = 0; i < a.size(); ++i)
				// std::cout << a.at(i) << ", ";
		// std::cout << std::endl << std::endl;
	// }
	
	// std::cout << std::endl << "reserved: " << alphaPartSize << std::endl;
	// std::cout << std::endl << "sum: " << sum << std::endl;

	// #pragma omp parallel for
	// for (size_t i = 0; i < nParts; ++i)
		// find(A_mtx, alphaArr.at(i), betaArr.at(i), gammaArr.at(i), deltaArr.at(i));

	// std::cout << "finding indices finished!" << std::endl;

	// std:: cout << "sizes: " << betaArr.at(0).at(0).size() << ", " << gammaArr.at(0).at(0).size() << ", " << deltaArr.at(0).at(0).size() << std::endl;
	
	// this->A_mtxArr.resize(this->nParts, std::vector<Mtx_CSR>(s));
	// this->L_mtxArr.resize(this->nParts, std::vector<Mtx_CSR>(s));
	// this->U_mtxArr.resize(this->nParts, std::vector<Mtx_CSR>(s));

	// #pragma omp parallel for	
	// for (size_t i = 0; i < this->nParts; ++i) {
		// for (size_t j = 0; j < s; ++j) {
			// extract(A_mtx, &A_mtxArr.at(i).at(j), gammaArr.at(i).at(j), deltaArr.at(i).at(j));
			// extract(A_mtx, &L_mtxArr.at(i).at(j), gammaArr.at(i).at(j), gammaArr.at(i).at(j));
			// extract(A_mtx, &U_mtxArr.at(i).at(j), betaArr.at(i).at(j), betaArr.at(i).at(j));
		// }
	// }
	
	// ksp->print(&U_mtxArr.at(1).at(0));		
	// std:: cout << "sizes: " << betaArr.at(0).at(0).size() << ", " << gammaArr.at(0).at(0).size() << ", " << deltaArr.at(0).at(0).size() << std::endl;

	mkl_free(perm);
	mkl_free(iperm);
	mkl_free(ia);
	mkl_free(ja);
	mkl_free(b_row_ptr);
	mkl_free(b_col_indx);

	return SPARSE_STATUS_SUCCESS;
}

sparse_status_t ILU0_ca::destroy_MtxArrs() {

	#pragma omp parallel for
	for (size_t i = 0; i < this->nParts; ++i) {
		for (size_t j = 0; j < ksp->getS(); ++j) {
			ksp->destroyMtx(&this->A_mtxArr.at(i).at(j));
			ksp->destroyMtx(&this->L_mtxArr.at(i).at(j));
			ksp->destroyMtx(&this->U_mtxArr.at(i).at(j));
		}
	}
	
	return SPARSE_STATUS_SUCCESS;
}

sparse_status_t ILU0_ca::find(Mtx_CSR *A_mtx,
															std::vector<size_t>& alpha, 
															std::vector<std::vector<size_t>>& betaArr,
															std::vector<std::vector<size_t>>& gammaArr,
															std::vector<std::vector<size_t>>& deltaArr)
{
	size_t s = ksp->getS();
	size_t n = A_mtx->n;
	deltaArr.at(0) = alpha;
	
	for (size_t i = 0; i < s; ++i) {
		neighborhood(A_mtx, // CSR MATRIX
								 deltaArr.at(i), // input set
								 betaArr.at(i), // result vertices for some level of order
								 n, // the reachabilty order
								 GRAPH_UPPER
		);
		
		neighborhood(A_mtx, // CSR MATRIX
								 betaArr.at(i), // input set
								 gammaArr.at(i), // result vertices for some level of order
								 n, // the reachabilty order
								 GRAPH_LOWER
		);
		
		neighborhood(A_mtx, // CSR MATRIX
								 gammaArr.at(i), // input set
								 deltaArr.at(i), // result vertices for some level of order
								 1, // the reachabilty order
								 GRAPH_COMPLETE
		);
				
		if (i < s - 1)
			deltaArr.at(i+1) = deltaArr.at(i);
	}
	return SPARSE_STATUS_SUCCESS;	
}
