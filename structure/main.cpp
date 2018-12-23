
////////////
//  TODO  //
////////////
/*
	add non preconditioner-class
	add residual history
	add command promt parameter reader
*/

#ifndef MAIN_HPP

#include "KSP.hpp"
#include "GMRES_ca.hpp"
#include "PCILU0_ca.hpp"
#include "PCNone.hpp"
#include "MmtReader.hpp"
#include <papi.h>

#define MAIN_HPP

int main() {
	
	float                         rtime, ptime, mflops;
	long long                     flpops;	

	double                 		   	rTol = 1e-11;                   // the relative (possibly preconditioned) residual norm || A*x_{k + 1} - b || / || A*x_0 - b ||
                                                                // == || r_{k+1} || / || r_0 ||
	double                    		aTol = 1e-50;                   // the absolute (possibly preconditioned) residual norm || A*x_{k + 1} - b || == || r_{k+1} ||
	double                        dTol = 1e+4;                    // the divergence tolerance, amount (possibly preconditioned) residual norm can increase
	size_t                        maxit = 100000;                   // maximum number of iterations to use
	sparse_status_t               stat;
	sparse_matrix_t               A_mkl;                          // n x n matrix A
	double                        *b;
	double                        *x;
	double                        *tx;
	KSP						 							  ksp;							 						  // linear solver context	
	PCILU0_ca                     ilu0;                           // PCType
	PCNone                        pcnone;                         // PCType
	MmtReader											mmtReader;
	const size_t                  t = 12;                         // number of 'outer iterations' before restart
	const size_t									s = 5;													// step-size ('inner iterations')
	// Basis                         basis = MONOMIAL;
	Basis                         basis = NEWTON;
	GMRES_ca											gmres(s, t, basis);             // KSPType
	
	sparse_index_base_t           indexing;	
	size_t                        n, m;
	size_t                        *rows_start;
	size_t                        *rows_end;
	size_t                        *col_indx;
	double                        *values;
	struct matrix_descr           descr;
	mkl_set_num_threads(2);
	// read matrix
	// stat = mmtReader.read_matrix_from_file("../matrix_market/e05r0000.mtx", &A_mkl);
	stat = mmtReader.read_matrix_from_file("../matrix_market/watt1.mtx", &A_mkl);
	// stat = mmtReader.read_matrix_from_file("../matrix_market/sparse9x9complex.mtx", &A_mkl);
	// stat = mmtReader.read_matrix_from_file("../matrix_market/bcsstk18.mtx", &A_mkl);
	// stat = mmtReader.read_matrix_from_file("../matrix_market/bcsstk08.mtx", &A_mkl);
	// stat = mmtReader.read_matrix_from_file("../matrix_market/xenon2.mtx", &A_mkl);
	// stat = mmtReader.read_matrix_from_file("../matrix_market/bmw7st_1.mtx", &A_mkl);
	// stat = mmtReader.read_matrix_from_file("../matrix_market/goodwin.mtx", &A_mkl);
	// stat = mmtReader.read_matrix_from_file("../matrix_market/dwb512.mtx", &A_mkl);
	// stat = mmtReader.read_matrix_from_file("../matrix_market/1138_bus.mtx", &A_mkl);
	// stat = mmtReader.read_matrix_from_file("../matrix_market/nasa4704.mtx", &A_mkl);
	// stat = mmtReader.read_matrix_from_file("../matrix_market/mini_test.mtx", &A_mkl); // too small, to work properly
	// stat = mmtReader.read_matrix_from_file("../matrix_market/CA-ILU(0).mtx", &A_mkl); // too small, to work properly
	
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;

	stat = mkl_sparse_d_export_csr(A_mkl, &indexing, &n, &m, &rows_start, &rows_end, &col_indx, &values); // n is needed
	
	b = (double *) mkl_malloc(n * sizeof(double), 64);if(b == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	x = (double *) mkl_calloc(n, sizeof(double), 64);if(x == NULL){return SPARSE_STATUS_ALLOC_FAILED;}
	tx = (double *) mkl_malloc(n * sizeof(double), 64);if(tx == NULL){return SPARSE_STATUS_ALLOC_FAILED;}

	
	for (size_t i = 0; i < n; ++i)
		tx[i] = 0.1*(i%3 + 1);
		// tx[i] = 1;

	// for (size_t i = 0; i < n; ++i)
		// x[i] = 0.1*(i%3 + 1);

	mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A_mkl, descr, tx, 0, b);	

	stat = ksp.setOperator(&A_mkl);
	stat = ksp.setKSPType(&gmres);
	// stat = ksp.setPCType(&ilu0);
	stat = ksp.setPCType(&pcnone);

	stat = ksp.setOptions(rTol, aTol, dTol, maxit);

	stat = ksp.setUp();
	
	if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
		exit(1);	

	ksp.solve(x, b);

	if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
		exit(1);

	PAPI_shutdown();

	std::cout << "runtime ksp-solve: " << std::scientific << rtime << " seconds." << std::endl;

	std::cout << "solution x:" << std::endl;
	for (size_t i = 0; i < (n < 10 ? n : 10); ++i)
		std::cout << std::scientific << x[i] << "\n";
	if (n > 10)
		std::cout << ".\n.\n.\n";
	for (size_t i = n - 10; i < n; ++i) // size_t can't be negative
		std::cout << std::scientific << x[i] << "\n";

	stat = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A_mkl, descr, x, 0, tx);	

	for (size_t i = 0; i < n; ++i) {
		x[i] = b[i] - tx[i];
	}
	double r_knrm = cblas_dnrm2(n, x, 1);

	double r_0nrm = cblas_dnrm2(n, b, 1);

	double rRes = r_knrm / r_0nrm;	
	std::cout << "\n============= final main rel. res.: ";
	printf("%e\n",rRes);
	std::cout << "r_knrm: " << r_knrm << ", r_0nrm: " << r_0nrm << std::endl;

	stat = mkl_sparse_destroy(A_mkl);
	mkl_free(x);
	mkl_free(tx);
	mkl_free(b);
	mkl_free_buffers();
	return stat;
}

#endif