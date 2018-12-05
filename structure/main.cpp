

/**************/
/*  END TODO  */
/**************/

/*********/
/* notes */
/*********/


#ifndef MAIN_HPP

#include "KSP.hpp"
#include "GMRES_ca.hpp"
#include "ILU0_ca.hpp"
#include "MmtReader.hpp"
#include <papi.h>

#define MAIN_HPP

int main() {
	
	float                         rtime, ptime, mflops;
	long long                     flpops;	
	
	std::cout.setf(std::ios_base::fixed);
	std::cout.precision(5);

	double                 		   	rTol = 1e-11;                   // the relative (possibly preconditioned) residual norm || A*x_{k + 1} - b || / || A*x_0 - b || == || r_{k+1} || / || r_0 ||
	double                    		aTol = 1e-50;                   // the absolute (possibly preconditioned) residual norm || A*x_{k + 1} - b || == || r_{k+1} ||
	double                        dTol = 1e+4;                    // the divergence tolerance, amount (possibly preconditioned) residual norm can increase
	size_t                        maxit = 100;                    // maximum number of iterations to use
	sparse_status_t               stat;
	sparse_matrix_t               A_mkl;                          // n x n matrix A
	// sparse_matrix_t               M_mkl;                          // n x n preconditioned matrix M == ilu0(A)
	// Mtx_CSR                       A_mtx;                         	
	KSP						 							  ksp;							 						  // linear solver context	
	GMRES_ca											gmres;													// KSPType
	ILU0_ca												ilu0;														// PCType
	MmtReader											mmtReader;
	// const size_t                  n = 4;                          // dim(A)
	const size_t                  t = 12;                         // number of 'outer iterations' before restart
	const size_t									s = 5;													// step-size ('inner iterations')
	
	if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
		exit(1);	
	
	// read matrix
	mmtReader.read_matrix_from_file("../matrix_market/sparse9x9complex.mtx", &A_mkl);
	// mmtReader.read_matrix_from_file("../matrix_market/bmw7st_1.mtx", &A_mkl);
	// mmtReader.read_matrix_from_file("../matrix_market/nasa4704.mtx", &A_mkl);
	// mmtReader.read_matrix_from_file("../matrix_market/mini_test.mtx", &A_mkl);

	if (stat);

	if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
		exit(1);

	PAPI_shutdown();

	std::cout << "runtime readfile: " << rtime << " seconds." << std::endl;
	
	ksp.setOperator(&A_mkl);
	ksp.setKSPType(&gmres);
	ksp.setPCType(&ilu0);
	
	ksp.setOptions(s, t, rTol, aTol, dTol, maxit);
	
	ksp.setUp();
	// what happens:
	// compute structure at_plus_a(A_mtx) ->output: b_rowptr and partition A with metis 
	// sort A_mtx with amml(s)
	// permute A_mtx -> col_indx out of order
	// create A_mkl and sort col_indx (maybe create own algorithm)
	// export A_mtx from A_mkl again
	// factorize A_mtx to get ILU(0)-values
	// find subsets alpha_p p == #threads depending on matrix size
	// openmp: find s-step dependencies for each alpha_p -> beta_p -> gamma_p -> delta_p 
	// 
	
	// ksp.solve(x, b);
	
	mkl_sparse_destroy(A_mkl);
	
	return 0;
}

#endif