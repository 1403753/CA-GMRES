

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

#define MAIN_HPP

int main() {
	
	std::cout.setf(std::ios_base::fixed);
	std::cout.precision(5);

	double                 		   	rTol = 1e-11;                   // the relative (possibly preconditioned) residual norm || A*x_{k + 1} - b || / || A*x_0 - b || == || r_{k+1} || / || r_0 ||
	double                    		aTol = 1e-50;                   // the absolute (possibly preconditioned) residual norm || A*x_{k + 1} - b || == || r_{k+1} ||
	double                        dTol = 1e+4;                    // the divergence tolerance, amount (possibly preconditioned) residual norm can increase
	size_t                        maxit = 100;                    // maximum number of iterations to use
	sparse_status_t               stat;
	// sparse_matrix_t               A_mkl;                          // n x n matrix A
	// sparse_matrix_t               M;                              // n x n preconditioned matrix M == ilu0(A)
	Mtx_CSR                       A_mtx;                         	
	KSP						 							  ksp;							 						  // linear solver context	
	GMRES_ca											gmres;													// KSPType
	ILU0_ca												ilu0;														// PCType
	// const size_t                  n = 4;                          // dim(A)
	const size_t                  t = 12;                         // number of 'outer iterations' before restart
	const size_t									s = 5;													// step-size ('inner iterations')
	
	// read matrix
	stat = MmtReader::read_matrix_from_file("../matrix_market/sparse9x9complex.mtx", &A_mtx);
	if (stat);
	// size_t row_ptr[n] = {0,1,2,3};
	// size_t col_indx[n] = {0,1,2,3};
	// double values[n] = {1,1,1,1};
	
	// for (size_t i = 0; i < n + 1; ++i)
		// std::cout << A_mtx.row_ptr[i] << ", ";
	// std::cout << std::endl;
	
	// Mtx_CSR A_mtx{4, 4, row_ptr, col_indx, values};
	
	ksp.setOperators(&A_mtx, &A_mtx);
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
	
	return 0;
}

#endif