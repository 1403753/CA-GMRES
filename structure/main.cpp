

/**************/
/*  END TODO  */
/**************/

/*********/
/* notes */
/*********/


#ifndef MAIN_HPP

#include "ca_iluz_gmres.hpp"
// #include "matrix_reader.hpp"
#include <omp.h>

#include <metis.h>

#define MAIN_HPP


int main() {
	
	std::cout.setf(std::ios_base::fixed);
	std::cout.precision(5);
	

	double                 		   	rTol;                           // the relative (possibly preconditioned) residual norm || A*x_{k + 1} - b || / || A*x_0 - b || == || r_{k+1} || / || r_0 ||
	double                    		aTol;                           // the absolute (possibly preconditioned) residual norm || A*x_{k + 1} - b || == || r_{k+1} ||
	double                        dTol;                           // the divergence tolerance, amount (possibly preconditioned) residual norm can increase
	size_t                        maxit;                          // maximum number of iterations to use
	sparse_status_t               stat;
	// sparse_matrix_t               A_mkl;                          // n x n matrix A
	// sparse_matrix_t               M;                              // n x n preconditioned matrix M == ilu0(A)
	// Mtx_CSR      								 A_mtx;                         	
	// ca_iluz_gmres 								ksp;														// linear solver context	
	size_t                        n;                              // dim(A)
	const size_t                  t = 12;                         // number of 'outer iterations' before restart
	const size_t									s = 5;													// step-size ('inner iterations')
	
	
	
	// read matrix
	// stat = matrix_reader<double>::read_matrix_from_file("../matrix_market/sparse9x9complex.mtx", &A_mtx);
	
	size_t row_ptr[4] = {0,1,2,3};
	size_t col_indx[4] = {0,1,2,3};
	double values[4] = {1,1,1,1};
	
	Mtx_CSR A_mtx{4, 4, row_ptr, col_indx, values};
	
	ca_iluz_gmres ksp{A_mtx};
	
	std::cout << ksp.pub << " hehe\n";
	
	// ksp.SetOptions(s,      // step-size
								 // t,      // restart length
								 // rTol,   
								 // aTol,  
								 // dTol,   
								 // maxit	  
							  // );
	
	// ksp.SetUp();
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