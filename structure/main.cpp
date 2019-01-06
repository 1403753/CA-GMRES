
////////////
//  TODO  //
////////////
/*
	clean up
	add init SpMV time measurement for CA-GMRES
*/

#ifndef MAIN_HPP

#include "KSP_.hpp"
#include "GMRES_ca.hpp"
#include "GMRES.hpp"
#include "PCILU0_ca.hpp"
#include "PCNone.hpp"
#include "MmtReader.hpp"
#include <papi.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>

#define MAIN_HPP

int main(int argc, char **args) {
	
	float                         gmres_ca_rtime;
	double                 		   	rTol = 1e-11;                   // the relative (possibly preconditioned) residual norm || A*x_{k + 1} - b || / || A*x_0 - b ||
                                                                // == || r_{k+1} || / || r_0 ||
	double                    		aTol = 1e-50;                   // the absolute (possibly preconditioned) residual norm || A*x_{k + 1} - b || == || r_{k+1} ||
	double                        dTol = 1e+4;                    // the divergence tolerance, amount (possibly preconditioned) residual norm can increase
	size_t                        maxit = 1000;                   // maximum number of iterations to use
	sparse_status_t               stat;
	sparse_matrix_t               A_mkl;                          // n x n matrix A
	double                        *b;
	double                        *x;
	double                        *tx;
	double                        *r;
	KSP_                          ksp;							 						  // linear solver context	
	PCILU0_ca                     ilu0;                           // PCType
	PCNone                        pcnone;                         // PCType
	MmtReader											mmtReader;
	const size_t									s = 5;													// step-size ('inner iterations')
	const size_t                  t = 12;                          // number of 'outer iterations' before restart
	// Basis                         basis = MONOMIAL;
	// Basis                         basis = NEWTON;
	GMRES_ca											gmres_ca(s, t, NEWTON);         // KSPType
	GMRES                         gmres(60);                      // KSPType
	
	sparse_index_base_t           indexing;	
	size_t                        n, m;
	size_t                        *rows_start;
	size_t                        *rows_end;
	size_t                        *col_indx;
	double                        *values;
	struct matrix_descr           descr;
	std::ofstream                 file;
	
	mkl_set_num_threads(8);
	// read matrix
	// stat = mmtReader.read_matrix_from_file("../matrix_market/dmat/dmat3.mtx", &A_mkl);
	// stat = mmtReader.read_matrix_from_file("../matrix_market/watt1.mtx", &A_mkl);
	// stat = mmtReader.read_matrix_from_file("../matrix_market/e05r0000.mtx", &A_mkl);
	// stat = mmtReader.read_matrix_from_file("../matrix_market/sparse9x9complex.mtx", &A_mkl);
	// stat = mmtReader.read_matrix_from_file("../matrix_market/bcsstk18.mtx", &A_mkl);
	// stat = mmtReader.read_matrix_from_file("../matrix_market/bcsstk08.mtx", &A_mkl);
	// stat = mmtReader.read_matrix_from_file("../matrix_market/xenon2.mtx", &A_mkl);
	stat = mmtReader.read_matrix_from_file("../matrix_market/pwtk.mtx", &A_mkl);
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
	r = (double *) mkl_malloc(n * sizeof(double), 64);if(r == NULL){return SPARSE_STATUS_ALLOC_FAILED;}

	gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(rng, time(NULL)); // Seed with time	

	for (size_t i = 0; i < n; ++i)
		tx[i] = gsl_ran_flat(rng, -1, 1) + std::sin(2*M_PI*i/n);

	mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1, A_mkl, descr, tx, 0, b);	

	stat = ksp.setOperator(&A_mkl);
	stat = ksp.setKSPType(&gmres_ca);
	stat = ksp.setPCType(&pcnone);

	stat = ksp.setOptions(rTol, aTol, dTol, maxit, true);

	stat = ksp.setUp();

	std::cout << " solve CA-GMRES" << std::endl;
	
	ksp.solve(x, b);

  // gmres_ca NEWTON

	stat = ksp.setKSPType(&gmres);
	
	std::fill(x, x + n, 0);

	std::cout << " solve std. GMRES" << std::endl;
	
	ksp.solve(x, b);

	gmres_ca_rtime = gmres_ca.getSDO() + gmres_ca.getBCGS() + gmres_ca.getTSQR() + gmres_ca.getMGS() + gmres_ca.getSpMV();
	
	file.open("../gnuplot/gmres.dat");
	file << 1 << " " << gmres.getMGS() / gmres_ca_rtime << " " << gmres.getSpMV() / gmres_ca_rtime << std::endl;
	file.close();

	file.open("../gnuplot/gmres_ca.dat");
	file << 1 << " " << gmres_ca.getSDO() / gmres_ca_rtime << " " << gmres_ca.getBCGS() / gmres_ca_rtime << " " << gmres_ca.getTSQR() / gmres_ca_rtime << " " << gmres_ca.getMGS() / gmres_ca_rtime << " " << gmres_ca.getSpMV() / gmres_ca_rtime << std::endl;
	file.close();

	stat = mkl_sparse_destroy(A_mkl);


	// next matrix

	stat = mmtReader.read_matrix_from_file("../matrix_market/bmw7st_1.mtx", &A_mkl);

	stat = ksp.setOperator(&A_mkl);
	stat = ksp.setKSPType(&gmres);
	stat = ksp.setPCType(&pcnone);

	stat = ksp.setOptions(rTol, aTol, dTol, maxit, true);

	stat = ksp.setUp();
	
	std::cout << " solve std. GMRES" << std::endl;

	ksp.solve(x, b);

  // gmres_ca NEWTON

	stat = ksp.setKSPType(&gmres_ca);
	
	std::fill(x, x + n, 0);

	std::cout << " solve CA-GMRES" << std::endl;
	
	ksp.solve(x, b);

	gmres_ca_rtime = gmres_ca.getSDO() + gmres_ca.getBCGS() + gmres_ca.getTSQR() + gmres_ca.getMGS() + gmres_ca.getSpMV();
	
	file.open("../gnuplot/gmres.dat", std::ios_base::app);
	file << 2 << " " << gmres.getMGS() / gmres_ca_rtime << " " << gmres.getSpMV() / gmres_ca_rtime << std::endl;
	file.close();

	file.open("../gnuplot/gmres_ca.dat", std::ios_base::app);
	file << 2 << " " << gmres_ca.getSDO() / gmres_ca_rtime << " " << gmres_ca.getBCGS() / gmres_ca_rtime << " " << gmres_ca.getTSQR() / gmres_ca_rtime << " " << gmres_ca.getMGS() / gmres_ca_rtime << " " << gmres_ca.getSpMV() / gmres_ca_rtime << std::endl;
	file.close();

	stat = mkl_sparse_destroy(A_mkl);
	

	// next matrix
	
	stat = mmtReader.read_matrix_from_file("../matrix_market/xenon2.mtx", &A_mkl);

	stat = ksp.setOperator(&A_mkl);
	stat = ksp.setKSPType(&gmres);
	stat = ksp.setPCType(&pcnone);

	stat = ksp.setOptions(rTol, aTol, dTol, maxit, true);

	stat = ksp.setUp();
	
	std::cout << " solve std. GMRES" << std::endl;

	ksp.solve(x, b);

  // gmres_ca NEWTON

	stat = ksp.setKSPType(&gmres_ca);
	
	std::fill(x, x + n, 0);

	std::cout << " solve CA-GMRES" << std::endl;
	
	ksp.solve(x, b);

	gmres_ca_rtime = gmres_ca.getSDO() + gmres_ca.getBCGS() + gmres_ca.getTSQR() + gmres_ca.getMGS() + gmres_ca.getSpMV();
	
	file.open("../gnuplot/gmres.dat", std::ios_base::app);
	file << 3 << " " << gmres.getMGS() / gmres_ca_rtime << " " << gmres.getSpMV() / gmres_ca_rtime << std::endl;
	file.close();

	file.open("../gnuplot/gmres_ca.dat", std::ios_base::app);
	file << 3 << " " << gmres_ca.getSDO() / gmres_ca_rtime << " " << gmres_ca.getBCGS() / gmres_ca_rtime << " " << gmres_ca.getTSQR() / gmres_ca_rtime << " " << gmres_ca.getMGS() / gmres_ca_rtime << " " << gmres_ca.getSpMV() / gmres_ca_rtime << std::endl;
	file.close();

	
	stat = mkl_sparse_destroy(A_mkl);
	
	gsl_rng_free(rng);
	
	mkl_free(x);
	mkl_free(tx);
	mkl_free(b);
	mkl_free(r);
	mkl_free_buffers();

	return stat;
}

#endif