//============================================================================
// Name        : ca-gmres.cpp
// Author      : Robert Ernstbrunner
// Version     :
// Copyright   : Your copyright notice
// Description : ca-gmres in C++
//============================================================================

// #include "gmres.hpp"
#ifndef GMRES_HPP_

#include "arnoldi_ca.hpp"
#include "spmv.hpp"
#include "tsqr.hpp"
#include "matrix_reader.hpp"
#include <papi.h>

#include "_hypre_utilities.h"
#include "HYPRE_struct_ls.h"

#define GMRES_HPP_

#define s 3

int main() {
	
	std::cout << "int size: " << CHAR_BIT*sizeof(MKL_INT) << std::endl;
	
	// SimpleTest();
  // HYPRE_StructSolver   solver;

	
	std::vector<double, mkl_allocator<double>> test;
	test.reserve(2);
		
	std::cout.setf(std::ios_base::fixed);
	std::cout.precision(3);
	
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(rng, time(NULL));
	
	float							rtime, ptime, mflops;
	long long 				flpops;
	double 						*V, **V_col_ptr, *r, beta;
	sparse_status_t 	stat;	
	sparse_matrix_t 	A;
	MatrixInfo				minfo;
	size_t						indx;

	
	if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
		exit(1);
	
	// stat = matrix_reader<double>::read_matrix_from_file("../matrix_market/matlab_example.mtx", &A, &minfo);
	stat = matrix_reader<double>::read_matrix_from_file("../matrix_market/mini_test.mtx", &A, &minfo);
	// stat = matrix_reader<double>::read_matrix_from_file("../matrix_market/goodwin.mtx", &A, &minfo);
	// stat = matrix_reader<double>::read_matrix_from_file("../matrix_market/nasa4704.mtx", &A, &minfo);
	// stat = matrix_reader<double>::read_matrix_from_file("../matrix_market/bmw7st_1.mtx", &A, &minfo);
	
	if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
			exit(1);
	
	PAPI_shutdown();		
	
	printf("runtime readfile: %f\n", rtime);

	if(stat != SPARSE_STATUS_SUCCESS) throw std::invalid_argument("MatCreate failed");

	const size_t m = minfo.rows;

	r = (double *)mkl_calloc(m, sizeof(double), 64);if(r == NULL){return 1;}
	V_col_ptr = (double **)mkl_malloc(m * sizeof(double *), 64);if(V_col_ptr == NULL){return 1;}
	V = (double *)mkl_calloc(m * (s+1), sizeof(double), 64);if(V == NULL){return 1;}

	for(size_t i = 0; i < m; ++i)
		r[i] = 1;
	
	indx = 0;
	for(size_t i = 0; i < m * (s+1); i += m)
		V_col_ptr[indx++] = &V[i];

	
/*
	compute L2-norm beta
*/
// double cblas_dnrm2 (const MKL_INT n, const double *x, const MKL_INT incx);
	beta = cblas_dnrm2 (m, r, 1);
/*
	copy and scale residual vector into V matrix
*/
	cblas_dcopy (m, r, 1, V_col_ptr[0], 1);
// void	cblas_dscal (const int N, const double alpha, double *X, const int incX)
	cblas_dscal (m, 1 / beta, V_col_ptr[0], 1);
	
	if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
			exit(1);
	
	for(size_t i = 0; i < s; ++i) {

		spmv::mv (A, V_col_ptr[i], &V_col_ptr[i+1], s);
		
		beta = cblas_dnrm2 (m, V_col_ptr[i + 1], 1);
		cblas_dscal (m, 1 / beta, V_col_ptr[i + 1], 1);
		
		// printf("\n============= V:\n");
		// for(size_t j = 0; j < m; ++j) {
			// for(size_t k = 0; k < s+1; ++k) {
				// std::cout << V[k*m + j] << '\t';
			// }
			// std::cout << std::endl;
		// }
	}
	
	if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
			exit(1);
	
	PAPI_shutdown();	

	// printf("\n============= V col major:\n");
	// for(size_t j = 0; j < m*(s+1); ++j) {
		// std::cout << V[j] << ", ";
	// }
	// std::cout << std::endl;

	
	/*
		transpose matrix in place
	*/
	// void mkl_dimatcopy (const char ordering, const char trans, size_t rows, size_t cols, const double alpha, double * AB, size_t lda, size_t ldb);
	// mkl_dimatcopy('R', 'T', s+1, m, 1, V, m, s+1);

	// printf("\n============= V row major:\n");
	// for(size_t j = 0; j < m*(s+1); ++j) {
		// std::cout << V[j] << ", ";
	// }
	// std::cout << std::endl;

	// if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
			// exit(1);

	tsqr::qr(&V, m, s+1);
	
	// if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) < PAPI_OK)
		// exit(1);
	
	// PAPI_shutdown();	
	
	// arnoldi_ca::givens_rotations();
	// arnoldi_ca::modified_leja_ordering();

	stat = mkl_sparse_destroy(A);
	gsl_rng_free(rng);
	mkl_free(V);
	mkl_free(r);
	mkl_free(V_col_ptr);
	printf("runtime SpMV: %f\n", rtime);
	
	return 0;
}

/* in place matrix / vector creation

	size_t rows_ptr[] = {0, 1, 2, 3, 4, 5};
	size_t col_indx[] = {0, 0, 1, 0, 0};
	double values[] = {3, 2, 4, 9, 10};	
	double beta;
	double r[M] = {3,2,4,9,10};
	// double r[M] = {1,1,1,1,1};
  // 5 x 3 Matrix testcode:	
		// MKL_INT rows_ptr[] = {0, 3, 5, 8, 11, 13};
		// MKL_INT col_indx[] = {0, 1, 3, 0, 1, 2, 3, 4, 0, 2, 3, 1, 4};
		// double values[] = 	 {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	
	stat = mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO, M, N, rows_ptr, rows_ptr + 1, col_indx, values);
	if(stat != SPARSE_STATUS_SUCCESS)
		throw std::invalid_argument("MatCreate failed");
	
	V_col_ptr = (double **)mkl_malloc(M * sizeof(double *), 64);if(V_col_ptr == NULL){return 1;}
	V = (double *)mkl_calloc(M*N, sizeof(double), 64);if(V == NULL){return 1;}
	
	size_t idx = 0;
	for(size_t i = 0; i < M*N; i += M) {
		V_col_ptr[idx++] = &V[i];
	}
*/

#endif

/*********************/
/* ARPACK NOT NEEDED */
/*********************/

// arpack includes 
// #include "blas1c.h"
// #include "arssym.h"

// template<class T>
// class SymMatrix {
// /*

  // This simple class exemplifies how to create a matrix class 
  // that can be used by ARPACK++.
  // Basically, SymMatrix is required to have a member function 
  // that calculates the matrix-vector product SymMatrix*v, 
  // where v is a vector with elements of type T.

// */

 // private:

  // int  n;

 // public:

  // int  ncols() { return n; }

  // void MultMv(T* v, T* w)

  // /*
    // Matrix vector subroutine
    // where the matrix is the one dimensional discrete Laplacian on
    // the interval [0,1] with zero Dirichlet boundary conditions. 
  // */

  // {
  
    // int  j;
    // T    h2;

    // const T two = 2.0;

    // w[0] =  two*v[0] - v[1];
    // for (j=1; j<n-1; j++) {
      // w[j] = - v[j-1] + two*v[j] - v[j+1];
    // }
    // w[n-1] = - v[n-2] + two*v[n-1];

    // Scaling the vector w by (1 / h^2) using fortran blas routine scal.

    // h2 = T((n+1)*(n+1));
    // scal(n, h2, w, 1L);

  // } // MultMv

  // SymMatrix(int  nval) { n = nval; }

// }; // SymMatrix.


// template<class DOUBLE, class EIGPROB>
// void Solution(SymMatrix<DOUBLE> &A, EIGPROB &Prob)
// /*
  // This function prints eigenvalues and eigenvetors on standard "std::cout" 
  // stream and exemplifies how to retrieve information from ARPACK++ classes.
// */

// {

  // int   i, n, nconv, mode;
  // double *Ax;
	// double *ResNorm;

  // /*
     // ARPACK++ includes some functions that provide information
     // about the problem. For example, GetN furnishes the dimension
     // of the problem and ConvergedEigenvalues the number of 
     // eigenvalues that attained the required accuracy. GetMode 
     // indicates if the problem was solved in regular, 
     // shift-and-invert or other mode.
  // */

  // n     = Prob.GetN();
  // nconv = Prob.ConvergedEigenvalues();
  // mode  = Prob.GetMode();

  // std::cout << std::endl << std::endl << "Testing ARPACK++ class ARSymEig \n";
  // std::cout << "Real symmetric eigenvalue problem: A*x - lambda*x" << std::endl;
  // switch (mode) {
  // case 1:
    // std::cout << "Regular mode" << std::endl << std::endl;
    // break;
  // case 3: 
    // std::cout << "Shift and invert mode" << std::endl << std::endl;
  // }

  // std::cout << "Dimension of the system            : " << n             << std::endl;
  // std::cout << "Number of 'requested' eigenvalues  : " << Prob.GetNev() << std::endl;
  // std::cout << "Number of 'converged' eigenvalues  : " << nconv         << std::endl;
  // std::cout << "Number of Arnoldi vectors generated: " << Prob.GetNcv() << std::endl;
  // std::cout << std::endl;

  // /*
    // EigenvaluesFound is a boolean function that indicates
    // if the eigenvalues were found or not. Eigenvalue can be
    // used to obtain one of the "converged" eigenvalues. There
    // are other functions that return eigenvectors elements,
    // Schur vectors elements, residual vector elements, etc.
  // */

  // if (Prob.EigenvaluesFound()) {
    // std::cout << "Eigenvalues:" << std::endl;
    // for (i=0; i<nconv; i++) {
      // std::cout << "  lambda[" << (i+1) << "]: " << Prob.Eigenvalue(i) << std::endl;
    // }
    // std::cout << std::endl;
  // }

  // /*
    // EigenvectorsFound indicates if the eigenvectors are
    // available. RawEigenvector is one of the functions that
    // provide raw access to ARPACK++ output data. Other functions
    // of this type include RawEigenvalues, RawEigenvectors, 
    // RawSchurVector, RawResidualVector, etc.
  // */

  // if (Prob.EigenvectorsFound()) {

    // Printing the residual norm || A*x - lambda*x ||
    // for the nconv accurately computed eigenvectors.

    // Ax      = (double *)mkl_malloc(n * sizeof(double), 64);if(Ax == NULL){}
    // ResNorm = (double *)mkl_malloc(nconv * sizeof(double), 64);if(ResNorm == NULL){}

    // for (i=0; i<nconv; i++) {
      // A.MultMv(Prob.RawEigenvector(i),Ax);
      // axpy(n, -Prob.Eigenvalue(i), Prob.RawEigenvector(i), 1, Ax, 1);
      // ResNorm[i] = nrm2(n, Ax, 1) / fabs(Prob.Eigenvalue(i));
    // }

    // for (i=0; i<nconv; i++) {
      // std::cout << "||A*x(" << (i+1) << ") - lambda(" << (i+1);
      // std::cout << ")*x(" << (i+1) << ")||: " << ResNorm[i] << "\n";
    // }
    // std::cout << "\n";

    // mkl_free(Ax);
    // mkl_free(ResNorm);

  // }

// } // Solution


// void SimpleTest()
// {

  // int  nconv;

  // Creating a double precision 100x100 matrix.

  // SymMatrix<double> A(100);

  // Defining what we need: the four eigenvectors of A with smallest magnitude.

  // ARSymStdEig<double, SymMatrix<double> >
    // dprob(A.ncols(), 4, &A, &SymMatrix<double>::MultMv, "SM");

  // /*
    // It is possible to pass other parameters directly to the constructor
    // of class ARSymStdEig in order to define a problem. The list
    // of parameters includes, among other values, the maximum number of 
    // iterations allowed and the the relative accuracy used to define the
    // stopping criterion. Alternatively, it is also possible to use function 
    // DefineParameters to set ARPACK++ variables after declaring dprob as an 
    // object of ARSymStdEig class using the default constructor.
  // */

  // Finding eigenvectors.

  // nconv = dprob.FindEigenvectors();

  // /*
    // FindEigenvalues, FindArnoldiBasis and FindSchurVectors
    // are alternatives for FindEigenvectors. However, FindSchurVectors
    // can only be used with non-symmetric real and complex problems.
    // Some other functions like Eigenvectors, Eigenvalues and
    // EigenValVectors could also be used. These functions store the
    // desired data in user supplied vectors and matrices.
  // */

  // Printing solution.

  // Solution(A, dprob);

// } // SimpleTest
