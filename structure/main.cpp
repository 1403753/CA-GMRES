
////////////
//  TODO  //
////////////
/*
	add non preconditioner-class
	add residual history
	add command promt parameter reader
*/

#ifndef MAIN_HPP

#include "KSP_ca.hpp"
#include "GMRES_ca.hpp"
#include "GMRES.hpp"
#include "PCILU0_ca.hpp"
#include "PCNone.hpp"
#include "MmtReader.hpp"
#include <papi.h>

#include <petscksp.h>

static char help[] = "Reads a PETSc matrix and vector from a file and solves a linear system.\n";
PetscErrorCode petsc(int argc, char **args);

#define MAIN_HPP

int main(int argc, char **args) {
	
	float                         rtime, ptime, mflops;
	long long                     flpops;	

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
	KSP_ca                        ksp;							 						  // linear solver context	
	PCILU0_ca                     ilu0;                           // PCType
	PCNone                        pcnone;                         // PCType
	MmtReader											mmtReader;
	const size_t                  t = 20;                         // number of 'outer iterations' before restart
	const size_t									s = 5;													// step-size ('inner iterations')
	// Basis                         basis = MONOMIAL;
	Basis                         basis = NEWTON;
	GMRES_ca											gmres(s, t, basis);             // KSPType
	// GMRES                         gmres(100);                     // KSPType
	
	sparse_index_base_t           indexing;	
	size_t                        n, m;
	size_t                        *rows_start;
	size_t                        *rows_end;
	size_t                        *col_indx;
	double                        *values;
	struct matrix_descr           descr;
	mkl_set_num_threads(48);
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
	stat = ksp.setPCType(&ilu0);
	// stat = ksp.setPCType(&pcnone);

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
	

	
	/////////////
	//  PETSc  //
	/////////////
	
	// petsc(argc, args);

	
	
	stat = mkl_sparse_destroy(A_mkl);
	mkl_free(x);
	mkl_free(tx);
	mkl_free(b);
	mkl_free_buffers();

	return stat;
}

PetscErrorCode petsc(int argc, char **args) {
  Vec            x,b,u;      /* approx solution, RHS, exact solution */
  Mat            A;            /* linear system matrix */
  KSP            ksp;         /* KSP context */
  KSP            *subksp;     /* array of local KSP contexts on this processor */
  PC             pc;           /* PC context */
  PC             subpc;        /* PC context for subdomain */
  PetscReal      norm;         /* norm of solution error */
  PetscErrorCode ierr;
  PetscInt       i,j,Ii,J,*blks,m = 4,n;
  PetscMPIInt    rank,size;
  PetscInt       its,nlocal,first,Istart,Iend;
  PetscScalar    v,one = 1.0,none = -1.0;
  PetscBool      isbjacobi;

  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  n    = m+2;

  /* -------------------------------------------------------------------
         Compute the matrix and right-hand-side vector that define
         the linear system, Ax = b.
     ------------------------------------------------------------------- */

  /*
     Create and assemble parallel matrix
  */
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m*n,m*n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(A,5,NULL,5,NULL);CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(A,5,NULL);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  for (Ii=Istart; Ii<Iend; Ii++) {
    v = -1.0; i = Ii/n; j = Ii - i*n;
    if (i>0)   {J = Ii - n; ierr = MatSetValues(A,1,&Ii,1,&J,&v,ADD_VALUES);CHKERRQ(ierr);}
    if (i<m-1) {J = Ii + n; ierr = MatSetValues(A,1,&Ii,1,&J,&v,ADD_VALUES);CHKERRQ(ierr);}
    if (j>0)   {J = Ii - 1; ierr = MatSetValues(A,1,&Ii,1,&J,&v,ADD_VALUES);CHKERRQ(ierr);}
    if (j<n-1) {J = Ii + 1; ierr = MatSetValues(A,1,&Ii,1,&J,&v,ADD_VALUES);CHKERRQ(ierr);}
    v = 4.0; ierr = MatSetValues(A,1,&Ii,1,&Ii,&v,ADD_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /*
     Create parallel vectors
  */
  ierr = VecCreate(PETSC_COMM_WORLD,&u);CHKERRQ(ierr);
  ierr = VecSetSizes(u,PETSC_DECIDE,m*n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(u);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&b);CHKERRQ(ierr);
  ierr = VecDuplicate(b,&x);CHKERRQ(ierr);

  /*
     Set exact solution; then compute right-hand-side vector.
  */
  ierr = VecSet(u,one);CHKERRQ(ierr);
  ierr = MatMult(A,u,b);CHKERRQ(ierr);

  /*
     Create linear solver context
  */
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

  /*
     Set operators. Here the matrix that defines the linear system
     also serves as the preconditioning matrix.
  */
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);

  /*
     Set default preconditioner for this program to be block Jacobi.
     This choice can be overridden at runtime with the option
        -pc_type <type>
  */
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCBJACOBI);CHKERRQ(ierr);


  /* -------------------------------------------------------------------
                   Define the problem decomposition
     ------------------------------------------------------------------- */

  /*
     Call PCBJacobiSetTotalBlocks() to set individually the size of
     each block in the preconditioner.  This could also be done with
     the runtime option
         -pc_bjacobi_blocks <blocks>
     Also, see the command PCBJacobiSetLocalBlocks() to set the
     local blocks.

      Note: The default decomposition is 1 block per processor.
  */
  ierr = PetscMalloc1(m,&blks);CHKERRQ(ierr);
  for (i=0; i<m; i++) blks[i] = n;
  ierr = PCBJacobiSetTotalBlocks(pc,m,blks);CHKERRQ(ierr);
  ierr = PetscFree(blks);CHKERRQ(ierr);


  /* -------------------------------------------------------------------
               Set the linear solvers for the subblocks
     ------------------------------------------------------------------- */

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Basic method, should be sufficient for the needs of most users.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     By default, the block Jacobi method uses the same solver on each
     block of the problem.  To set the same solver options on all blocks,
     use the prefix -sub before the usual PC and KSP options, e.g.,
          -sub_pc_type <pc> -sub_ksp_type <ksp> -sub_ksp_rtol 1.e-4
  */

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Advanced method, setting different solvers for various blocks.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     Note that each block's KSP context is completely independent of
     the others, and the full range of uniprocessor KSP options is
     available for each block. The following section of code is intended
     to be a simple illustration of setting different linear solvers for
     the individual blocks.  These choices are obviously not recommended
     for solving this particular problem.
  */
  ierr = PetscObjectTypeCompare((PetscObject)pc,PCBJACOBI,&isbjacobi);CHKERRQ(ierr);
  if (isbjacobi) {
    /*
       Call KSPSetUp() to set the block Jacobi data structures (including
       creation of an internal KSP context for each block).

       Note: KSPSetUp() MUST be called before PCBJacobiGetSubKSP().
    */
    ierr = KSPSetUp(ksp);CHKERRQ(ierr);

    /*
       Extract the array of KSP contexts for the local blocks
    */
    ierr = PCBJacobiGetSubKSP(pc,&nlocal,&first,&subksp);CHKERRQ(ierr);

    /*
       Loop over the local blocks, setting various KSP options
       for each block.
    */
    for (i=0; i<nlocal; i++) {
      ierr = KSPGetPC(subksp[i],&subpc);CHKERRQ(ierr);
      if (!rank) {
        if (i%2) {
          ierr = PCSetType(subpc,PCILU);CHKERRQ(ierr);
        } else {
          ierr = PCSetType(subpc,PCNONE);CHKERRQ(ierr);
          ierr = KSPSetType(subksp[i],KSPBCGS);CHKERRQ(ierr);
          ierr = KSPSetTolerances(subksp[i],1.e-6,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
        }
      } else {
        ierr = PCSetType(subpc,PCJACOBI);CHKERRQ(ierr);
        ierr = KSPSetType(subksp[i],KSPGMRES);CHKERRQ(ierr);
        ierr = KSPSetTolerances(subksp[i],1.e-6,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
      }
    }
  }

  /* -------------------------------------------------------------------
                      Solve the linear system
     ------------------------------------------------------------------- */

  /*
    Set runtime options
  */
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

  /*
     Solve the linear system
  */
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);

  /* -------------------------------------------------------------------
                      Check solution and clean up
     ------------------------------------------------------------------- */

  /*
     Check the error
  */
  ierr = VecAXPY(x,none,u);CHKERRQ(ierr);
  ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g iterations %D\n",(double)norm,its);CHKERRQ(ierr);

  /*
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
  */
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = VecDestroy(&u);CHKERRQ(ierr);  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}

#endif