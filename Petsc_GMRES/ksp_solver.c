#include "ksp_solver.h"
#define RUNS 1

PetscErrorCode compute(PetscReal *avg_rtime, PetscInt *its, PetscReal **a, PetscInt *na, PetscReal *rnorm) {
  Vec            x, b;		 /* approx solution, RHS */
  Mat            A;        /* linear system matrix */
  KSP            ksp;      /* linear solver context */
  PetscRandom    rctx;     /* random number generator context */
  PetscReal      *a_temp;	 /* residual histories */
  PetscInt       i, n;
  PetscErrorCode ierr;
  PetscBool      flg = PETSC_FALSE;
  PetscScalar    A_min, A_max;
	char           matrixInputFile[PETSC_MAX_PATH_LEN];
	MatrixInfo		 minfo;

	float rtime, ptime, mflops;
	long long flpops;
	float irtime, iptime, imflops;
  long long iflpops;	
	
	ierr = PetscOptionsGetString(NULL, PETSC_NULL, "-mmt", matrixInputFile, PETSC_MAX_PATH_LEN-1, &flg);CHKERRQ(ierr);
		
	/* Check if files where provided, raise an error if not */
	if(!flg){
		SETERRQ(PETSC_COMM_WORLD,1,"No matrix file neither vector one where provided, what am I supposed to convert ? Void ???\n");
	}
	
	MMTgetMatrixReal(matrixInputFile, &A, &minfo, &A_min, &A_max);
	n = minfo.n;

  /* A is symmetric. Set symmetric flag to enable ICC/Cholesky preconditioner */
	flg = PETSC_FALSE;
	ierr = PetscOptionsGetBool(NULL,NULL,"-sym",&flg,NULL);CHKERRQ(ierr);
  if(flg) {ierr = MatSetOption(A, MAT_SYMMETRIC, PETSC_TRUE);CHKERRQ(ierr);}
  
  ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
	ierr = VecSetSizes(b, PETSC_DECIDE, n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(b);CHKERRQ(ierr);
  ierr = VecDuplicate(b,&x);CHKERRQ(ierr);

  /*
     By default we use a RHS with all elements of 1.0;
		 Alternatively, using the runtime option
     '-random_b' forms a RHS with random components
		 with range [A_min, A_max].
  */
	flg = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL, NULL,"-random_b", &flg, NULL);CHKERRQ(ierr);
  if (flg) {		
    ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);CHKERRQ(ierr);
		ierr = PetscRandomSetSeed(rctx, time(NULL));
		ierr = PetscRandomSeed(rctx);CHKERRQ(ierr);
    ierr = PetscRandomSetInterval(rctx, A_min, A_max);CHKERRQ(ierr);     
		ierr = PetscRandomSetFromOptions(rctx);CHKERRQ(ierr);
		ierr = VecSetRandom(b,rctx);CHKERRQ(ierr);
    ierr = PetscRandomDestroy(&rctx);CHKERRQ(ierr);
  } else {
    ierr = VecSet(b, 1.0);CHKERRQ(ierr);
  }
		ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
		ierr = VecAssemblyEnd(b);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the linear solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

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
     Set linear solver defaults for this problem (optional).
     - By extracting the KSP and PC contexts from the KSP context,
       we can then directly call any KSP and PC routines to set
       various options.
     - The following two statements are optional; all of these
       parameters could alternatively be specified at runtime via
       KSPSetFromOptions().  All of these defaults can be
       overridden at runtime, as indicated below.
  */

  ierr = KSPSetTolerances(ksp, 1.e-2/((n+1)*(n+1)), 1.e-50, PETSC_DEFAULT, PETSC_DEFAULT);CHKERRQ(ierr);

  /*
    Set runtime options, e.g.,
        -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines.
  */
	
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
	
	/*
		 Set preconditioner and KSP type manually
	*/
	/*
	PC pc;
	KSPGetPC(ksp, &pc);
	PCSetType(pc, "cholesky");
	KSPSetPC(ksp, pc);
	KSPSetType(ksp, "preonly");
	*/

	ierr = KSPSetResidualHistory(ksp, PETSC_NULL, PETSC_DECIDE, PETSC_TRUE);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	*avg_rtime = 0;
	for (i = 0; i < RUNS; ++i) {	
		if (PAPI_flops(&irtime, &iptime, &iflpops, &imflops) != PAPI_OK)
			exit(1);
		
		ierr = KSPSolve(ksp, b, x);CHKERRQ(ierr);
		
		if (PAPI_flops(&rtime, &ptime, &flpops, &mflops) != PAPI_OK)
			exit(1);
		
		PAPI_shutdown();
		*avg_rtime += (PetscReal)rtime / RUNS;
	}

	/*
		Extract information
	*/
	ierr = KSPGetResidualHistory(ksp, &a_temp, na);
	ierr = KSPGetIterationNumber(ksp, its);CHKERRQ(ierr);
	ierr = KSPGetResidualNorm(ksp, rnorm);CHKERRQ(ierr);
	
	/*
		Allocate memory for 'a' and copy from 'a_temp' so that values 
		survive after call to KSPDestroy().
	*/
	PetscMalloc(*na * sizeof(a[0]), a);
	
	for (i = 0; i < *na; ++i) {
		(*a)[i] = a_temp[i];
	}

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	
  /*
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
  */
	
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);
	ierr = MatDestroy(&A);CHKERRQ(ierr);

	return(ierr);
}
