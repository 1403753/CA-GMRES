/*
  Include "petscksp.h" so that we can use KSP solvers.  Note that this file
  automatically includes:
     petscsys.h    - base PETSc routines   petscvec.h - vectors
     petscmat.h 	 - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners
*/
#ifndef _MMT_READER_H
#include "petscksp.h"
#define _MMT_READER_H


typedef struct _MatrixInfo{
	int n;
	int m;
	int nnz;
} MatrixInfo;

PetscErrorCode MMTgetMatrixReal(char *fin, Mat *A, MatrixInfo *minfo, PetscReal *A_min, PetscReal *A_max);

#endif //_MMT_READER_H