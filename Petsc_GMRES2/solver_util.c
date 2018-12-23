#include <petscmat.h>
#include <petscksp.h>
#include <petscvec.h>
#include <petscpctypes.h>
#include <petscsys.h>

#include "solver_util.h"

#define MY_COMM PETSC_COMM_SELF

int load_matrix(Mat* A, char* file){
  PetscViewer     fd;                     /* viewer */
  PetscErrorCode  ierr;

  ierr = PetscViewerBinaryOpen(MY_COMM,file,FILE_MODE_READ,&fd);CHKERRQ(ierr);

  ierr = MatCreate(MY_COMM,A); CHKERRQ(ierr);
  ierr = MatSetType(*A,MATSEQAIJ); CHKERRQ(ierr);
  ierr = MatLoad(*A,fd); CHKERRQ(ierr);

  ierr = PetscViewerDestroy(&fd); CHKERRQ(ierr);

  return 0;
}

int linsys_solve(Mat* A, Vec* b, Vec* x, KSPType kspmethod, PCType pcmethod) {
  PetscErrorCode ierr;
  KSP ksp;
  PC pc;

  ierr = KSPCreate(MY_COMM, &ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, *A, *A); CHKERRQ(ierr);
  ierr = KSPSetType(ksp, kspmethod); CHKERRQ(ierr);
  ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
  ierr = PCSetType(pc, pcmethod); CHKERRQ(ierr);

  ierr = KSPSolve(ksp, *b, *x); CHKERRQ(ierr);

  ierr = KSPDestroy(&ksp); CHKERRQ(ierr);

  return 0;
}

PetscReal residual(Mat* A, Vec* b, Vec* x) {
  PetscErrorCode ierr;
  Vec y;
  PetscReal result;

  ierr = VecCreate(MY_COMM, &y); CHKERRQ(ierr);
  ierr = VecDuplicate(*b, &y); CHKERRQ(ierr);

  ierr = MatMult(*A, *x, y); CHKERRQ(ierr);
  ierr = VecAXPY(y, -1, *b); CHKERRQ(ierr);

  ierr = VecNorm(y, NORM_2, &result); CHKERRQ(ierr);

  ierr = VecDestroy(&y); CHKERRQ(ierr);

  return result;
}




int create_vec_random_from_mat(Mat* A, Vec* b, PetscRandom* r){
  PetscErrorCode ierr;
  PetscReal min, max;
  PetscInt rows, cols;

  //get dimensions of matrix
  ierr = MatGetSize(*A, &rows, &cols); CHKERRQ(ierr);

  //create b vector depending on matrix size
  ierr = VecCreate(MY_COMM, b); CHKERRQ(ierr);
  ierr = VecSetSizes(*b, PETSC_DECIDE, cols); CHKERRQ(ierr); //vector, local size, global size
  ierr = VecSetType(*b, VECSTANDARD); CHKERRQ(ierr);

  //create help vector to store the matrix rowmax / rowmin values
  Vec help;
  ierr = VecCreate(MY_COMM, &help); CHKERRQ(ierr);
  ierr = VecSetSizes(help, PETSC_DECIDE, rows); CHKERRQ(ierr);
  ierr = VecSetType(help, VECSTANDARD); CHKERRQ(ierr);

  //get minimum and maximum matrix entries -> will be stored in min and max
  ierr = MatGetRowMax(*A, help, PETSC_NULL); CHKERRQ(ierr); //matrix, vector to store max values, vector to store indices of max values
  ierr = VecMax(help, PETSC_NULL, &max); CHKERRQ(ierr); //vector, index of max value, max value
  ierr = MatGetRowMin(*A, help, PETSC_NULL); CHKERRQ(ierr);
  ierr = VecMin(help, PETSC_NULL, &min); CHKERRQ(ierr);

  //test output
  /*
  printf("Max of help vector: %f\n", max);
  printf("Min of help vector: %f\n", min);
  */


  //set random generator interval
  ierr = PetscRandomSetInterval(*r, min, max); CHKERRQ(ierr);

  //set vector entries to random values
  ierr = VecSetRandom(*b, *r); CHKERRQ(ierr);

  //free memory
  ierr = VecDestroy(&help); CHKERRQ(ierr);


  return 0;
}
