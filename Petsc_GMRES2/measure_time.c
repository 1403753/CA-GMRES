#include <stdio.h>
#include <papi.h>
#include <petscksp.h>
#include <petscmat.h>

#include "measure_time.h"

#define MY_COMM PETSC_COMM_SELF

//struct performance_info measure_time(Mat* A, Vec* b, Vec* x, KSPType kspmethod, PCType pcmethod, struct performance_info *info) {
int measure_time(Mat* A, Vec* b, Vec* x, KSPType kspmethod, PCType pcmethod, struct performance_info *info, PetscReal* residualHistory, int* history_length) {
  PetscErrorCode ierr;

  //setting up the KSP
  KSP ksp;
  PC pc;
  ierr = KSPCreate(MY_COMM, &ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, *A, *A); CHKERRQ(ierr);
  ierr = KSPSetType(ksp, kspmethod); CHKERRQ(ierr);
  ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
  ierr = PCSetType(pc, pcmethod); CHKERRQ(ierr);

  ierr = KSPSetResidualHistory(ksp, residualHistory, *history_length, PETSC_TRUE); CHKERRQ(ierr);

  //PAPI performance measurement
  int retval = 0;
  float real_time = 0;
  float proc_time = 0;
  long long flpins = 0;
  float mflops = 0;

  retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops);
  if (retval != PAPI_OK) {
    printf("PAPI error");
    exit(1);
  }

  //function call that will be measured
  ierr = KSPSolve(ksp, *b, *x); CHKERRQ(ierr);

  //PAPI performance measurement
  retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops);
  if (retval != PAPI_OK) {
    printf("PAPI error");
    exit(1);
  }

  //struct performance_info info;
  info->real_time = real_time;
  info->mflops = mflops;

  PAPI_shutdown();

  ierr = KSPGetResidualHistory(ksp, &residualHistory, history_length);

  //free KSP
  ierr = KSPDestroy(&ksp); CHKERRQ(ierr);

  //return info;
  return 0;
}
