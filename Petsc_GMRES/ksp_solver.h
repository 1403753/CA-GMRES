#ifndef _KSP_SOLVER_H
#include "mmt_reader.h"
#include <time.h>
#include <papi.h>
#define _KSP_SOLVER_H

PetscErrorCode compute(PetscReal *avg_rtime, PetscInt *its, PetscReal **a, PetscInt *na, PetscReal *rnorm);

#endif //_KSP_SOLVER_H