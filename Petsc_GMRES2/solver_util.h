#include <petscmat.h>
#include <petscksp.h>
#include <petscpctypes.h>

int load_matrix(Mat* A, char* file);
int linsys_solve(Mat* A, Vec* b, Vec* x, KSPType kspmethod, PCType pcmethod);
PetscReal residual(Mat* A, Vec* b, Vec* x);
int create_vec_random_from_mat(Mat* A, Vec* b, PetscRandom* r);
