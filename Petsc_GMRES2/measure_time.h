#include <time.h>

#include <petscmat.h>
#include <petscksp.h>
#include <petscpctypes.h>

struct performance_info {
  float real_time;
  float mflops;
};

//struct performance_info measure_time(void (*fun) (Mat*, Vec*, Vec*), Mat* A, Vec* b, Vec* x);
int measure_time(Mat* A, Vec* b, Vec* x, KSPType kspmethod, PCType pcmethod, struct performance_info *info, PetscReal* residualHistory, int* history_length);
