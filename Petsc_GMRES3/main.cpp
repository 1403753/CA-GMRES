#include <petscksp.h>
#include <iostream>

static char help[] = "Reads a PETSc matrix and vector from a file and solves a linear system.\n\
This version first preloads and solves a small system, then loads \n\
another (larger) system and solves it as well.  This example illustrates\n\
preloading of instructions with the smaller system so that more accurate\n\
performance monitoring can be done with the larger one (that actually\n\
is the system of interest).  See the 'Performance Hints' chapter of the\n\
users manual for a discussion of preloading.  Input parameters include\n\
  -f0 <input_file> : first file to load (small system)\n\
  -f1 <input_file> : second file to load (larger system)\n\n\
  -nearnulldim <0> : number of vectors in the near-null space immediately following matrix\n\n\
  -trans  : solve transpose system instead\n\n";

int main(int argc,char **args)
{
  // KSP            ksp;             /* linear solver context */
  Mat            A;               /* matrix */
  // Vec            x,b,u;           /* approx solution, RHS, exact solution */
  PetscViewer    fd;              /* viewer */
  char           file[PETSC_MAX_PATH_LEN];     /* input file name */
  PetscBool      flg,trans=PETSC_FALSE;
  PetscErrorCode ierr;
  // PetscInt       its,num_numfac,m,n,M,nearnulldim = 0;
  // PetscReal      norm;
  PetscMPIInt    rank;
  char           initialguessfilename[PETSC_MAX_PATH_LEN];

  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

  /*
     Determine files from which we read the two linear systems
     (matrix and right-hand-side vector).
  */
  ierr = PetscOptionsGetString(NULL,NULL,"-f",file,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);  
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate binary file with the -f option");

  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file,FILE_MODE_READ,&fd);CHKERRQ(ierr);

	
	ierr  = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr  = MatSetType(A,MATMPIAIJ);CHKERRQ(ierr);
	
	
  ierr = PetscFinalize();
  return ierr;
}
