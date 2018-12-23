static char help[] = "PETSC program.\n\n";

#include <petscsys.h>
#include <petscmat.h>
#include <papi.h>
#include <time.h>
#include <math.h>

#include "solver_util.h"
#include "measure_time.h"

#define MY_COMM PETSC_COMM_SELF

void write_data_to_file(float value, char* matrix_name, char* value_type) {
		const char* solver_type = "indirect_solvers_pc";

		char filename[100] = "../data/";
		strcat(filename, solver_type);
		strcat(filename, "_");
		strcat(filename, value_type);
		strcat(filename, ".dat");

		FILE *f = fopen(filename, "a");
		if (!f) {
			printf("Could not open file \n");
			return;
		}
		fprintf(f, "%s      %f\n", matrix_name, value);
		fclose(f);
}

void write_residual_history_to_file(PetscReal* residual_history, int history_length, char* matrix_name) {
	char filename[100] = "residual_history_";
	strcat(filename, matrix_name);
	strcat(filename, ".dat");
	FILE *f = fopen(filename, "w");
	if (!f) {
		printf("Could not open file \n");
		return;
	}

	for (int i = 0; i < history_length; ++i) {
		fprintf(f, "%d    %f\n", i, log(residual_history[i]));
	}
}

int check_matrix(Mat mat, Vec b, Vec x, KSPType kspmethod, PCType pcmethod, char* matrix_name) {

	PetscErrorCode ierr;
	float sum_runtime = 0;
	int number_its = 5;

	PetscReal residual_history[10000]; //10000 is the default max number of iterations
	int history_length;

	for (int i = 0; i < number_its; ++i) {
		history_length = 10000;
		struct performance_info info;
		ierr = measure_time(&mat, &b, &x, kspmethod, pcmethod, &info, residual_history, &history_length); // Residual History is not necessary for direct solvers
		sum_runtime = sum_runtime + info.real_time;
	}

	float linsys_residual = residual(&mat, &b, &x);
	float avg_runtime = sum_runtime / number_its;

	write_data_to_file(avg_runtime, matrix_name, "runtime");
	write_data_to_file(linsys_residual, matrix_name, "residual");

	write_residual_history_to_file(residual_history, history_length, matrix_name);

	// temprorary
	printf("Avg runtime %s %f \n", matrix_name, avg_runtime);
	printf("Residual %s: %f\n", matrix_name, linsys_residual);

	return ierr;
}

int main(int argc,char **argv)
{
	PetscErrorCode ierr;

	ierr = PetscInitialize(&argc,&argv,(char *)0,help); CHKERRQ(ierr);

	Mat nasa, bc, gemat, goodwin;
	Vec b_nasa, b_bc, b_gemat, b_goodwin; //right hand side vectors
	Vec x_nasa, x_bc, x_gemat, x_goodwin; //solution vectors

	//
	//load matrices from files
	//
	ierr = load_matrix(&nasa, "../nasa4704/nasa4704_petsc"); CHKERRQ(ierr);
	ierr = load_matrix(&bc, "../bcsstk18/bcsstk18_petsc"); CHKERRQ(ierr);
	ierr = load_matrix(&gemat, "../gemat11/gemat11_petsc"); CHKERRQ(ierr);
	ierr = load_matrix(&goodwin, "../goodwin/goodwin_petsc"); CHKERRQ(ierr);


	//
	//create and initialize the right-hand-side vectors
	//

	//create random generator that should be used for random initialization
  PetscRandom rand;
  ierr = PetscRandomCreate(MY_COMM, &rand); CHKERRQ(ierr);
  ierr = PetscRandomSetType(rand, PETSCRAND); CHKERRQ(ierr);
	ierr = PetscRandomSetSeed(rand, time(NULL)); CHKERRQ(ierr);
	ierr = PetscRandomSeed(rand); CHKERRQ(ierr);

	//create right hand side vectors
	ierr = create_vec_random_from_mat(&nasa, &b_nasa, &rand); CHKERRQ(ierr);
	ierr = create_vec_random_from_mat(&bc, &b_bc, &rand); CHKERRQ(ierr);
	ierr = create_vec_random_from_mat(&gemat, &b_gemat, &rand); CHKERRQ(ierr);
	ierr = create_vec_random_from_mat(&goodwin, &b_goodwin, &rand); CHKERRQ(ierr);

	//some test output
	/*
	PetscScalar min, max;
	ierr = VecMax(b_nasa, NULL, &max); CHKERRQ(ierr); //vector, index of max value, max value
	ierr = VecMin(b_nasa, NULL, &min); CHKERRQ(ierr);
	printf("Max vector entry: %f\n", max);
	printf("Min vector entry: %f\n", min);
	ierr = VecMax(b_goodwin, NULL, &max); CHKERRQ(ierr); //vector, index of max value, max value
	ierr = VecMin(b_goodwin, NULL, &min); CHKERRQ(ierr);
	printf("Max vector entry: %f\n", max);
	printf("Min vector entry: %f\n", min);
	*/



	//
	//create solution vectors
	//

	//duplicating the right-hand-side vectors (because length and type have to match)
	ierr = VecDuplicate(b_nasa, &x_nasa); CHKERRQ(ierr);
	ierr = VecDuplicate(b_bc, &x_bc); CHKERRQ(ierr);
	ierr = VecDuplicate(b_gemat, &x_gemat); CHKERRQ(ierr);
	ierr = VecDuplicate(b_goodwin, &x_goodwin); CHKERRQ(ierr);

	//
	//solve linear system
	//
	//linsys_solve(Mat* A, Vec* b, Vec* x, KSPType kspmethod, PCType pcmethod)

	/*
	TO DO:
	get residual history ( KSPSetResidualHistory(), KSPGetResidualHistory() )
	output results to files

	Fix LU factorisation
	*/

	ierr = check_matrix(nasa, b_nasa, x_nasa, KSPCG, PCJACOBI, "nasa_jacobi"); CHKERRQ(ierr); //conjugate gradient
	ierr = check_matrix(nasa, b_nasa, x_nasa, KSPCG, PCICC, "nasa_incompl-cholesky"); CHKERRQ(ierr); //conjugate gradient
	ierr = check_matrix(bc, b_bc, x_bc, KSPCG, PCJACOBI, "bcsstk18_jacobi"); CHKERRQ(ierr);
	ierr = check_matrix(bc, b_bc, x_bc, KSPCG, PCICC, "bcsstk18_incompl_cholesky"); CHKERRQ(ierr);

	ierr = check_matrix(gemat, b_gemat, x_gemat, KSPBCGS, PCJACOBI, "gemat_jacobi"); CHKERRQ(ierr); //bi conjugate gradient stabilized
//	ierr = check_matrix(gemat, b_gemat, x_gemat, KSPBCGS, PCILU, "gemat_incompl_lu"); CHKERRQ(ierr); //bi conjugate gradient stabilized
	ierr = check_matrix(goodwin, b_goodwin, x_goodwin, KSPBCGS, PCJACOBI, "goodwin_jacobi"); CHKERRQ(ierr);
//	ierr = check_matrix(goodwin, b_goodwin, x_goodwin, KSPBCGS, PCILU, "goodwin_incompl_lu"); CHKERRQ(ierr);



	//
	//free memory
	//

  ierr = PetscRandomDestroy(&rand); CHKERRQ(ierr);

	ierr = VecDestroy(&b_nasa); CHKERRQ(ierr);
	ierr = VecDestroy(&b_bc); CHKERRQ(ierr);
	ierr = VecDestroy(&b_gemat); CHKERRQ(ierr);
	ierr = VecDestroy(&b_goodwin); CHKERRQ(ierr);

	ierr = VecDestroy(&x_nasa); CHKERRQ(ierr);
	ierr = VecDestroy(&x_bc); CHKERRQ(ierr);
	ierr = VecDestroy(&x_gemat); CHKERRQ(ierr);
	ierr = VecDestroy(&x_goodwin); CHKERRQ(ierr);

	ierr = MatDestroy(&nasa); CHKERRQ(ierr);
	ierr = MatDestroy(&bc); CHKERRQ(ierr);
	ierr = MatDestroy(&gemat); CHKERRQ(ierr);
	ierr = MatDestroy(&goodwin); CHKERRQ(ierr);

	ierr = PetscFinalize(); CHKERRQ(ierr);

	return 0;
}
