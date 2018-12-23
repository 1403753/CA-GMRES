static char help[] = "Solves a linear system with KSP.\n";

#include "main.h"

int main(int argc, char ** args){
	PetscErrorCode ierr;
	PetscReal 		 avg_rtime, rnorm, *a;
	PetscInt 			 i, its, na;
	char 					 fname[PETSC_MAX_PATH_LEN];
	char 					 fname2[PETSC_MAX_PATH_LEN];
	char 					 number[3];
	FILE					 *fp;
	
	if (argc < 2) {
    printf("missing filename\n");
    exit(1);
	}
	if (strlen(args[1]) >= sizeof(fname)) {
    printf("filename too long: %s\n", args[1]);
    exit(1);
	}
	/*
		process name input
	*/
	strcpy(fname, args[1]);
	strcpy(fname2, args[1]);
	memcpy(number, fname, 2);
	number[2] = '\0';
	fname2[0] = '_';
	fname2[1] = '_';
	
  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;

	/*
		solve with KSP and obtain values 
	*/
	ierr = compute(&avg_rtime, &its, &a, &na, &rnorm);CHKERRQ(ierr);
	
	/*
		Write data to files
	*/
	ierr = PetscFOpen(PETSC_COMM_SELF, fname, "w", &fp);CHKERRQ(ierr);
  if (!fp) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER, "Cannot open file");
	
	for(i = 0; i < na; ++i) {
		ierr = PetscFPrintf(PETSC_COMM_SELF, fp,"%ld, %f\n", i + 1, a[i]);CHKERRQ(ierr);
	}
	ierr = PetscFClose(PETSC_COMM_SELF, fp);CHKERRQ(ierr);
	
	ierr = PetscFOpen(PETSC_COMM_SELF, fname2, "a", &fp);CHKERRQ(ierr);
  if (!fp) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file");
	
	ierr = PetscFPrintf(PETSC_COMM_SELF, fp, "%s, %e,\t%ld,\t%e\n", number, avg_rtime, its, rnorm);CHKERRQ(ierr);
	ierr = PetscFClose(PETSC_COMM_SELF, fp);CHKERRQ(ierr);

	ierr = PetscFree(a);CHKERRQ(ierr);
	
	ierr = PetscFinalize();CHKERRQ(ierr);	
	return(ierr);
}



