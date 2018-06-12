/*
 * matrix_reader.cpp
 *
 *  Created on: 06.05.2018
 *      Author: Robert
 */
#include "matrix_reader.hpp"
matrix_reader::matrix_reader() {
	
}

void matrix_reader::read_matrix_from_file() {
	/*
	char buf[PETSC_MAX_PATH_LEN];
	PetscInt       i, n, nnz, nz, nztemp, nzmax, col, row;
	PetscScalar	   value;
	PetscReal	   	 r_value;
	FILE*          file;
	PetscInt 	   	 ierr;
	nz=1;
	nztemp=-1;
	nzmax=1;
	
	ierr = PetscFOpen(PETSC_COMM_WORLD,fin,"r",&file);CHKERRQ(ierr);
		
	// Print the first line of the file
	fgets(buf, PETSC_MAX_PATH_LEN-1, file);
	fscanf(file, "%d %d %d\n", &n, &n, &nnz);
	
	if(nnz<=0) SETERRQ(PETSC_COMM_WORLD,1,"Matrix Market Converter : you must verify the format of entry file\n");

	minfo->n=n;
	minfo->m=n;
	minfo->nnz=nnz;

	PetscPrintf(PETSC_COMM_WORLD,"Matrix properties : m = %d, n = %d, nnz = %d\n", n, n, nnz);CHKERRQ(ierr);
	
	ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, n, n, PETSC_DEFAULT, PETSC_NULL, A);CHKERRQ(ierr);
	
	for (i = 0; i < n; ++i) {
		ierr = MatSetValue(*A, i, i, 0, INSERT_VALUES);CHKERRQ(ierr);
	}
	
	fscanf(file,"%d %d %le\n",&row,&col,&r_value);
	row = row-1; col = col-1 ;
	if(nztemp!=col){
		nz=1;
		nztemp=col;
	} else {
		nz++;
	}
	if(nz>nzmax) {
		nzmax=nz;
	}

	value=(PetscScalar)r_value;
	*A_min = value;
	*A_max = value;
	ierr = MatSetValues(*A,1,&row,1,&col,&value,INSERT_VALUES);CHKERRQ(ierr);
		
	//Matrix reading
	for (i = 1; i < nnz; ++i) {
		fscanf(file,"%d %d %le\n",&row,&col,&r_value);
		row = row-1; col = col-1 ;
		if(nztemp!=col){
			nz=1;
			nztemp=col;
		} else {
			nz++;
		}
		if(nz>nzmax){
			nzmax=nz;
		}
		value = (PetscScalar)r_value;
		if(value < *A_min) *A_min = value;
		if(value > *A_max) *A_max = value;
		ierr = MatSetValues(*A,1,&row,1,&col,&value,INSERT_VALUES);CHKERRQ(ierr);
	}
	PetscPrintf(PETSC_COMM_WORLD,"Maximum number of NNZ on a line : %d\n",nz);
	PetscPrintf(PETSC_COMM_WORLD,"A_min: %f\nA_max: %f\n", *A_min, *A_max);
	
	
	//Matrix assembly
	PetscPrintf(PETSC_COMM_WORLD,"Assembling matrix within PETSc.\n");
	ierr = MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(	*A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"Finished matrix assembly.\n");
	
	fclose(file);
	*/
}


matrix_reader::~matrix_reader() {
	
}

