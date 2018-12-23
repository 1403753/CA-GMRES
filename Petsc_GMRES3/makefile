PETSC_DIR=/usr/lib/petsc
MPI_DIR=/usr/include/openmpi
FLAGS=-Wall -std=iso9899:2011 -pedantic

all: ex10 ex27

ex10: ex10.o
	gcc -o ex10 ex10.o -lpetsc -lmpi -lpapi -lm

ex10.o: ex10.c
	gcc $(FLAGS) -g -c ex10.c -I$(PETSC_DIR)/include -I$(MPI_DIR)

ex27: ex27.o
	gcc -o ex27 ex27.o -lpetsc -lmpi -lpapi -lm

ex27.o: ex27.c
	gcc $(FLAGS) -g -c ex27.c -I$(PETSC_DIR)/include -I$(MPI_DIR)
	
ex1: ex10
	# mpirun -n 0 --mca orte_base_help_aggregate 0 ./ex10 -f0 e05r0000_petsc -ksp_type cg -pc_type asm -pc_asm_type basic -sub_pc_type icc -mat_type sbaij 
	mpirun -n 0 --mca btl_base_warn_component_unused 0 ./ex10 -f0 e05r0000_petsc -ksp_type gmres -pc_type none -mat_type mpiaij -ksp_rtol 1e-11

ex2: ex27
	mpiexec -n 1 ./ex27 -f e05r0000_petsc -ksp_converged_reason

clear:
	rm *.o
	rm execute
