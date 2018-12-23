PETSC_DIR=/usr/lib/petsc
MPI_DIR=/usr/include/openmpi
FLAGS=-Wall -std=iso9899:2011 -pedantic



all: main.o solver_util.o measure_time.o
	gcc -o execute main.o solver_util.o measure_time.o -lpetsc -lmpi -lpapi -lm

main.o: main.c solver_util.h measure_time.h
	gcc $(FLAGS) -g -ggdb -c main.c -I$(PETSC_DIR)/include -I$(MPI_DIR)

solver_util.o: solver_util.c solver_util.h
	gcc $(FLAGS) -g -ggdb -c solver_util.c -I$(PETSC_DIR)/include -I$(MPI_DIR)

measure_time.o: measure_time.c measure_time.h
	gcc $(FLAGS) -g -ggdb -c measure_time.c -I$(PETSC_DIR)/include -I$(MPI_DIR)

run: execute
	./execute -ksp_converged_reason

plot:
	gnuplot plot_history.gp

clear:
	rm *.o
	rm execute
