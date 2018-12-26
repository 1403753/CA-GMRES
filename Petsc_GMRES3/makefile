PETSC_DIR = /usr/lib/petsc
MPI_DIR = /usr/include/openmpi
FLAGS = -Wall -std=iso9899:2011 -pedantic
CC = gcc
XX = g++
CXXFLASGS = -Wall -pedantic -std=c++14
PETSCFLAGS = --mca btl_base_warn_component_unused 0

all: ex10 ex27 main

main: main.o
	$(CXX) -o main main.o -lpetsc -lmpi -lpapi -lm

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -g -c main.cpp -I$(PETSC_DIR)/include -I$(MPI_DIR)

run: main
	mpirun -n 4 $(PETSCFLAGS) ./main -f e05r0000_petsc
	
ex10: ex10.o
	$(CC) -o ex10 ex10.o -lpetsc -lmpi -lpapi -lm

ex10.o: ex10.c
	$(CC) $(FLAGS) -g -c ex10.c -I$(PETSC_DIR)/include -I$(MPI_DIR)

ex27: ex27.o
	$(CC) -o ex27 ex27.o -lpetsc -lmpi -lpapi -lm

ex27.o: ex27.c
	$(CC) $(FLAGS) -g -c ex27.c -I$(PETSC_DIR)/include -I$(MPI_DIR)



ex1: ex10
	mpirun -n 0 --mca btl_base_warn_component_unused 0 ./ex10 -f0 e05r0000_petsc \
	-ksp_type gmres -pc_type none -mat_type mpiaij -ksp_rtol 1e-11 -f1 e05r0000_petsc \
	-ksp_converged_reason

ex2: ex27
	mpiexec -n 1 ./ex27 -f e05r0000_petsc -ksp_converged_reason

clean:
	rm -rf *.o ex1 ex2 main
