CC = g++
CFLAGS = -Wall -pedantic -std=c++14
INC = -I/home/a01403753/local/include -DBLASINT=int64_t -DMKL_ILP64 -m64 -I$(MKLROOT)/include -fopenmp
LIB = /home/a01403753/local/lib

PETSCFLAGS = --mca btl_base_warn_component_unused 0

PETSC_DIR = /usr/lib/petsc
MPI_DIR = /usr/include/openmpi

DFLAGS = -g -O0 
#-ggdb

SRC         := $(wildcard *.cpp)
OBJ         := $(SRC:.cpp=.o)
#DEP         := $(OBJ:.o=.tpp)

MKLROOT = /home/a01403753/intel/mkl

all: main

main: $(OBJ)
	$(CC) -o $@ $^ -lgsl \
	-Wl,--start-group \
	$(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a \
	$(MKLROOT)/lib/intel64/libmkl_gnu_thread.a \
	$(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group \
	-lgomp -lpthread -lm -ldl \
	-lpapi \
	-lgfortran \
	-lpetsc -lmpi
	# -L$(LIB) -l:libmetis.a

main.o: main.cpp
	$(CC) $(CFLAGS) $(INC) -c $< $(DFLAGS) -I$(PETSC_DIR)/include -I$(MPI_DIR)
	
%.o: %.cpp %.hpp
	$(CC) $(CFLAGS) $(INC) -c $< $(DFLAGS)
	
run:
	mpirun -n 1 ./main $(PETSCFLAGS)

burn: main
	./main
	gnuplot ../gnuplot/residuals.gnu
	
vurn: main
	valgrind ./main

	
clean:
	rm -rf *.o main
