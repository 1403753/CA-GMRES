CC = g++
CFLAGS = -Wall -pedantic -std=c++14 -I$(HYPRE)/include -DBLASINT=int64_t -DMKL_ILP64 -m64 -I$(MKLROOT)/include -fopenmp
DFLAGS = -ggdb -g -O0
OBJ = gmres_ca.o gmres.o tsqr.o  arnoldi_ca.o matrix_reader.o block_cgs.o
DEP = $(wildcard *.tpp)

MKLROOT = /home/a01403753/intel/mkl
HYPRE = /home/a01403753/local/hypre

all: gmres_ca

gmres_ca: $(OBJ)
	$(CC) -o $@ $^ -lgsl \
	-Wl,--start-group \
	$(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a \
	$(MKLROOT)/lib/intel64/libmkl_gnu_thread.a \
	$(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group \
	-fopenmp -lpthread -lm -ldl \
	-lpapi \
	-lgfortran -L$(HYPRE)/lib -lHYPRE

gmres_ca.o: gmres_ca.cpp
	$(CC) $(CFLAGS) -c $< $(DFLAGS)
	
%.o: %.cpp %.hpp
	$(CC) $(CFLAGS) -c $< $(DFLAGS)
	
run:
	./gmres_ca

burn: gmres_ca
	./gmres_ca

vurn: gmres_ca
	valgrind ./gmres_ca

clean:
	rm -rf *.o gmres_ca
