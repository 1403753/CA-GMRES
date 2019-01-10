CC = g++
CFLAGS = -Wall -pedantic -std=c++14
INC = -I/home/a01403753/local/include -DBLASINT=int64_t -DMKL_ILP64 -m64 -I$(MKLROOT)/include -fopenmp
LIB = /home/a01403753/local/lib

PETSCFLAGS = --mca btl_base_warn_component_unused 0

PETSC_DIR = /usr/lib/petsc
MPI_DIR = /usr/include/openmpi

DFLAGS = -g -O0
# -ggdb

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
	-lgfortran
#	-lpetsc -L$(LIB) -l:libmetis.a

main.o: main.cpp
	$(CC) $(CFLAGS) $(INC) -c $< $(DFLAGS)
	
%.o: %.cpp %.hpp
	$(CC) $(CFLAGS) $(INC) -c $< $(DFLAGS)
	
run:
	mpirun -n 1 ./main $(PETSCFLAGS)


plot: main
	./main 8 speedup 10 5 12 '\shortstack{\footnotesize\,Runtime\,per\,kernel,\,relative\,to\,CA-GMRES(5,12),\,for\,all\,test\,matrices,\\\footnotesize\,using\,8\,threads\,and\,restart\,length\,60}' \
	dwb512
	# pwtk
	# bcsstk18 watt1 bmw7st1 xenon2 
	gnuplot ../gnuplot/runtimes_matrices.gnu

	#./main 2 residuals dmat1 120 15 8 20 6 \\shortstack{\\footnotesize\\,1e+04\\,$\\times$\\,1e+04\\,diagonal,\\,logspace\\,eigs,\\,cond\\,1e+05\\\\\\footnotesize\\,Residual\\,2-norm,\\,log\\,scale}
	#gnuplot ../gnuplot/residuals_dmat.gnu
	#./main 2 residuals dmat2 120 15 8 20 6 \\shortstack{\\footnotesize\\,1e+04\\,$\\times$\\,1e+04\\,diagonal,\\,logspace\\,eigs,\\,cond\\,1e+10\\\\\\footnotesize\\,Residual\\,2-norm,\\,log\\,scale}
	#gnuplot ../gnuplot/residuals_dmat.gnu
	#./main 2 residuals dmat3 120 15 8 20 6 \\shortstack{\\footnotesize\\,1e+04\\,$\\times$\\,1e+04\\,diagonal,\\,logspace\\,eigs,\\,cond\\,1e+15\\\\\\footnotesize\\,Residual\\,2-norm,\\,log\\,scale}
	#gnuplot ../gnuplot/residuals_dmat.gnu
	

burn: main
	./main
	# gnuplot ../gnuplot/runtimes_matrices.gnu
	gnuplot ../gnuplot/residuals_dmat.gnu
	# gnuplot ../gnuplot/residuals_dmat1_2.gnu
	# gnuplot ../gnuplot/residuals_dmat2.gnu
	# gnuplot ../gnuplot/residuals_dmat2_2.gnu
	# gnuplot ../gnuplot/residuals_dmat3.gnu
	# gnuplot ../gnuplot/residuals_dmat3_2.gnu
	# gnuplot ../gnuplot/residuals_bmw.gnu
	
vurn: main
	valgrind ./main

	
clean:
	rm -rf *.o main
