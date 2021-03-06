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
	
	#############
	#  threads  #
	#############

	./main 1 threads 10 5 4 '\shortstack{\footnotesize\,Runtime\,per\,kernel,\,relative\,to\,CA-GMRES(5,4),\,for\,the\,unpreconditioned\,matrix\\\footnotesize\,Xenon2\,and\,restart\,length\,20\,for\,various\,numbers\,of\,threads}' \
	xenon2 2>&1 | tee log.txt
	gnuplot ../gnuplot/runtimes_threads.gnu
	
	# ./main 1 threads 100 5 6 '\shortstack{\footnotesize\,Runtime\,per\,kernel,\,relative\,to\,CA-GMRES(5,6),\,for\,ILU(0)\,preconditioned\,matrix\\\footnotesize\,Sherman3\,and\,restart\,length\,30\,for\,various\,numbers\,of\,threads}' \
	# sherman3 2>&1 | tee log.txt
	# gnuplot ../gnuplot/runtimes_threads.gnu

	#############
	#  speedup  #
	#############

	# ./main 2 speedup 10 5 12 '\shortstack{\footnotesize\,Runtime\,per\,kernel,\,relative\,to\,CA-GMRES(5,12),\,for\,all\,test\,matrices,\\\footnotesize\,using\,8\,threads\,and\,restart\,length\,60}' \
	# bcsstk18 
	# dwb512
	# pwtk
	# watt1 bmw7st1 xenon2 

	# gnuplot ../gnuplot/runtimes_matrices.gnu

	###############
	#  residuals  #
	###############

	##############
	#  SHERMAN3  #
	##############	

	# ./main 2 residuals sherman3 30 3 10 5 6 '\shortstack{Sherman3}'
	# gnuplot ../gnuplot/sherman3/residuals_sherman3.gnu

	#######################
	#  SHERMAN3 + ILU(0)  #
	#######################

	# ./main 2 residuals sherman3 30 3 10 5 6 '\shortstack{Sherman3}'
	# gnuplot ../gnuplot/sherman3/residuals_sherman3_ilu0.gnu

	#########
	#  bmw  #
	#########

	# ./main 2 residuals bmw7st1 20 5 4 10 2 '\shortstack{bmw}'
	# gnuplot ../gnuplot/bmw/residuals_bmw.gnu
	
	###########
	#  WATT1  #
	###########

	# ./main 2 residuals watt1 30 3 10 5 6 '\shortstack{Watt1}'
	# gnuplot ../gnuplot/Watt1_pcnone_s3_t10_s5_t6/residuals_watt1.gnu

	####################
	#  WATT1 + ILU(0)  #
	####################

	# ./main 2 residuals watt1 30 5 6 15 2 '\shortstack{Watt1}'
	# gnuplot ../gnuplot/Watt1_pcnone_s3_t10_s5_t6/residuals_watt1_ilu0.gnu

	############
	#  Xenon2  #
	############

	# ./main 2 residuals xenon2 20 5 4 10 2 '\shortstack{Xenon2}'
	# gnuplot ../gnuplot/Xenon2/residuals_xenon2.gnu

	#####################
	#  Xenon2 + ILU(0)  #
	#####################

	# ./main 2 residuals xenon2 20 5 4 10 2 '\shortstack{Xenon2}'
	# gnuplot ../gnuplot/Xenon2/residuals_xenon2_ilu0.gnu	
	
	##########
	#  pwtk  #
	##########

	# ./main 2 residuals pwtk 20 5 4 10 2 '\shortstack{pwtk}'
	# gnuplot ../gnuplot/pwtk/residuals_pwtk.gnu

	##############
	#  bcsstk18  #
	##############
	
	# ./main 2 residuals bcsstk18 60 5 12 10 6 '\shortstack{bcsstk08}'
	# gnuplot ../gnuplot/bcsstk18/residuals_bcsstk18.gnu
	
	##########
	#  dmat  #
	##########	
	
	# ./main 2 residuals dmat1 60 15 4 20 3 '\shortstack{\footnotesize\,1e+04\,$\times$\,1e+04\,diagonal,\,logspace\,eigs,\,cond\,1e+05\\\footnotesize\,Residual\,2-norm,\,log\,scale}'
	# gnuplot ../gnuplot/dmat/residuals_dmat.gnu
	# ./main 2 residuals dmat2 60 15 4 20 3 '\shortstack{\footnotesize\,1e+04\,$\times$\,1e+04\,diagonal,\,logspace\,eigs,\,cond\,1e+10\\\footnotesize\,Residual\,2-norm,\,log\,scale}'
	# gnuplot ../gnuplot/dmat/residuals_dmat.gnu
	# ./main 2 residuals dmat3 60 15 4 20 3 '\shortstack{\footnotesize\,1e+04\,$\times$\,1e+04\,diagonal,\,logspace\,eigs,\,cond\,1e+15\\\footnotesize\,Residual\,2-norm,\,log\,scale}'
	# gnuplot ../gnuplot/dmat/residuals_dmat.gnu

burn: main
	./main

	
vurn: main
	valgrind ./main

	
clean:
	rm -rf *.o main
