CC = g++
CFLAGS = -Wall -pedantic -std=c++14 -I$(HYPRE)/include -DBLASINT=int64_t -DMKL_ILP64 -m64 -I$(MKLROOT)/include -fopenmp
# CFLAGS = -Wall -pedantic -std=c++14 -I$(HYPRE)/include -DMKL_ILP64 -m64 -I$(MKLROOT)/include -fopenmp
# CFLAGS = -Wall -pedantic -std=c++14 -I$(INC_ARPACK) -I$(HYPRE)/include -DBLASINT=int64_t -DMKL_ILP64 -m64 -I$(MKLROOT)/include -fopenmp
# CFLAGS = -Wall -pedantic -std=c++14 -I$(INC_ARPACK) -I$(HYPRE)/include -m64 -I$(MKLROOT)/include -fopenmp
DFLAGS = -ggdb -g -O0
OBJ = gmres.o spmv.o tsqr.o  arnoldi_ca.o matrix_reader.o block_cgs.o

MKLROOT = /home/a01403753/intel/mkl
# INC_ARPACK = /home/a01403753/local/include/arpackpp
HYPRE = /home/a01403753/local/hypre

##INC_VIENNACL = /home/a01403753/BSC/ViennaCL_rep/ViennaCL-1.7.1

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
	# -llapack \
	# -L/home/a01403753/local/lib -l:libarpack.so.2.0.0 \


gmres.o: gmres.cpp
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
	
# -L/usr/lib/x86_64-linux-gnu/libsuperlu.so -lsuperlu \
# -L/usr/lib/x86_64-linux-gnu/libarpack.so.2 -l:libarpack.so.2 \

#spmv.o: spmv.cpp spmv.hpp
#	$(CC) $(CFLAGS) -c $^
#-I$(INC_VIENNACL)

#tsqr.o: tsqr.cpp tsqr.hpp
#	$(CC) $(CFLAGS) -c $^

#block_cgs.o: block_cgs.cpp block_cgs.hpp
#	$(CC) $(CFLAGS) -c $^

#arnoldi_ca.o: arnoldi_ca.cpp arnoldi_ca.hpp
#	$(CC) $(CFLAGS) -c $^
#-I$(INC_LAPACK)

#gmres.o: gmres.cpp
#	$(CC) $(CFLAGS) -c $^ $(DFLAGS)
#-I$(INC_VIENNACL)