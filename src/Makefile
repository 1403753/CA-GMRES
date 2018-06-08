CC = g++
CFLAGS = -Wall -pedantic -std=c++14 -DMKL_ILP64 -m64 -I$(MKLROOT)/include -I$(INC_PAPI)
DFLAGS = -ggdb -g -O0
DEPS = mpk.hpp tsqr.hpp block_cgs.hpp arnoldi_ca.hpp spmv.hpp
OBJ =  spmv.o mpk.o tsqr.o block_cgs.o arnoldi_ca.o gmres.o 
INC_LAPACK = /usr/lib/x86_64-linux-gnu/lapack
INC_VIENNACL = /home/a01403753/BSC/ViennaCL_rep/ViennaCL-1.7.1
MKLROOT = /home/a01403753/intel/mkl
LIB_PAPI = /usr/lib/x86_64-linux-gnu
INC_PAPI = /usr/include

all: gmres_ca

gmres_ca: $(OBJ)
	$(CC) -o $@ $^ -lgsl \
	-Wl,--start-group \
	$(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a \
	$(MKLROOT)/lib/intel64/libmkl_gnu_thread.a \
	$(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group \
	-lgomp -lpthread -lm -ldl \
	-L$(LIB_PAPI) -lpapi

%.o: %.cpp %.hpp
	$(CC) $(CFLAGS) -c $< $(DFLAGS)

run:
	./gmres_ca
	./gmres_ca
	./gmres_ca
	./gmres_ca
	./gmres_ca
	./gmres_ca

burn: gmres_ca
	./gmres_ca

vurn: gmres_ca
	valgrind ./gmres_ca

clean:
	rm -rf *.o gmres_ca *.gch



#spmv.o: spmv.cpp spmv.hpp
#	$(CC) $(CFLAGS) -c $^
#-I$(INC_VIENNACL)

#mpk.o: mpk.cpp mpk.hpp
#	$(CC) $(CFLAGS) -c $^

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