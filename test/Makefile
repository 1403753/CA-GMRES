CC = g++
CFLAGS = -Wall -pedantic -std=c++14 
LIB = /home/a01403753/local/lib
INCLUDE = /home/a01403753/local/include
INC = -DBLASINT=int64_t -DMKL_ILP64 -m64 -I$(MKLROOT)/include -fopenmp

DFLAGS = -g -O0 
#-ggdb

SRC         := $(wildcard *.cpp)
OBJ         := $(SRC:.cpp=.o)
#DEP         := $(OBJ:.o=.tpp)

MKLROOT = /home/a01403753/intel/mkl

all: metis

metis: $(OBJ)
	$(CC) -o $@ $^ \
	-lpapi \
	-lgfortran \
	-L$(LIB) -l:libmetis.a \
	-L$(LIB) -l:libigraph.a \
	-L$(LIB) -lcxsparse \
	-Wl,--start-group \
	$(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a \
	$(MKLROOT)/lib/intel64/libmkl_gnu_thread.a \
	$(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group \
	-lgomp -lpthread -lm -ldl

metis.o: metis.cpp
	$(CC) $(CFLAGS) $(INC) -I$(INCLUDE) -I$(INCLUDE)/igraph -c $< $(DFLAGS)

%.o: %.cpp %.hpp
	$(CC) $(CFLAGS) $(INC) -c $< $(DFLAGS)

run:
	./metis

burn: metis
	./metis

vurn: metis
	valgrind ./metis

clean:
	rm -rf *.o metis
