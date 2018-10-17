CC = g++
CFLAGS = -Wall -pedantic -std=c++14 -I$(HYPRE)/include -DBLASINT=int64_t -DMKL_ILP64 -m64 -I$(MKLROOT)/include -fopenmp

DFLAGS = -ggdb -g -O0

SRC         := $(wildcard *.cpp)
OBJ         := $(SRC:.cpp=.o)
#DEP         := $(OBJ:.o=.tpp)

MKLROOT = /home/a01403753/intel/mkl
HYPRE = /home/a01403753/local/hypre

all: main

main: $(OBJ)
	$(CC) -o $@ $^ -lgsl \
	-Wl,--start-group \
	$(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a \
	$(MKLROOT)/lib/intel64/libmkl_gnu_thread.a \
	$(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group \
	-fopenmp -lpthread -lm -ldl \
	-lpapi \
	-lgfortran -L$(HYPRE)/lib -lHYPRE
	
main.o: main.cpp
	$(CC) $(CFLAGS) -c $< $(DFLAGS)
	
%.o: %.cpp %.hpp
	$(CC) $(CFLAGS) -c $< $(DFLAGS)
	
run:
	./main

burn: main
	./main

vurn: main
	valgrind ./main

clean:
	rm -rf *.o main
