CC = g++
CFLAGS = -pedantic -Wall -std=c++14

all: dmat_generator

dmat_generator: dmat_generator.o
	$(CC) -o $@ $^
	
dmat_generator.o: dmat_generator.cpp
	$(CC) $(CFLAGS) -c $<
	
run:
	./dmat_generator
	