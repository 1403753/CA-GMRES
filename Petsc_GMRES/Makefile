CFLAGS = -Wall -pedantic -std=iso9899:2011
DEPS = main.h ksp_solver.h mmt_reader.h
OBJ = main.o ksp_solver.o mmt_reader.o
DIR_LIB = /home/a01403753/exercises/lib
LIB_PAPI = $(DIR_LIB)/PAPI/lib
INC_PAPI = $(DIR_LIB)/PAPI/include

all: main

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

main: $(OBJ) $(DEPS) chkopts 
	-${CLINKER} -I${PETSC_DIR}/include -I$(INC_PAPI) $(CFLAGS) -o $@ $(OBJ) -L$(LIB_PAPI) -l:libpapi_.a ${PETSC_KSP_LIB}
	${RM} $(OBJ)

mmt_reader.o: mmt_reader.c mmt_reader.h chkopts
	-${PETSC_COMPILE} $(CFLAGS) -c $<

ksp_solver.o: ksp_solver.c ksp_solver.h chkopts
	-${PETSC_COMPILE} $(CFLAGS) -c $<
	
test:
	./main 01test.txt -mmt mmt/1138_bus.mtx -sym -random_b -ksp_type preonly -pc_type cholesky -ksp_final_residual -ksp_converged_reason
	./main 02test.txt -mmt mmt/1138_bus.mtx -sym -random_b -ksp_type preonly -pc_type lu -ksp_final_residual -ksp_converged_reason
	./main 03test.txt -mmt mmt/1138_bus.mtx -sym -random_b -ksp_type bcgs -pc_type none -ksp_final_residual -ksp_converged_reason
	./main 04test.txt -mmt mmt/1138_bus.mtx -sym -random_b -ksp_type bcgs -pc_type jacobi -ksp_final_residual -ksp_converged_reason
	gnuplot gnuplot_test_hw10.txt
	
run:
# printf "1a) nasa4704\n"
	./main 01out.txt -mmt mmt/nasa4704.mtx -sym -random_b -ksp_type preonly -pc_type cholesky -ksp_final_residual -ksp_converged_reason
# printf "1a) bcsstk18\n"
	./main 02out.txt -mmt mmt/bcsstk18.mtx -sym -random_b -ksp_type preonly -pc_type cholesky -ksp_final_residual -ksp_converged_reason
# printf "1b) gemat11\n"
	./main 03out.txt -mmt mmt/gemat11.mtx -random_b -ksp_type preonly -pc_type lu -ksp_final_residual -ksp_converged_reason
# printf "1b ) goodwin\n"
	./main 04out.txt -mmt mmt/goodwin.mtx -random_b -ksp_type preonly -pc_type lu -ksp_final_residual -ksp_converged_reason
# printf "2a) nasa4704\n"
	./main 05out.txt -mmt mmt/nasa4704.mtx -sym -random_b -ksp_type cg -pc_type none -ksp_final_residual -ksp_converged_reason
# printf "2a) bcsstk18\n"
	./main 06out.txt -mmt mmt/bcsstk18.mtx -sym -random_b -ksp_type cg -pc_type none -ksp_final_residual -ksp_converged_reason
# printf "2b) gemat11\n"
	./main 07out.txt -mmt mmt/gemat11.mtx -random_b -ksp_type bcgs -pc_type none -ksp_final_residual -ksp_converged_reason
# printf "2b) goodwin\n"
	./main 08out.txt -mmt mmt/goodwin.mtx -random_b -ksp_type bcgs -pc_type none -ksp_final_residual -ksp_converged_reason
# printf "3a) nasa4704\n"
	./main 09out.txt -mmt mmt/nasa4704.mtx -sym -random_b -ksp_type cg -pc_type jacobi -ksp_final_residual -ksp_converged_reason
# printf "3a) bcsstk18\n"
	./main 10out.txt -mmt mmt/bcsstk18.mtx -sym -random_b -ksp_type cg -pc_type icc -ksp_final_residual -ksp_converged_reason
# printf "3a) nasa4704\n"
	./main 11out.txt -mmt mmt/nasa4704.mtx -sym -random_b -ksp_type cg -pc_type icc -ksp_final_residual -ksp_converged_reason
# printf "3a) bcsstk18\n"
	./main 12out.txt -mmt mmt/bcsstk18.mtx -sym -random_b -ksp_type cg -pc_type jacobi -ksp_final_residual -ksp_converged_reason
# printf "3b) gemat11\n"
	./main 13out.txt -mmt mmt/gemat11.mtx -random_b -ksp_type bcgs -pc_type jacobi -ksp_final_residual -ksp_converged_reason
# printf "3b) goodwin\n"
	./main 14out.txt -mmt mmt/goodwin.mtx -random_b -ksp_type bcgs -pc_type ilu -ksp_final_residual -ksp_converged_reason
# printf "3b) gemat11\n"
	./main 15out.txt -mmt mmt/gemat11.mtx -random_b -ksp_type bcgs -pc_type ilu -ksp_final_residual -ksp_converged_reason
# printf "3b) goodwin\n"
	./main 16out.txt -mmt mmt/goodwin.mtx -random_b -ksp_type bcgs -pc_type jacobi -ksp_final_residual -ksp_converged_reason
	
	gnuplot gnuplot_hw10.txt
	
clean_all:
	rm  -rf *.o main