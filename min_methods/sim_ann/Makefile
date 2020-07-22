FC=gfortran
FFLAGS=-O2
#FFLAGS= -g -O0 -Wall -Wextra -pedantic -fcheck=all -fbacktrace
all: grid main
grid:
	$(FC) $(FFLAGS) -c annealing.f90 -llapack
main:
	$(FC) $(FFLAGS) annealing.o -llapack

clean:
	rm -rf *.o *.mod *.out 
