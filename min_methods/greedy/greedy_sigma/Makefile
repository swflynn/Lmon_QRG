FC=gfortran
FFLAGS=-O2
#FFLAGS= -g -O0 -Wall -Wextra -pedantic -fcheck=all -fbacktrace
all: grid main
grid:
	$(FC) $(FFLAGS) -c greedy_sigma.f90 -llapack
main:
	$(FC) $(FFLAGS) greedy_sigma.o -llapack

clean:
	rm -rf *.o *.mod *.out 
