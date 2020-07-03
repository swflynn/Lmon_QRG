FC=gfortran
FFLAGS=-O2
#FFLAGS= -g -O0 -Wall -Wextra -pedantic -fcheck=all -fbacktrace
all: Tip Lmon main
Tip:
	$(FC) $(FFLAGS) -c TIP4P.f90
Lmon:
	$(FC) $(FFLAGS) -c main_Lmon_grid.f90 -llapack
main:
	$(FC) $(FFLAGS) main_Lmon_grid.o TIP4P.o -llapack

clean:
	rm -rf *.o *.mod *.out 
