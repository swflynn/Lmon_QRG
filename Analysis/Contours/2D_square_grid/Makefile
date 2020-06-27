OBJS=Tip Lmon main
FC=gfortran
#FFLAGS=-O2
FFLAGS= -g -O0 -Wall -Wextra -pedantic -fcheck=all -fbacktrace
all: Tip Lmon main
Tip:
	$(FC) $(FFLAGS) -c TIP4P.f90
Lmon:
	$(FC) $(FFLAGS) -c unif_grid_pot.f90 -llapack
main:
	$(FC) $(FFLAGS) unif_grid_pot.o TIP4P.o -llapack

clean:
	rm -rf *.o *.mod *.out $(OBJS)
