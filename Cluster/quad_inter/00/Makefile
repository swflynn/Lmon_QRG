#==============================================================================#
FC=gfortran
#FFLAGS=-O2
FFLAGS= -g -O0 -Wall -Wextra -pedantic -fcheck=all -fbacktrace
#==============================================================================#
all: spectra
#==============================================================================#
#                       Compile all Fortran Files
#==============================================================================#
%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $^
#==============================================================================#
spectra: TIP4P.o gen_hermite_rule.o main.o
	$(FC) -o $@ $^ -llapack
#==============================================================================#
clean:
	rm -rf *.o *.mod spectra
