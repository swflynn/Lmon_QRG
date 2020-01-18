OBJS = Tip SSobol SS_std Lmon main
FC = gfortran

all: SSobol Tip SS_std Lmon main

Tip:
	$(FC) -c TIP4P.f90

SSobol:
	$(FC) -c ssobol.f

SS_std:
	$(FC) -c ssobol_stdnormal.f90

Lmon:
	$(FC) -c -fopenmp water_Lmon.f90 -llapack

main:
	$(FC) -O -Wall -fcheck=all -g -fbacktrace -fopenmp water_Lmon.o TIP4P.o ssobol.o ssobol_stdnormal.o  mbpol/libmbpol.a   -lstdc++ -llapack

clean:
	rm -rf *.o *.mod *.out $(OBJS)
