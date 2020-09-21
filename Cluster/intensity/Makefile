#==============================================================================#
#                          TIP4P/Mbpol Compile/Flags
#The MBPOL dipole surface uses the MKL libraries, compile with intel on HPC
#==============================================================================#
FC=ifort
FFLAGS=-O2
MKLFLAGS=-L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
#==============================================================================#
all: intensity
#==============================================================================#
#                          Compile Fortran Files
#==============================================================================#
%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $^
#==============================================================================#
#                               DGB Cluster
#==============================================================================#
intensity: mb-md/libmb-md.a gen_hermite_rule.o main_Lmon_intensity.o
	$(FC) -o $@ $^ src_mb-mu/mb-md/libmb-md.a $(FFLAGS) $(MKLFLAGS) -lstdc++ -llapack
#==============================================================================#
clean:
	rm -rf *.o *.mod intensity
