!=================20==================40==================60==================80
!                    nD Vibrational EigenSpectra monomer
!==============================================================================!
!Vibrational EigenSpectra calculations for a given potentials QRG
!All equations written for general d-dimensional case
!Basis Functions:
!             phi(r):=(2*alpha_i/pi)^(d/4)*exp[-alpha_i(r-r_i)^2]
!Generates alpha (inverse gaussian width i.e. small alpha=broad gaussian)
!using nearest neighbor (assumes a Quasi-Regular Grid as input)
!Needs gen_hermite_rule.f90 code (Gauss-Hermite-Quadriture) for potential eval.
!Needs LLAPACK to solve the generalized eigenvalue problem
!==============================================================================!
!       Modified:
!   20 April 2020
!       Author:
!   Shane Flynn
!==============================================================================!
!Begin Developing monomer calculations
!Grid is generated in d1 subspace
!==============================================================================!
program main
use mon_spectra_mod
!==============================================================================!
!d                  ==>Full System dimensionality (3*Natoms)
!Natoms             ==>Number of Atoms
!Natoms             ==>Number of Molecules
!atom_type
!Mass
!sqrt_mass
!force
!Lwork              ==>LLAPACK dsygv
!E0_xyz             ==>Initial energy xyz
!potential          ==>Potential name
!d1                 ==>Monomer dimensionality
!d2                 ==>Monomer subspace
!NG                 ==>Number of grid points
!GH_order           ==>Order for the Gauss-Hermite Quadriture
!alpha0             ==>Flat Scaling Parameter for Gaussian Widths
!alpha(NG)          ==>Inverse Gaussian Widths
!x(d,NG)            ==>All grid points coordinates (x^i=x_1,x_2,...x_d)
!V_i                ==>Potential Energy evaluation V(x_i)
!eigenvalues(NG)    ==>Eigenvalues of the matrix
!Smat(NG,NG)        ==>Overlap Matrix
!Hmat(NG,NG)        ==>Hamiltonian Matrix
!==============================================================================!
implicit none
character(len=20)::potential
integer::d,d1,Natoms,Nmol,Lwork0,NG,Lwork,GH_order
double precision::E0_xyz,alpha0,V_i,freq_cutoff,freq_replace
character(len=2),allocatable,dimension(:)::atom_type
double precision,allocatable,dimension(:)::mass,sqrt_mass,q0,forces,omega,alpha
double precision,allocatable,dimension(:)::eigenvalues
double precision,allocatable,dimension(:,:)::Hess_Mat,U,x,Smat,Hmat
!==============================================================================!
double precision,parameter::bohr=0.52917721092
double precision,parameter::autocm=2.194746313D5
double precision,parameter::autokcalmol=627.5096
double precision,parameter::melectron=1822.88839
!==============================================================================!
!                             Read Input Data File
!==============================================================================!
read(*,*) d1
read(*,*) NG
read(*,*) GH_order
read(*,*) alpha0
read(*,*) potential
read(*,*) freq_cutoff
read(*,*) freq_replace
!==============================================================================!
call get_xyz_parameters(d,Natoms,Nmol)
!==============================================================================!
allocate(atom_type(Natoms),mass(Natoms),sqrt_mass(d),q0(d),forces(d))
allocate(Hess_Mat(d,d),omega(d),U(d,d),x(d1,NG),alpha(NG),eigenvalues(NG))
allocate(Smat(NG,NG),Hmat(NG,NG))
Lwork0=max(1,3*d-1)                                          !LLAPACK Suggestion
Lwork=max(1,3*NG-1)                                          !LLAPACK Suggestion
!==============================================================================!
call read_input_geometry(d,Natoms,atom_type,mass,sqrt_mass,q0)
call convert_to_atomic_units(potential,d,Nmol,q0,E0_xyz,forces,freq_cutoff,&
  freq_replace,bohr,autocm,autokcalmol)
call Get_Hessian(potential,d,Nmol,q0,Hess_Mat)
call Mass_Scale_Hessian(d,Hess_Mat,sqrt_mass)
call Frequencies_Scaled_Hess(d,Hess_mat,omega,U,Lwork0)
call replace_small_frequencies(d,omega,freq_cutoff,freq_replace)
call read_grid(d,NG,x)
call generate_alphas(d,NG,alpha0,alpha,x)
!==============================================================================!
call overlap_matrix(d,NG,x,alpha,Smat)
call overlap_eigenvalues(NG,Smat,eigenvalues,Lwork)
call get_hamiltonian(potential,d,Nmol,NG,x,alpha,Smat,GH_order,Hmat,q0,&
  sqrt_mass,U,forces)
call hamiltonian_eigenvalues(NG,Smat,Hmat,eigenvalues,Lwork)
call write_out(potential,d,d1,Nmol,NG,alpha0,GH_order,Natoms,E0_xyz)
!==============================================================================!
end program main
