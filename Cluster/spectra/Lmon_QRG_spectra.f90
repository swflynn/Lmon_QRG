!=================20==================40==================60==================80
!                     Lmon-QRG Vibrational EigenSpectra
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
!d1                 ==>Monomer dimensionality (9 for water by definition)
!d2                 ==>Monomer Subspace dimensionality
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
integer,parameter::d1=9                             !monomer size, :=9 for water
double precision,parameter::bohr=0.52917721092
double precision,parameter::autocm=2.194746313D5
double precision,parameter::autokcalmol=627.5096
!==============================================================================!
character(len=20)::potential
character(len=2),allocatable,dimension(:)::atom_type
integer::d,d2,Natoms,monomer_number,Lwork0,Lwork,NG,GH_order
double precision::E0,alpha0,V
double precision,allocatable,dimension(:)::mass,sqrt_mass,x0,forces,omega,alpha
double precision,allocatable,dimension(:)::eigenvalues
double precision,allocatable,dimension(:,:)::Hess_Mat,U,r,Smat,Hmat
!==============================================================================!
!                             Read Input Data File
!==============================================================================!
read(*,*) d2
read(*,*) potential
read(*,*) NG
read(*,*) GH_order
read(*,*) alpha0
read(*,*) monomer_number
!==============================================================================!
call get_xyz_parameters(d,Natoms)
!==============================================================================!
allocate(atom_type(Natoms),mass(Natoms),sqrt_mass(d),x0(d),forces(d))
allocate(Hess_Mat(d1,d1),omega(d1),U(d1,d1),r(d1,NG),alpha(NG),eigenvalues(NG))
allocate(Smat(NG,NG),Hmat(NG,NG))
Lwork0=max(1,3*d1-1)                                          !LLAPACK Suggestion
Lwork=max(1,3*NG-1)                                          !LLAPACK Suggestion
!==============================================================================!
call read_input_geometry(d,Natoms,atom_type,mass,sqrt_mass,x0)
call convert_to_atomic_units(potential,d,x0,E0,forces,bohr,autocm,&
  autokcalmol)
call Get_Hessian(potential,d,d1,x0,Hess_Mat)
call Mass_Scale_Hessian(d,d1,Hess_Mat,sqrt_mass)
call Frequencies_Scaled_Hess(d1,Hess_mat,omega,U,Lwork0)
call read_grid(d1,d2,NG,r)
call generate_alphas(d1,NG,alpha0,alpha,r)
!==============================================================================!
call overlap_matrix(d1,d2,NG,r,alpha,Smat)
call overlap_eigenvalues(NG,Smat,eigenvalues,Lwork)
call get_hamiltonian(potential,d,d1,d2,NG,r,alpha,Smat,GH_order,Hmat,x0,&
  sqrt_mass,U)
call hamiltonian_eigenvalues(NG,Smat,Hmat,eigenvalues,Lwork)
call write_out(potential,d,d1,d2,Natoms,x0,NG,alpha0,GH_order,E0)
!==============================================================================!
end program main
