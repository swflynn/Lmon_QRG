!=================20==================40==================60==================80
!                         nD Vibrational EigenSpectra
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
!   15 April 2020
!       Author:
!   Shane Flynn
!==============================================================================!
program main
use Lmon_spectra_mod
!==============================================================================!
!d                  ==>Full System dimensionality (3*Natoms)
!Natoms             ==>
!atom_type
!Mass
!sqrt_mass
!force
!lwork              ==>LLAPACK dsygv
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
integer::d,Natoms,Lwork0,NG,Lwork,GH_order
double precision::E0_xyz,alpha0,V_i
character(len=2),allocatable,dimension(:)::atom_type
double precision,allocatable,dimension(:)::mass,sqrt_mass,q0,omega,alpha
double precision,allocatable,dimension(:)::eigenvalues
double precision,allocatable,dimension(:,:)::Hess_Mat,U,x,Smat,Hmat
!==============================================================================!
!                             Read Input Data File
!==============================================================================!
!read(*,*) d1
!read(*,*) d2
read(*,*) NG
read(*,*) GH_order
read(*,*) alpha0
read(*,*) potential
!==============================================================================!
call get_xyz_parameters(d,Natoms)
write(*,*) 'd ==> ',d
!write(*,*) 'Natoms==> ',Natoms
!==============================================================================!
allocate(atom_type(Natoms),mass(Natoms),sqrt_mass(d),q0(d),Hess_Mat(d,d))
allocate(omega(d),U(d,d),x(d,NG),alpha(NG),eigenvalues(NG),Smat(NG,NG))
allocate(Hmat(NG,NG))
Lwork0=max(1,3*d-1)                                          !LLAPACK Suggestion
Lwork=max(1,3*NG-1)                                          !LLAPACK Suggestion
!==============================================================================!
call read_input_geometry(d,Natoms,atom_type,mass,sqrt_mass,q0)
call potentials(potential,d,q0,E0_xyz)
call Get_Hessian(potential,d,q0,Hess_Mat)
call Mass_Scale_Hessian(d,Hess_Mat,sqrt_mass)
call Frequencies_Scaled_Hess(d,Hess_mat,omega,U,Lwork0)
call read_grid(d,NG,x)
call generate_alphas(d,NG,alpha0,alpha,x)
!==============================================================================!
call overlap_matrix(d,NG,x,alpha,Smat)
call overlap_eigenvalues(NG,Smat,eigenvalues,Lwork)
call get_hamiltonian(potential,d,NG,x,alpha,Smat,GH_order,Hmat)
call hamiltonian_eigenvalues(NG,Smat,Hmat,eigenvalues,Lwork)
call write_out(potential,d,NG,alpha0,GH_order,Natoms,E0_xyz)
!==============================================================================!
end program main
