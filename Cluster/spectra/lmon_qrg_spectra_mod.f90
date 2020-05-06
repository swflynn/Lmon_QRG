!=================20==================40==================60==================80
!              nD Vibrational EigenSpectra monomer implementation
!==============================================================================!
!begin developing Lmon transformations start by testing a single monomer with
!the tip4p potential and compare against our previour results
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
!   21 April 2020
!       Author:
!   Shane Flynn
!==============================================================================!
module mon_spectra_mod
implicit none
!==============================================================================!
contains
!==============================================================================!
subroutine get_xyz_parameters(d,Natoms)
!==============================================================================!
!Read xyz file to determine d, Natoms for allocations.
!==============================================================================!
!Natoms             ==>Number of atoms in the system
!d                  ==>Total System Dimensionality  (d:=3*Natoms)
!==============================================================================!
implicit none
integer::d,Natoms
open(16,File='cage_tip4p.xyz')
read(16,*) Natoms
d=3*Natoms
close(16)
end subroutine get_xyz_parameters
!==============================================================================!
subroutine read_input_geometry(d,Natoms,atom_type,mass,sqrt_mass,x0)
!==============================================================================!
!Read in initial minimized geometry
!==============================================================================!
!d                  ==>Total System Dimensionality  (d:=3*Natoms)
!Natoms             ==>Number of atoms in the system
!mass(Natoms)       ==>Mass
!sqrt_mass(d)       ==>Square Root Mass
!q(d)               ==>Coordinates (entire system)
!atom_type(Natoms)  ==>Element abbreviation from xyz
!==============================================================================!
implicit none
integer::d,Natoms,i
character(len=2)::atom_type(Natoms)
double precision::mass(Natoms),sqrt_mass(d),x0(d)
open(17,File='cage_tip4p.xyz')
read(17,*)
read(17,*)
do i=1,Natoms
  read(17,*) atom_type(i),x0(3*i-2:3*i)                 !input is xyz therefore 3
  mass(i)=Atom_mass(atom_type(i))
  sqrt_mass(3*i-2:3*i)=sqrt(mass(i))
enddo
close(17)
end subroutine read_input_geometry
!==============================================================================!
function Atom_Mass(atom)
!==============================================================================!
!Compute mass for each atom
!==============================================================================!
implicit none
double precision::Atom_Mass
character(len=2)::atom
double precision,parameter::melectron=1822.88839
double precision,parameter::Hmass=1.00782503223*melectron
double precision,parameter::Omass=15.99491461957*melectron
if(atom=='H'.or.atom=='h')then
  Atom_mass=Hmass
elseif(atom=='O'.or.atom=='o')then
  Atom_mass=Omass
else
  write(*,*) 'atom ', atom, ' is not recognized'
  stop 'Check Atom_Mass Function'
endif
end function Atom_Mass
!==============================================================================!
subroutine convert_to_atomic_units(potential,d,x0,E0,forces,bohr,autocm,&
  autokcalmol)
!==============================================================================!
!Use atomic units for calculations (Assumes input geometry is in Angstroms)
!==============================================================================!
!potential        ==>potential name
!d                  ==>Total System Dimensionality  (d:=3*Natoms)
!q0(d)              ==>Initial Coordinates (entire system)
!E0                 ==>Potential Energy of q0
!forces(d)          ==>Forces from potential
!freq_cutoff        ==>Replace small frequencies below this cutoff
!freq_replace       ==>Replace value for small frequencies below the cutoff
!==============================================================================!
implicit none
character(len=20)::potential
integer::d
double precision::x0(d),E0,forces(d),bohr,autocm,autokcalmol
x0=x0/bohr
call water_potential(potential,d,x0,E0,forces)
end subroutine convert_to_atomic_units
!==============================================================================!
subroutine water_potential(potential,d,x,V,forces)
!==============================================================================!
!Call the water potential
!To call the water potentials you need to pass in the number of water atoms in
!the system. d=3*Natoms, Nmol=Natoms/3 ==> d/9=Nmol
!==============================================================================!
!potential          ==>potential name
!d                  ==Total system dimensionality
!x(d)               ==>coordinates
!V                  ==>Potential Energy evaluation V(x)
!forces(d)          ==>Forces from potential
!==============================================================================!
use iso_c_binding
use TIP4P_module
implicit none
character(len=20)::potential
integer::d
double precision::x(d),forces(d),V
!==============================================================================!
if(potential=='tip4p'.or.potential=="TIP4P") then
  call TIP4P(d/9,x,V,forces)
else
  stop 'Cannot Identify Potential, Check "water_potential" Subroutine'
endif
end subroutine water_potential
!==============================================================================!
subroutine Get_Hessian(potential,d,d1,x0,Hess_Mat)
!==============================================================================!
!Numerically evaluate the Hessian
!==============================================================================!
!potential          ==>potential name
!d                  ==>Total System Dimensionality  (d:=3*Natoms)
!x0(d)              ==>Initial Coordinates (entire system)
!force(d)           ==>Forces from potential
!E0                 ==>Potential Energy of q0
!s                  ==>Perturbation Parameter
!Hess_Mat(d1,d1)    ==>Numerical Hessian
!==============================================================================!
implicit none
integer::d,d1,i,j
double precision::Hess_Mat(d1,d1),x0(d),x1(d),force0(d),force1(d),E0
double precision,parameter::s=1d-6
character(len=20)::potential
x1=x0
call water_potential(potential,d,x1,E0,force0)
do i=1,d1
  x1(i)=x0(i)+s
  call water_potential(potential,d,x1,E0,force1)
  x1(i)=x0(i)
  do j=1,d1
    Hess_Mat(i,j)=(force0(j)-force1(j))/s
  enddo
enddo
end subroutine Get_Hessian
!==============================================================================!
subroutine Mass_Scale_Hessian(d,d1,Hess_Mat,sqrt_mass)
!==============================================================================!
!Symmetrize and Mass-Scale the Hessian
!==============================================================================!
!d                  ==>Total System Dimensionality  (d:=3*Natoms)
!Hess_Mat(d,d)      ==>Numerical Hessian
!sqrt_mass(d)       ==>Square Root Mass
!==============================================================================!
implicit none
integer::d,d1,i,j
double precision::Hess_Mat(d1,d1),sqrt_mass(d)
!==============================================================================!
do i=1,d1
  do j=1,i
    if(i.ne.j) Hess_Mat(i,j)=(Hess_Mat(i,j)+Hess_Mat(j,i))/2
    Hess_Mat(i,j)=Hess_Mat(i,j)/(sqrt_mass(i)*sqrt_mass(j))
    if(i.ne.j) Hess_Mat(j,i)=Hess_Mat(i,j)
    enddo
enddo
end subroutine Mass_Scale_Hessian
!==============================================================================!
subroutine Frequencies_Scaled_Hess(d1,Hess_mat,omega,U,Lwork)
!==============================================================================!
!Compute Eigenvalues and Eigenvectors for the mass-scaled hessian
!Uses the LLAPACK real symmetric eigen-solver (dsygev)
!==============================================================================!
!d                  ==>Total System Dimensionality  (d:=3*Natoms)
!Hess_Mat(d,d)      ==>Numerical Hessian
!omega(d)           ==>Eigenvalues
!U(d,d)             ==>Eigenvectors
!     LLAPACK(dsyev):
!v                  ==>Compute both Eigenvalues and Eigenvectors
!u                  ==>Use Upper-Triangle of matrix
!Lwork              ==>Allocation size
!==============================================================================!
implicit none
integer::d1,i,info,Lwork
double precision::Hess_mat(d1,d1),omega(d1),U(d1,d1),work(max(1,Lwork))
double precision::temp(d1),temp2(d1,d1)
U=Hess_mat
call dsyev('v','u',d1,U,d1,omega,work,Lwork,info)
open(18,File='freq_hess.dat')
write(18,*) 'Frequencies from the Mass-Scaled Hessian:'
do i=d1,1,-1
  omega(i)=sign(sqrt(abs(omega(i))),omega(i))
  write(18,*) omega(i), 'normalized = 1?', sum(U(:,i)**2)
enddo
close(18)
!==============================================================================!
!Subspace Needs Largest Eigenvalues: llapack outputs small to large ==>re-order
!==============================================================================!
temp=omega
temp2=U
do i=1,d1
  omega(i)=temp(d1+1-i)
  U(:,i)=temp2(:,d1+1-i)
enddo
end subroutine Frequencies_Scaled_Hess
!==============================================================================!
subroutine read_grid(d1,d2,NG,r)
!==============================================================================!
!Reads in grid points, each line contains coordinates for a d2-dimesional point
!need to check and see if this works, currently not scaling by omega
!==============================================================================!
!d2                 ==>Monomer Subspace dimensionality (x^i=x_1,x_2,...x_d2)
!NG                 ==>Number of grid points
!r(d2,NG)           ==>All grid points coordinates (x^i=x_1,x_2,...x_d)
!==============================================================================!
implicit none
integer::d1,d2,NG,i
double precision::r(d1,NG)
r=0                                                           !set (d2+1:d1) = 0
open(19,File='grid.dat')
do i=1,NG
  read(19,*) r(1:d2,i)          !grid is generated in d2 subspace, set rest to 0
enddo
close(19)
end subroutine read_grid
!==============================================================================!
subroutine generate_alphas(d1,NG,alpha0,alpha,r)
!==============================================================================!
!Symmetric gaussians (alpha_i is the same for each dimension of the gaussian)
!Uses the nearest neighbor to determine alpha
!==============================================================================!
!d1                 ==>Monomer dimensionality (x^i=x_1,x_2,...x_d1)
!NG                 ==>Number of grid points
!alpha0             ==>Flat Scaling Parameter for Gaussian Widths
!alpha(NG)          ==>Inverse Gaussian Widths
!r(d1,NG)           ==>All grid points coordinates (x^i=x_1,x_2,...x_d1)
!==============================================================================!
implicit none
integer::d1,NG,i,j
double precision::alpha0,alpha(NG),r(d1,NG),r2
!==============================================================================!
do i=1,NG
  alpha(i)=1d20                          !large distance for initial placeholder
  do j=1,NG
    if(j.ne.i) then
      r2=sum((r(:,i)-r(:,j))**2)                    !distance between gridpoints
      if(r2<alpha(i)) alpha(i)=r2
    endif
  enddo
  alpha(i)=alpha0/alpha(i)
enddo
open(unit=20,file='alphas.dat')
do i=1,NG
  write(20,*) alpha(i)
enddo
close(20)
end subroutine generate_alphas
!==============================================================================!
subroutine overlap_elements(d1,d2,r_i,r_j,alpha_i,alpha_j,S_ij)
!==============================================================================!
!Compute i-jth element of the overlap matrix
!scale overlap by the subspace dimensionality
!==============================================================================!
!n                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_n)
!x_i(n)           ==>i-th grid points coordinates (x^i=x_1,x_2,...x_n)
!alpha_i          ==>i-th grid points gaussian width parameters
!Sij              ==>i-j element of the overlap matrix
!==============================================================================!
integer::d1,d2
double precision::alpha_i,alpha_j,S_ij,r_i(d1),r_j(d1),aij,r2
aij=alpha_i*alpha_j/(alpha_i+alpha_j)
r2=sum((r_i(:)-r_j(:))**2)
S_ij=(2*sqrt(alpha_i*alpha_j)/(alpha_i+alpha_j))**(0.5*d2)*exp(-aij*r2)
end subroutine overlap_elements
!==============================================================================!
subroutine overlap_matrix(d1,d2,NG,r,alpha,Smat)
!==============================================================================!
!Compute the overlap matrix (symmetric, positive-definite)
!Used to check the stability of the basis; see overlap_eigenvalues subroutine
!==============================================================================!
!d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!NG               ==>Number of grid points
!x(d,NG)          ==>All grid points coordinates (x={x^i})
!alpha(NG)        ==>Gaussian width parameter
!Smat(NG,NG)      ==>Overlap Matrix
!==============================================================================!
integer::d1,d2,NG,i,j
double precision::alpha(NG),Smat(NG,NG),r(d1,NG)
do i=1,NG
  do j=i,NG
    call overlap_elements(d1,d2,r(:,i),r(:,j),alpha(i),alpha(j),Smat(i,j))
    Smat(j,i)=Smat(i,j)
  enddo
enddo
end subroutine overlap_matrix
!==============================================================================!
subroutine overlap_eigenvalues(NG,Smat,eigenvalues,Lwork,RCN)
!==============================================================================!
!Compute the Eigenvalues for the overlap matrix to check numerical stability
!Must be positive definite; check the Recriprocal Condition Number
!Lwork, work, and info defined according to LLAPACK suggestions
!==============================================================================!
!NG               ==>Number of grid points
!Smat(NG,NG)      ==>Overlap Matrix
!eigenvalues(NG)  ==>Eigenvalues of the matrix
!==============================================================================!
implicit none
integer::i,NG,info,Lwork
double precision::eigenvalues(NG),Smat(NG,NG),work(max(1,Lwork)),RCN
call dsyev('v','u',NG,Smat,NG,eigenvalues,work,Lwork,info)              !LLAPACK
RCN=eigenvalues(1)/eigenvalues(NG)
open(21,file='overlap.dat')
do i=1,NG
  write(21,*) eigenvalues(i)
enddo
close(21)
end subroutine overlap_eigenvalues
!==============================================================================!
subroutine kinetic_elements(d1,d2,r_i,r_j,alpha_i,alpha_j,T_ij)
!==============================================================================!
!Compute i-jth element of the kinetic matrix
!scale according to subspace dimensionality
!==============================================================================!
!d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x_i(d)           ==>i-th grid points coordinates (x^i=x_1,x_2,...x_d)
!alpha_i          ==>i-th grid points gaussian width parameter
!T_ij             ==>i-j element of the kinetic matrix
!==============================================================================!
implicit none
integer::d1,d2
double precision::r_i(d1),r_j(d1),alpha_i,alpha_j,aij,r2,T_ij
aij=alpha_i*alpha_j/(alpha_i+alpha_j)
r2=sum((r_i(:)-r_j(:))**2)
T_ij=aij*(d2-2.*aij*r2)
end subroutine kinetic_elements
!==============================================================================!
subroutine scale_evaluate_water_potential(potential,d,d1,x0,rr,E,U,sqrt_mass)
!==============================================================================!
!Evaluate potential in d-dimensionality. Lmon moves 1 monomer of the cluster
!We only have dynamics over d2 subspace of the d1 monomer
!==============================================================================!
!x0(d)        ==>initial cluster configuration
!E            ==>Potential Energy Evalation
!U            ==>Normal mode eigenvectors
!rr is the scaled r_ij coordinate to evaluate
!==============================================================================!
implicit none
character(len=20)::potential
integer::d,d1
double precision::x0(d),x1(d),rr(d1),U(d1,d1),forces(d),sqrt_mass(d),E
x1=x0
x1(1:d1)=x1(1:d1)+matmul(U,rr)/sqrt_mass(1:d1)
call water_potential(potential,d,x1,E,forces)
end subroutine scale_evaluate_water_potential
!==============================================================================!
subroutine get_hamiltonian(potential,d,d1,d2,NG,r,alpha,Smat,GH_order,Hmat,x0,&
  sqrt_mass,U)
!==============================================================================!
!Compute the Hamiltonian Matrix
!Use Gauss-Hermite quadrature for the potential
!Potential is evaluated in Cartesian space
!==============================================================================!
!potential        ==>Potential name
!d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!NG               ==>Number of grid points
!x(d,NG)          ==>All grid points coordinates (x={x^i})
!alpha(NG)        ==>Gaussian width parameter
!Smat(NG,NG)      ==>Overlap Matrix
!GH_order         ==>Gauss-Hermite quadrature order
!Hmat(NG.NG)      ==>Hamiltonian Matrix
!z(GH_order)      ==>Quadrature points
!w(GH_order)      ==>Quadrature points weights
!==============================================================================!
implicit none
character(len=20)::potential
integer::d,d1,d2,NG,GH_order,i,j,k,ll
double precision,parameter::pi=4.*atan(1d0)
double precision::r(d1,NG),alpha(NG),Smat(NG,NG),Hmat(NG,NG),z(GH_order)
double precision::w(GH_order),r_ij(d1),V_ij,l(d1),rr(d1),E,x0(d)
double precision::sqrt_mass(d),U(d1,d1)
call cgqf(GH_order,6,0d0,0d0,0d0,1d0,z,w)                             !GHQ
w=w/sqrt(pi)
!==============================================================================!
do i=1,NG
  do j=i,NG
    call overlap_elements(d1,d2,r(:,i),r(:,j),alpha(i),alpha(j),Smat(i,j))
    call kinetic_elements(d1,d2,r(:,i),r(:,j),alpha(i),alpha(j),Hmat(i,j))
!==============================================================================!
!              Compute Potential Energy with Quadriture
!==============================================================================!
    V_ij=0d0
    r_ij(:)=(alpha(i)*r(:,i)+alpha(j)*r(:,j))/(alpha(i)+alpha(j))
    l(:)=1
    do ll=1,GH_order**d
      do k=1,d1
        rr(k)=z(l(k))
      enddo
      rr=r_ij+rr/sqrt(alpha(i)+alpha(j))
      call scale_evaluate_water_potential(potential,d,d1,x0,rr,E,U,sqrt_mass)
      do k=1,d1
        E=E*w(l(k))
      enddo
      V_ij=V_ij+E
      do k=1,d1
        l(k)=mod(l(k),float(GH_order))+1
        if(l(k).ne.1) exit
      enddo
    enddo
    Hmat(i,j)=(Hmat(i,j)+V_ij)*Smat(i,j)
    Hmat(j,i)=Hmat(i,j)
  enddo
enddo
end subroutine get_hamiltonian
!==============================================================================!
subroutine hamiltonian_eigenvalues(NG,Smat,Hmat,eigenvalues,Lwork,autocm,alpha0)
!==============================================================================!
!Compute the Eigenvalues for the Hamiltonian matrix
!Lwork, work, and info defined according to LLAPACK suggestions
!==============================================================================!
!NG               ==>Number of grid points
!Smat(NG,NG)      ==>Overlap Matrix
!Hmat(NG.NG)      ==>Hamiltonian Matrix
!eigenvalues(NG)  ==>Matrix eigenvalues
!==============================================================================!
implicit none
integer::NG,info,Lwork,itype,i
double precision::Smat(NG,NG),Hmat(NG,NG),eigenvalues(NG),work(max(1,Lwork))
double precision::autocm,alpha0
itype=1
call dsygv(itype,'n','u',NG,Hmat,NG,Smat,NG,eigenvalues,work,Lwork,info)
open(unit=22,file='eigenvalues_atomic.dat')
open(unit=23,file='eigenvalues_cm-1.dat')
open(unit=24,file='fundamentals_cm-1.dat')
open(unit=25,file='alpha_100_atomic.dat')
open(unit=26,file='alpha_100_cm-1.dat')
do i=1,NG
  write(22,*) eigenvalues(i)
  write(23,*) eigenvalues(i)*autocm
  write(25,*) alpha0, eigenvalues(i)
  write(26,*) alpha0, eigenvalues(i)*autocm
enddo
write(24,*) alpha0, (eigenvalues(2)-eigenvalues(1))*autocm
write(24,*) alpha0, (eigenvalues(3)-eigenvalues(1))*autocm
write(24,*) alpha0, (eigenvalues(4)-eigenvalues(1))*autocm
close(22)
close(23)
close(24)
end subroutine hamiltonian_eigenvalues
!==============================================================================!
subroutine write_out(potential,d,d1,d2,Natoms,x0,NG,alpha0,GH_order,E0,RCN)
!==============================================================================!
!Simulation details
!==============================================================================!
!alpha_method       ==>Alpha generation method (constant, nearest, distribution)
!potential          ==>Potential name !d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!NG               ==>Number of grid points
!alpha0           ==>Flat scaling parameter for optimizing widths
!GH_order         ==>Gauss-Hermite quadrature order
!==============================================================================!
implicit none
integer::d,d1,d2,Natoms,NG,GH_order
double precision::alpha0,E0,x0(d),RCN
character(len=20)::potential
open(unit=99,file='out')
write(99,*) 'potential ==> ', potential
write(99,*) 'system dimensionality ==> ', d
write(99,*) 'monomer dimensionality ==> ', d1
write(99,*) 'monomer subspace dimensionality ==> ', d2
write(99,*) 'Number of atoms==> ', Natoms
write(99,*) 'Number of Gaussians==> ', NG
write(99,*) 'alpha0==>', alpha0
write(99,*) 'GH Order==>', GH_order
write(99,*) 'E0==>', E0
write(99,*) 'Initial configuration==> ', x0
write(99,*) 'RCN==>', RCN
close(99)
end subroutine write_out
!==============================================================================!
end module mon_spectra_mod
