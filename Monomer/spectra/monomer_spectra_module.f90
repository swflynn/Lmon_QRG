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
subroutine get_xyz_parameters(d,Natoms,Nmol)
implicit none
integer::d,Natoms,Nmol
open(16,File='tip_geometry.xyz')
read(16,*) Natoms
d=3*Natoms
Nmol=Natoms/3
close(16)
end subroutine get_xyz_parameters
!==============================================================================!
subroutine read_input_geometry(d,Natoms,atom_type,mass,sqrt_mass,q)
implicit none
integer::d,Natoms,i
character(len=2)::atom_type(Natoms)
double precision::mass(Natoms),sqrt_mass(d),q(d)
open(16,File='tip_geometry.xyz')
read(16,*)
read(16,*)
do i=1,Natoms
  read(16,*) atom_type(i),q(3*i-2:3*i)
  mass(i)=Atom_mass(atom_type(i))
  sqrt_mass(3*i-2:3*i)=sqrt(mass(i))
enddo
close(16)
end subroutine read_input_geometry
!==============================================================================!
function Atom_Mass(atom)
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
subroutine convert_to_atomic_units(potential,d,Nmol,q0,E0,forces,freq_cutoff,&
  freq_replace,bohr,autocm,autokcalmol)
!==============================================================================!
!Assumes input geometry is in Angstroms
!==============================================================================!
implicit none
character(len=20)::potential
integer::d,Nmol
double precision::q0(d),E0,forces(d),freq_cutoff,freq_replace,bohr,autocm
double precision::autokcalmol
q0=q0/bohr
freq_cutoff=freq_cutoff/autocm
freq_replace=freq_replace/autocm
call water_potential(potential,Nmol,q0,E0,forces)
write(*,*) 'E0 ==> ', E0*autocm, 'cm-1', E0*autokcalmol, 'kcal/mol'
end subroutine convert_to_atomic_units
!==============================================================================!
subroutine water_potential(potential,Nmol,x,V,forces)
!==============================================================================!
!Determines which Water Potential Energy Function to call
!==============================================================================!
!potential        ==>potential name
!Nmol             ==>Number of water molecules
!x_i(9*Nmol)      ==>coordinates
!V_i              ==>Potential Energy evaluation V(x_i)
!forces(9*Nmol)   ==>Forces computed by external water potential
!==============================================================================!
use iso_c_binding
use TIP4P_module
implicit none
character(len=20)::potential
integer::d,Nmol
double precision::x(9*Nmol),forces(9*Nmol),V
!==============================================================================!
if(potential=='tip4p') then
  call TIP4P(Nmol,x,V,forces)
else
  stop 'Cannot Identify Potential, Check "potentials" Subroutine'
endif
end subroutine water_potential
!==============================================================================!
subroutine Get_Hessian(potential,d,Nmol,q,Hess_Mat)
!==============================================================================!
!Numerically evaluate the Hessian
!==============================================================================!
!s          ==>Perturbation Parameter
!Hess_Mat   ==>(d1,d1): Symmetrized Mass-Scaled Hessian
!q(d)       ==>coordinates
!==============================================================================!
implicit none
integer::d,Nmol,i,j
double precision::Hess_Mat(d,d),q(d),r(d),force(d),force0(d),E
double precision,parameter::s=1d-6
character(len=20)::potential
r=q
call water_potential(potential,Nmol,r,E,force0)
do i=1,d
  r(i)=q(i)+s
  call water_potential(potential,Nmol,r,E,force)
  r(i)=q(i)
  do j=1,d
    Hess_Mat(i,j)=(force0(j)-force(j))/s
  enddo
enddo
end subroutine Get_Hessian
!==============================================================================!
subroutine Mass_Scale_Hessian(d,Hess_Mat,sqrt_mass)
!==============================================================================!
!                   Symmetrize and Mass-Scale the Hessian
!==============================================================================!
!s          ==>Perturbation Parameter
!Hess_Mat   ==>(d1,d1): Symmetrized Mass-Scaled Hessian
!x          ==>(d): XYZ coordinates
!==============================================================================!
implicit none
integer::d,i,j
double precision::Hess_Mat(d,d),sqrt_mass(d)
!==============================================================================!
do i=1,d
  do j=1,i
    if(i.ne.j) Hess_Mat(i,j)=(Hess_Mat(i,j)+Hess_Mat(j,i))/2
    Hess_Mat(i,j)=Hess_Mat(i,j)/(sqrt_mass(i)*sqrt_mass(j))
    if(i.ne.j) Hess_Mat(j,i)=Hess_Mat(i,j)
    enddo
enddo
end subroutine Mass_Scale_Hessian
!==============================================================================!
subroutine Frequencies_Scaled_Hess(d,Hess_mat,omega,U,Lwork0)
!==============================================================================!
!Compute Eigenvalues and Eigenvectors of Hessian
!Uses the LLAPACK real symmetric eigen-solver (dsygev)
!Hess_mat   ==>(d1,d1); Hessian Matrix
!omega      ==>(d1); Hessian Eigenvalues
!U          ==>(d1,d1); Hessian Eigenvectors
!     LLAPACK(dsyev):
!v          ==>Compute both Eigenvalues and Eigenvectors
!u          ==>Use Upper-Triangle of matrix
!==============================================================================!
implicit none
integer::d,i,info,Lwork0
double precision::Hess_mat(d,d),omega(d),U(d,d),work(max(1,lwork0)),temp(d)
double precision::temp2(d,d)
U=Hess_mat
call dsyev('v','u',d,U,d,omega,work,Lwork0,info)
write(*,*) 'Frequencies from the Mass-Scaled Hessian:'
do i=d,1,-1
  omega(i)=sign(sqrt(abs(omega(i))),omega(i))
  write(*,*) omega(i), 'normalized = 1?', sum(U(:,i)**2)
enddo
!==============================================================================!
!Subspace Needs Largest Eigenvalues: llapack outputs small to large ==>re-order
!==============================================================================!
temp=omega
temp2=U
do i=1,d
  omega(i)=temp(d+1-i)
  U(:,i)=temp2(:,d+1-i)
enddo
end subroutine Frequencies_Scaled_Hess
!==============================================================================!
subroutine replace_small_frequencies(d,omega,freq_cutoff,freq_replace)
!monomer is special case need to replace small frequencies
!==============================================================================!
implicit none
integer::d,k
double precision::omega(d),freq_cutoff,freq_replace
do k=1,d
  if(omega(k)<freq_cutoff)then
    omega(k)=freq_replace
  endif
enddo
write(*,*) 'omega replaced ==>', omega
end subroutine replace_small_frequencies
!==============================================================================!
subroutine read_grid(d,NG,x)
!==============================================================================!
!Reads in grid points, each line contains coordinates for a d-dimesional point
!==============================================================================!
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!NG                 ==>Number of grid points
!x(d,NG)            ==>All grid points coordinates (x^i=x_1,x_2,...x_d)
!==============================================================================!
implicit none
integer::d,NG,i
double precision::x(d,NG)
open(17,File='grid.dat')
do i=1,NG
  read(17,*) x(:,i)
enddo
close(17)
end subroutine read_grid
!==============================================================================!
subroutine generate_alphas(d,NG,alpha0,alpha,x)
!==============================================================================!
!Symmetric gaussians (alpha_i is the same for each dimension of the gaussian)
!Uses the nearest neighbor to determine alpha
!==============================================================================!
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!NG                 ==>Number of grid points
!alpha0             ==>Flat Scaling Parameter for Gaussian Widths
!alpha(NG)          ==>Inverse Gaussian Widths
!x(d,NG)            ==>All grid points coordinates (x^i=x_1,x_2,...x_d)
!==============================================================================!
implicit none
integer::d,NG,i,j
double precision::alpha0,alpha(NG),x(d,NG),r2
!==============================================================================!
do i=1,NG
  alpha(i)=1d20                          !large distance for initial placeholder
  do j=1,NG
    if(j.ne.i) then
      r2=sum((x(:,i)-x(:,j))**2)                    !distance between gridpoints
      if(r2<alpha(i)) alpha(i)=r2
    endif
  enddo
  alpha(i)=alpha0/alpha(i)
enddo
call write_alphas(NG,alpha)
end subroutine generate_alphas
!==============================================================================!
subroutine write_alphas(NG,alpha)
!==============================================================================!
!Write alphas to file
!==============================================================================!
!NG                 ==>Number of grid points
!alpha(NG)          ==>Inverse Gaussian Widths
!==============================================================================!
integer::NG,i
double precision::alpha(NG)
open(unit=18,file='alphas.dat')
do i=1,NG
  write(18,*) alpha(i)
enddo
close(18)
endsubroutine write_alphas
!==============================================================================!
subroutine overlap_elements(d,x_i,x_j,alpha_i,alpha_j,S_ij)
!==============================================================================!
!Compute i-jth element of the overlap matrix
!==============================================================================!
!d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x_i(d)           ==>i-th grid points coordinates (x^i=x_1,x_2,...x_d)
!alpha_i          ==>i-th grid points gaussian width parameters
!Sij              ==>i-j element of the overlap matrix
!==============================================================================!
integer::d
double precision::alpha_i,alpha_j,S_ij,x_i(d),x_j(d),aij,r2
aij=alpha_i*alpha_j/(alpha_i+alpha_j)
r2=sum((x_i(:)-x_j(:))**2)
S_ij=(2*sqrt(alpha_i*alpha_j)/(alpha_i+alpha_j))**(0.5*d)*exp(-aij*r2)
end subroutine overlap_elements
!==============================================================================!
subroutine overlap_matrix(d,NG,x,alpha,Smat)
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
integer::d,NG,i,j
double precision::alpha(NG),Smat(NG,NG),x(d,NG)
do i=1,NG
  do j=i,NG
    call overlap_elements(d,x(:,i),x(:,j),alpha(i),alpha(j),Smat(i,j))
    Smat(j,i)=Smat(i,j)
  enddo
enddo
end subroutine overlap_matrix
!==============================================================================!
subroutine overlap_eigenvalues(NG,Smat,eigenvalues,Lwork)
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
double precision::eigenvalues(NG),Smat(NG,NG),work(max(1,Lwork))
call dsyev('v','u',NG,Smat,NG,eigenvalues,work,Lwork,info)              !LLAPACK
write(*,*) 'Info (Overlap Matrix) ==> ', info
write(*,*) 'RCN ==> ', eigenvalues(1)/eigenvalues(NG)
open(unit=19,file='overlap.dat')
do i=1,NG
  write(19,*) eigenvalues(i)
enddo
close(19)
end subroutine overlap_eigenvalues
!==============================================================================!
subroutine kinetic_elements(d,x_i,x_j,alpha_i,alpha_j,T_ij)
!==============================================================================!
!Compute i-jth element of the kinetic matrix
!==============================================================================!
!d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x_i(d)           ==>i-th grid points coordinates (x^i=x_1,x_2,...x_d)
!alpha_i          ==>i-th grid points gaussian width parameter
!T_ij             ==>i-j element of the kinetic matrix
!==============================================================================!
implicit none
integer::d
double precision::x_i(d),x_j(d),alpha_i,alpha_j,aij,r2,T_ij
aij=alpha_i*alpha_j/(alpha_i+alpha_j)
r2=sum((x_i(:)-x_j(:))**2)
T_ij=aij*(d-2.*aij*r2)
end subroutine kinetic_elements
!==============================================================================!
subroutine get_hamiltonian(potential,d,Nmol,NG,x,alpha,Smat,GH_order,Hmat,q0,&
  sqrt_mass,U,forces)
!==============================================================================!
!Compute the Hamiltonian Matrix
!Use Gauss-Hermite quadrature for the potential
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
integer::d,Nmol,NG,GH_order,i,j,k,ll
double precision,parameter::pi=4.*atan(1d0)
double precision::x(d,NG),alpha(NG),Smat(NG,NG),Hmat(NG,NG),z(GH_order)
double precision::w(GH_order),x_ij(d),V_ij,l(d),rr(d),r2,q0(d),q1(d),forces(d)
double precision::sqrt_mass(d),U(d,d)
call cgqf(GH_order,6,0d0,0d0,0d0,1d0,z,w)                             !GHQ
w=w/sqrt(pi)
!==============================================================================!
do i=1,NG
  do j=i,NG
    call overlap_elements(d,x(:,i),x(:,j),alpha(i),alpha(j),Smat(i,j))
    call kinetic_elements(d,x(:,i),x(:,j),alpha(i),alpha(j),Hmat(i,j))
!==============================================================================!
!                    Compute Potential Energy with Quadriture
!==============================================================================!
    x_ij(:)=(alpha(i)*x(:,i)+alpha(j)*x(:,j))/(alpha(i)+alpha(j))
    V_ij=0d0
    l(:)=1
    do ll=1,GH_order**d
      do k=1,d
        rr(k)=z(l(k))
      enddo
      rr=x_ij+rr/sqrt(alpha(i)+alpha(j))
!==============================================================================!
!Scale to cartesian space
!==============================================================================!
      q1=q0
      q1=q1+matmul(U,rr)/sqrt_mass
      call water_potential(potential,Nmol,q1,r2,forces)
      do k=1,d
        r2=r2*w(l(k))
      enddo
      V_ij=V_ij+r2
      do k=1,d
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
subroutine hamiltonian_eigenvalues(NG,Smat,Hmat,eigenvalues,Lwork)
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
itype=1
call dsygv(itype,'n','u',NG,Hmat,NG,Smat,NG,eigenvalues,work,Lwork,info)
open(unit=20,file='eigenvalues.dat')
do i=1,NG
  write(20,*) eigenvalues(i)
enddo
close(20)
end subroutine hamiltonian_eigenvalues
!==============================================================================!
subroutine write_out(potential,d,d1,Nmol,NG,alpha0,GH_order,Natoms,E0_xyz)
!==============================================================================!
!Simulation details
!==============================================================================!
!alpha_method       ==>Alpha generation method (constant, nearest, distribution)
!potential          ==>Potential name
!d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!NG               ==>Number of grid points
!alpha0           ==>Flat scaling parameter for optimizing widths
!GH_order         ==>Gauss-Hermite quadrature order
!==============================================================================!
implicit none
integer::d,d1,Nmol,NG,GH_order,Natoms
double precision::alpha0,E0_xyz
character(len=20)::potential
open(unit=99,file='out')
write(99,*) 'potential ==> ', potential
write(99,*) 'dimensionality ==> ', d
write(99,*) 'subspace dimensionality ==> ', d1
write(99,*) 'Number of Molecules==> ', Nmol
write(99,*) 'Number of Gaussians==> ', NG
write(99,*) 'alpha0==>', alpha0
write(99,*) 'GH Order==>', GH_order
write(99,*) 'Natoms ==>', Natoms
write(99,*) 'E0_xyz==>', E0_xyz
close(99)
end subroutine write_out
!==============================================================================!
end module mon_spectra_mod
