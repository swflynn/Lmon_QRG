!=================20==================40==================60==================80
!                       3D Tilted Harmonic EigenSpectra
!==============================================================================!
!       Modified:
!   13 Feburary 2020
!       Author:
!   Shane Flynn
!==============================================================================!
!Develop/test normal mode coordinates and and subspace calculations using a
!tilted harmonic potential.
!==============================================================================!
module toy_mod
implicit none
!==============================================================================!
!                            Global Variables
!d              ==>total system dimensionality
!d1             ==>monomer dimensionality
!d2             ==>monomer subspace dimensionality
!==============================================================================!
double precision,parameter::Hmass=1d0
!==============================================================================!
integer::d,d1,d2,Natoms
character(len=2),allocatable::atom_type(:)
double precision,allocatable::mass(:),sqrt_mass(:)
!==============================================================================!
contains
!==============================================================================!
subroutine potential(x,energy)
!==============================================================================!
!Tilted harmonic potential for testing normal mode transformations
!V:=(x-y)^2 + (x+y)^2 + (z-x)^2
!==============================================================================!
!x              ==>(d) ith atoms coordinates
!V              ==>evaluate V(x_i)
!==============================================================================!
implicit none
double precision::x(d),energy
energy=(x(1)-x(2))**2+(x(1)+x(2))**2+(x(3)-x(1))**2
end subroutine potential
!==============================================================================!
function Atom_Mass(atom)
!==============================================================================!
implicit none
double precision::Atom_Mass
character(len=2)::atom
if(atom=='H'.or.atom=='h')then
  Atom_mass=Hmass
else
  write(*,*) 'Atom ', atom, ' is not recognized.'
  STOP 'Check Atom_Mass Function'
endif
end function Atom_Mass
!==============================================================================!
subroutine Toy_Force(x,forces)
!==============================================================================!
!Returns the Forces associated with Toy_Potential Subroutine
!Forces are hard-coded based on Toy_Potential
!==============================================================================!
implicit none
integer::i
double precision::x(d),forces(d)
forces(1)=-2.*(3.*x(1)-x(3))
forces(2)=-4.*x(2)
forces(3)=-2.*(x(3)-x(1))
end subroutine Toy_Force
!==============================================================================!
subroutine Toy_Hessian(x,Hess_Mat)
!==============================================================================!
!Numerically evaluate the Hessian using forces from "Toy_Force Subroutine"
!==============================================================================!
!s          ==>Perturbation Parameter
!Hess_Mat   ==>(d1,d1): Symmetrized Mass-Scaled Hessian
!x          ==>(d): XYZ coordinates
!==============================================================================!
implicit none
integer::i,j
double precision::Hess_Mat(d1,d1),x(d),r(d),force(d),force0(d)
double precision,parameter::s=1d-6
r=x
call Toy_Force(r,force0)
do i=1,d1
  r(i)=x(i)+s
  call Toy_Force(r,force)
  r(i)=x(i)
  do j=1,d1
    Hess_Mat(i,j)=(force0(j)-force(j))/s
  enddo
enddo
!==============================================================================!
!                   Symmetrize and Mass-Scale the Hessian
!==============================================================================!
do i=1,d1
  do j=1,i
    if(i.ne.j) Hess_Mat(i,j)=(Hess_Mat(i,j)+Hess_Mat(j,i))/2
    Hess_Mat(i,j)=Hess_Mat(i,j)/(sqrt_mass(i)*sqrt_mass(j))
    if(i.ne.j) Hess_Mat(j,i)=Hess_Mat(i,j)
    enddo
enddo
write(*,*) 'Toy Hessian ==> ', Hess_Mat
end subroutine Toy_Hessian
!==============================================================================!
subroutine Freq_Hess(Hess_mat,omega,U)
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
integer::i,info,lwork,Dimen
double precision::Hess_mat(d,d),omega(d),U(d,d)
double precision,allocatable::work(:)
lwork=max(1,3*d1-1)                              !suggested by LAPACK Developers
allocate(work(max(1,lwork)))                     !suggested by LAPACK Developers
U=Hess_mat
call dsyev('v','u',d1,U,d1,omega,work,lwork,info)
write(*,*) 'Frequencies from the Hessian:'
do i=d1,1,-1
  omega(i)=sign(sqrt(abs(omega(i))),omega(i))
  write(*,*) omega(i), 'normalized = 1?', sum(U(:,i)**2)
enddo
end subroutine Freq_Hess
!==============================================================================!
end module toy_mod
!==============================================================================!
program main
use toy_mod
!==============================================================================!
!grid_in        ==>Filename Containing Gridpoints
!NG             ==>Number of Gaussian Basis Functions (gridpoints)
!GH_order       ==>Number of Points for evaluating the potential (Gauss-Hermite)
!alpha0         ==>Flat Scaling Parameter for Gaussian Widths
!alpha          ==>(d) Gaussian Widths
!x              ==>(d) ith atoms coordinates
!x_ij           ==>i-jth gaussian center (product of Gaussians is a Gaussian)
!Hess           ==>
!Smat           ==>(NG,NG) Overlap Matrix
!Hmat           ==>(NG,NG) Hamiltonian Matrix
!==============================================================================!
implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!character(len=50)::theory_in
!!!integer::GH_order,j,k,ll
!!!double precision::Vij
!!!double precision,parameter::pi=4.*atan(1d0)
!!!double precision,allocatable,dimension(:)::eigenvalues,x_ij,z,w,rr,l
!!!double precision,allocatable,dimension(:)::theory
!!!double precision,allocatable,dimension(:,:)::x,Smat,Hmat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
character(len=50)::coord_in,grid_in
double precision::E0
integer::NG,i,j
logical::unif_grid
double precision::alpha0,aij,x2,rcn
double precision,allocatable,dimension(:)::q0,force,omega,alpha,eigenvalues
double precision,allocatable,dimension(:,:)::Hess,U,x,Smat
!==============================================================================!
!                       LLAPACK dsygv variables                                !
!==============================================================================!
integer::itype,info,lwork
double precision,allocatable::work(:)
!==============================================================================!
!                           Read Input Data File                               !
!==============================================================================!
read(*,*) coord_in
read(*,*) grid_in
read(*,*) d1
read(*,*) d2
read(*,*) NG
read(*,*) unif_grid
read(*,*) alpha0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!read(*,*) GH_order
!!!!read(*,*) alpha0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==============================================================================!
!                         Set Input Water Geometry
!==============================================================================!
open(17,File=coord_in)
read(17,*) Natoms
read(17,*)
d=3*Natoms
!==============================================================================!
allocate(atom_type(Natoms),mass(Natoms),sqrt_mass(d),q0(d),force(d),Hess(d,d))
allocate(omega(d),U(d,d),x(d,NG),alpha(NG),Smat(NG,NG),eigenvalues(NG))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!allocate(x(d2,NG),x_ij(d),rr(d))
!!!!allocate(Hmat(NG,NG),z(GH_order),w(GH_order),l(d),theory(NG))
!!!!write(*,*) 'Test 0; Successfully Read Input File'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==============================================================================!
!                         Input Configuration Energy
!==============================================================================!
do i=1,Natoms
  read(17,*) atom_type(i),q0(3*i-2:3*i)                  !xyz spacial dimensions
  mass(i)=Atom_mass(atom_type(i))
  write(*,*) 'mass', mass(i)
  sqrt_mass(3*i-2:3*i)=sqrt(mass(i))
  write(*,*) 'sqrt mass', sqrt_mass(i)
enddo
close(17)
write(*,*) 'd system dimensionality ==> ', d
write(*,*) 'd1 monomer dimensionality ==> ', d1
write(*,*) 'd2 monomer subspace dimensionality ==> ', d2
write(*,*) 'q0 ==>', q0
call potential(q0,E0)
write(*,*) 'E0 ==> ', E0
!==============================================================================!
! 			                Compute Hessian and Frequencies
!==============================================================================!
call Toy_Hessian(q0,Hess)
call Freq_Hess(Hess,omega,U)
write(*,*) 'Test 2; Successfully Computed Hessian'
!==============================================================================!
!                         Read GridPoints x(d2,NG)
!==============================================================================!
open(18,File=grid_in)
do i=1,NG
  read(18,*) x(:,i)
enddo
close(18)
write(*,*) 'Test 3; Successfully Read Grid Points'
!==============================================================================!
!                       Generate Gaussian Widths
!==============================================================================!
if(unif_grid.EQV..TRUE.)then
  alpha=alpha0
  write(*,*) 'Test 3; Uniform Grid, Alpha:=Constant', alpha0
else
  do i=1,NG
    alpha(i)=1d20                            !large distance for placeholder
    do j=1,NG
      if(j.ne.i) then
        x2=sum((x(:,i)-x(:,j))**2)          !distance between gridpoints
        if(x2<alpha(i)) alpha(i)=x2
      endif
    enddo
    alpha(i)=alpha0/alpha(i)
  enddo
  write(*,*) 'Test 3; Successfully Generated Alphas according to P(x)'
endif
!==============================================================================!
!                              Write Alphas to File
!==============================================================================!
open(unit=19,file='alphas.dat')
do i=1,NG
  write(19,*) alpha(i)
enddo
close(19)
!==============================================================================!
!                               Overlap Matrix (S)
!==============================================================================!
do i=1,NG
  do j=i,NG
    aij=alpha(i)*alpha(j)/(alpha(i)+alpha(j))
    x2=sum((x(:,i)-x(:,j))**2)
    Smat(i,j)=(2*sqrt(alpha(i)*alpha(j))/(alpha(i)+alpha(j)))**(0.5*d)&
    *exp(-aij*x2)
    Smat(j,i)=Smat(i,j)
  enddo
enddo
!==============================================================================!
!                   Check to see if S is positive definite
!   If removed you will need to allocate llapack arrays before using Hmat
!==============================================================================!
lwork=max(1,3*NG-1)
allocate(work(max(1,lwork)))
call dsyev('v','u',NG,Smat,NG,eigenvalues,work,Lwork,info)
write(*,*) 'Info (Initial Overlap Matrix) ==> ', info
rcn=eigenvalues(1)/eigenvalues(NG)
write(*,*) 'RCN ==> ', RCN
open(unit=20,file='overlap_eigenvalues.dat')
do i=1,NG
  write(20,*) eigenvalues(i)
enddo
close(20)
write(*,*) 'Test 4; Overlap Matrix is Positive Definite'
!!!!!==============================================================================!
!!!!!             Use Gauss Hermit quadrature to evaluate the potential matrix
!!!!!               z(GH-order) w(GH-order) --- quadrature points and weights
!!!!!==============================================================================!
!!!!call cgqf(GH_order,6,0d0,0d0,0d0,1d0,z,w)
!!!!w=w/sqrt(pi)
!!!!!==============================================================================!
!!!!!                   Solve Generalized Eigenvalue Problem
!!!!!==============================================================================!
!!!!do i=1,NG
!!!!  do j=i,NG
!!!!    aij=alpha(i)*alpha(j)/(alpha(i)+alpha(j))
!!!!    r2=sum((x(:,i)-x(:,j))**2)
!!!!    Smat(i,j)=(2*sqrt(alpha(i)*alpha(j))/(alpha(i)+alpha(j)))**(0.5*d)&
!!!!    *exp(-aij*r2)
!!!!    Smat(j,i)=Smat(i,j)
!!!!!==============================================================================!
!!!!!                          Kinetic Energy Matrix
!!!!!==============================================================================!
!!!!    Hmat(i,j)=aij*(d-2*aij*r2)
!!!!    x_ij(:)=(alpha(i)*x(:,i)+alpha(j)*x(:,j))/(alpha(i)+alpha(j))
!!!!    Vij=0d0
!!!!    l(:)=1
!!!!    do ll=1,GH_order**d
!!!!      do k=1,d
!!!!        rr(k)=z(l(k))
!!!!      enddo
!!!!      rr=x_ij+rr/sqrt(alpha(i)+alpha(j))
!!!!      r2=V(rr)
!!!!      do k=1,d
!!!!        r2=r2*w(l(k))
!!!!      enddo
!!!!      Vij=Vij+r2
!!!!      do k=1,d
!!!!        l(k)=mod(l(k),float(GH_order))+1
!!!!        if(l(k).ne.1) exit
!!!!      enddo
!!!!    end do
!!!!!==============================================================================!
!!!!!                      Hamiltonian = Kinetic + Potential
!!!!!==============================================================================!
!!!!    Hmat(i,j)=(Hmat(i,j)+Vij)*Smat(i,j)
!!!!    Hmat(j,i)=Hmat(i,j)
!!!!  enddo
!!!!enddo
!!!!itype=1
!!!!call dsygv(itype,'n','u',NG,Hmat,NG,Smat,NG,eigenvalues,work,Lwork,info)
!!!!write(*,*) 'info ==> ', info
!!!!open(unit=20,file='eigenvalues.dat')
!!!!write(20,*) alpha0, eigenvalues(:)
!!!!close(20)
!!!!!==============================================================================!
!!!!!                              Exact Eigenvalues
!!!!!==============================================================================!
!!!!open(21,File=theory_in)
!!!!do i=1,NG
!!!!  read(21,*) theory(i)
!!!!enddo
!!!!close(21)
!!!!open(unit=22,file='abs_error.dat')
!!!!open(unit=23,file='rel_error.dat')
!!!!open(unit=24,file='alpha_abs_error.dat')
!!!!open(unit=25,file='alpha_rel_error.dat')
!!!!do i=1,NG
!!!!  write(22,*) i, abs(theory(i)-eigenvalues(i))
!!!!  write(23,*) i, (eigenvalues(i)-theory(i))/theory(i)
!!!!  write(24,*) alpha0, abs(theory(i)-eigenvalues(i))
!!!!  write(25,*) alpha0, (eigenvalues(i)-theory(i))/theory(i)
!!!!enddo
!!!!close(22)
!!!!close(23)
!!!!close(24)
!!!!close(25)
!!!!open(unit=26,file='alpha_rel_150.dat')
!!!!do i=1,150
!!!!  write(26,*) alpha0, (eigenvalues(i)-theory(i))/theory(i)
!!!!enddo
!!!!close(26)
!!!!!==============================================================================!
!!!!!                                Output file                                   !
!!!!!==============================================================================!
!!!!open(99,file='simulation.dat')
!!!!write(99,*) 'dimensionality ==> ', d
!!!!write(99,*) 'NG ==> ', NG
!!!!write(99,*) 'alpha0==>', alpha0
!!!!write(99,*) 'RCN Overlap Matrix==>', RCN
!!!!write(99,*) 'GH Order==>', GH_order
!!!!close(99)
!!!!write(*,*) 'Hello Universe!'
end program main
