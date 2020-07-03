!=================20==================40==================60==================80
!              QRG_Lmon(d2) RoVibrational EigenSpectra using DGB
!==============================================================================!
!Vibrational EigenSpectra calculations for a given potentials QRG
!Generates a DGB for a monomer in a water cluster using the Lmon approximation.
!Basis Functions:
!             phi(r):=(2*alpha_i/pi)^(d/4)*exp[-alpha_i(r-r_i)^2]
!Generates alpha (inverse gaussian width i.e. small alpha=broad gaussian)
!using nearest neighbor (assumes a Quasi-Regular Grid as input)
!Needs gen_hermite_rule.f90 code (Gauss-Hermite-Quadriture) for potential eval.
!Needs LLAPACK to solve the generalized eigenvalue problem
!==============================================================================!
!    Modified:
!30 June 2020
!    Author:
!Shane Flynn
!==============================================================================!
module QRG_Lmon_Spectra
implicit none
!==============================================================================!
!potential            ==>Potential name
!NG                   ==>Number of Gaussians to generate
!d                    ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!d1                   ==>Monomer Space := 9 for warer
!d2                   ==>Monomer Subspace Lmon-d2
!r(d2,Npoints)        ==>All grid points coordinates (x^i=x_1,x_2,...x_d)
!V_i                  ==>Potential Energy evaluation V(x_i)
!x0(d)                ==>initial cluster configuration
!U(d1,d1)             ==>Normal mode eigenvectors
!==============================================================================!
!                            Global Parameters                                 !
!==============================================================================!
double precision,parameter::bohr=0.52917721092
double precision,parameter::autocm=2.194746313D5
double precision,parameter::autokcalmol=627.5096
double precision,parameter::melectron=1822.88839
double precision,parameter::Hmass=1.00782503223*melectron
double precision,parameter::Omass=15.99491461957*melectron
integer,parameter::d1=9
!==============================================================================!
!                            Global Variables                                  !
!==============================================================================!
integer::Natoms,NG,d,d2
character(len=2),allocatable::atom_type(:)
character(len=20)::potential
double precision,allocatable::sqrt_mass(:),mass(:),x0(:),U(:,:)
!==============================================================================!
contains
!==============================================================================!
function Atom_Mass(atom)
!==============================================================================!
!Compute mass of each atom (assumes water as input)
!==============================================================================!
implicit none
double precision::Atom_Mass
character(len=2)::atom
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
subroutine water_potential(x,V,forces)
!==============================================================================!
!To call the water potentials you need to pass in the number of water atoms in
!the system. d=3*Natoms, Nmol=Natoms/3 ==> d/9=Nmol
!==============================================================================!
!x(d)               ==>coordinates
!V                  ==>Potential Energy evaluation V(x)
!forces(d)          ==>Forces from potential
!==============================================================================!
use iso_c_binding
use TIP4P_module
implicit none
double precision::x(d),forces(d),V
!==============================================================================!
if(potential=='tip4p'.or.potential=='TIP4P') then
  call TIP4P(d/9,x,V,forces)
else
  stop 'Cannot Identify Potential, Check "potentials" Subroutine'
endif
end subroutine water_potential
!==============================================================================!
subroutine normal_to_cartesian(r_i,x,flag)
!==============================================================================!
!r_i(1:d2)    ==>x(1:d)     flag=.true.
!x(1:d)       ==>r_i(1:d2)  flag=.false.
!==============================================================================!
!x0(d)        ==>initial cluster configuration
!r_i(d2)      ==>coordinate we want to evaluate the potential at
!x(d)         ==>scaled coordinate (cartesian space) to call potential with
!==============================================================================!
implicit none
integer::i
double precision::x(d),r_i(d2),rr(d1)
logical::flag
if(flag) then
  rr=0
  do i=1,d2
    rr(1:d1)=rr(1:d1)+U(1:d1,i)*r_i(i)
  enddo
  x(1:d)=x0(1:d)
  x(1:d1)=x(1:d1)+rr(1:d1)/sqrt_mass(1:d1)
else
  rr=(x(1:d1)-x0(1:d1))*sqrt_mass(1:d1)  !cartesian-->mass-scaled coordinates
  r_i(1:d2)=0
  do i=1,d2
    r_i(i)=r_i(i)+sum(U(:,i)*rr(:))      !mass-scaled coordinates-->normal modes
  enddo
endif
end subroutine normal_to_cartesian
!==============================================================================!
subroutine Get_Hessian(Hess_Mat)
!==============================================================================!
!Numerically evaluate the Hessian
!==============================================================================!
!potential          ==>potential name
!d                  ==>Total System Dimensionality  (d:=3*Natoms)
!x0(d)              ==>Initial Configuration (entire system)
!force(d)           ==>Forces from potential
!E0                 ==>Potential Energy of x0
!s                  ==>Perturbation Parameter
!Hess_Mat(d1,d1)      ==>Numerical Hessian
!==============================================================================!
implicit none
integer::i,j
double precision::Hess_Mat(d1,d1),x1(d),force0(d),force1(d),E0
double precision,parameter::ss=1d-6
x1=x0
call water_potential(x1,E0,force0)
do i=1,d1
  x1(i)=x0(i)+ss
  call water_potential(x1,E0,force1)
  x1(i)=x0(i)
  do j=1,d1
    Hess_Mat(i,j)=(force0(j)-force1(j))/ss
  enddo
enddo
end subroutine Get_Hessian
!==============================================================================!
subroutine Mass_Scale_Hessian(Hess_Mat)
!==============================================================================!
!Symmetrize and Mass-Scale the Hessian
!==============================================================================!
!d                  ==>Total System Dimensionality  (d:=3*Natoms)
!Hess_Mat(d,d)      ==>Numerical Hessian
!sqrt_mass(d)       ==>Square Root Mass
!==============================================================================!
implicit none
integer::i,j
double precision::Hess_Mat(d1,d1)
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
subroutine reverse(N,A)
!==============================================================================!
! Reverse the elements in array A(N)
!==============================================================================!
integer::N,i
double precision::A(N),temp
do i=1,N/2
  temp=A(i)
  A(i)=A(N-i+1)
  A(N-i+1)=temp
enddo
end subroutine reverse
!==============================================================================!
subroutine Frequencies_Scaled_Hess(Hess_mat,omega)
!==============================================================================!
!Compute Eigenvalues and Eigenvectors for the mass-scaled hessian
!Uses the LLAPACK real symmetric eigen-solver (dsygev)
!==============================================================================!
!d                  ==>Total System Dimensionality  (d:=3*Natoms)
!Hess_Mat(d1,d1)    ==>Numerical Hessian
!omega(d)           ==>Eigenvalues
!U(d,d)             ==>Eigenvectors
!     LLAPACK(dsyev):
!v                  ==>Compute both Eigenvalues and Eigenvectors
!u                  ==>Use Upper-Triangle of matrix
!Lwork              ==>Allocation size
!==============================================================================!
implicit none
integer::i,info,Lwork
double precision::Hess_mat(d1,d1),omega(d1)
double precision,allocatable::work(:)
Lwork=max(1,3*d1-1)
allocate(work(max(1,Lwork)))
U=Hess_mat
call dsyev('v','u',d1,U,d1,omega,work,Lwork,info)
!==============================================================================!
!                   sqrt hessian matrix eigenvalues
!==============================================================================!
do i=1,d1
  if(omega(i)<0d0) write(*,*) 'Warning:  lambda(',i,')=',omega(i)
  omega(i)=sign(sqrt(abs(omega(i))),omega(i))
enddo
!==============================================================================!
!Subspace Needs Largest Eigenvalues: llapack outputs small to large ==>re-order
!==============================================================================!
call reverse(d1,omega)
do i=1,d1
  call reverse(d1,U(i,:))
enddo
open(18,File='freq_scaled_hess.dat')
do i=1,d1
  write(18,*) omega(i), 'normalized = 1?', sum(U(:,i)**2)
enddo
close(18)
end subroutine Frequencies_Scaled_Hess
!==============================================================================!
end module QRG_Lmon_Spectra
!==============================================================================!
program main_spectra
use QRG_Lmon_Spectra
!==============================================================================!
implicit none
character(len=50)::grid_in,coord_in
integer::monomer,Lwork,GH_order,info,i,j,k,ll,itype
double precision::E0,E,alpha0,V_ij,r2,a_ij
double precision,allocatable,dimension(:)::forces,omega,alpha,eigenvalues,z,w,rr
double precision,allocatable,dimension(:)::r_ij,l,work,x
double precision,allocatable,dimension(:,:)::Hess_Mat,Smat,Hmat,r
double precision,parameter::pi=4.*atan(1d0)
!==============================================================================!
!                             Read Input Data File
!==============================================================================!
read(*,*) d2
read(*,*) potential
read(*,*) coord_in
read(*,*) grid_in
read(*,*) NG
read(*,*) GH_order
read(*,*) alpha0
read(*,*) monomer               !which monomer in the cluster you are looking at
!==============================================================================!
!                                  Read xyz
!==============================================================================!
open(17,file=coord_in)
read(17,*) Natoms
read(17,*)
d=3*Natoms
!==============================================================================!
allocate(atom_type(Natoms),mass(Natoms),sqrt_mass(d),x0(d),forces(d),alpha(NG))
allocate(Hess_Mat(d1,d1),omega(d1),U(d1,d1),r(d1,NG),rr(d1),eigenvalues(NG))
allocate(Smat(NG,NG),Hmat(NG,NG),z(GH_order),w(GH_order),l(d1),r_ij(d1),x(d))
do i=1,Natoms
  read(17,*) atom_type(i),x0(3*i-2:3*i)                !input is xyz therefore 3
  mass(i)=Atom_mass(atom_type(i))
  sqrt_mass(3*i-2:3*i)=sqrt(mass(i))
enddo
close(17)
!==============================================================================!
!        Input coordinates are in angstroms, convert to atomic units
!==============================================================================!
x0=x0/bohr
open(18,File='cluster_initial.dat')
write(18,*) 'x0 in atomic units (xo/bohr)', x0
call water_potential(x0,E0,forces)
write(18,*) 'E0 (atomic) ==> ', E0
write(18,*) 'E0 (cm-1) ==> ', E0*autocm
write(18,*) 'E0 (kcal/mol) ==> ', E0*autokcalmol
close(18)
!==============================================================================!
call Get_Hessian(Hess_Mat)
call Mass_Scale_Hessian(Hess_Mat)
call Frequencies_Scaled_Hess(Hess_mat,omega)
!==============================================================================!
r=0                                                           !set (d2+1:d1) = 0
open(19,file=grid_in)
do i=1,NG
  read(19,*) r(1:d2,i)                         !grid is generated in d2 subspace
enddo
close(19)
!==============================================================================!
!         Symmetric gaussians. Use nearest neighbor to determine alpha
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
!==============================================================================!
!Compute Hamiltonian Matrix. Gauss-Hermite quadrature for potential
!==============================================================================!
call cgqf(GH_order,6,0d0,0d0,0d0,1d0,z,w)                                   !GHQ
w=w/sqrt(pi)
!==============================================================================!
do i=1,NG
  do j=i,NG
    a_ij=alpha(i)*alpha(j)/(alpha(i)+alpha(j))
    r2=sum((r(:,i)-r(:,j))**2)
!==============================================================================!
!               Compute i-jth element of the overlap matrix
!==============================================================================!
    Smat(i,j)=(2.*sqrt(alpha(i)*alpha(j))/(alpha(i)+alpha(j)))**(0.5*d2)&
    *exp(-a_ij*r2)
    Smat(j,i)=Smat(i,j)                   !Pass Full Smat to compute Eigenvalues
!==============================================================================!
!              Compute i-jth element of the kinetic matrix
!==============================================================================!
    Hmat(i,j)=a_ij*(d2-2.*a_ij*r2)
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
      call normal_to_cartesian(rr,x,.true.)
      call water_potential(x,E,forces)
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
!==============================================================================!
!               Eigenvalues of the Hamiltonian matrix
!==============================================================================!
Lwork=max(1,3*NG-1)                                          !LLAPACK Suggestion
allocate(work(max(1,Lwork)))
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
write(24,*) alpha0, (eigenvalues(5)-eigenvalues(1))*autocm
close(22)
close(23)
close(24)
!==============================================================================!
!                               Out File
!==============================================================================!
open(99,file='out')
write(99,*) 'd ==> ', d
write(99,*) 'd1 ==> ', d1
write(99,*) 'd2 ==> ', d2
write(99,*) 'potential ==> ', potential
write(99,*) 'Number of atoms==> ', Natoms
write(99,*) 'Number of Gaussians==> ', NG
write(99,*) 'alpha0==>', alpha0
write(99,*) 'GH Order==>', GH_order
write(99,*) 'E0==>', E0
write(99,*) 'Initial configuration==> ', x0
close(99)
!==============================================================================!
end program main_spectra
