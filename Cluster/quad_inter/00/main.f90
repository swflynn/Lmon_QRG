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
!This code is written using TIP4P potential only.
!==============================================================================!
!    Modified:
!26 October 2020
!    Author:
!Shane Flynn
!==============================================================================!
module QRG_Lmon_Spectra
implicit none
!==============================================================================!
!potential            ==>Potential name
!atom_type(Natoms)    ==>Atom name short-hand notation
!Natoms               ==>Number of atoms (line 1 in xyz file format)
!NG                   ==>Number of Gaussians to generate
!d                    ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!d1                   ==>Monomer Space := 9 for warer
!d2                   ==>Monomer Subspace Lmon-d2
!x0(d)                ==>Initial cluster configuration
!mass(Natoms)         ==>Atom masses
!sqrt_mass(d)         ==>Square root atom masses
!alpha(NG)            ==>Gaussian width parameter
!r(d2,Npoints)        ==>All grid points coordinates (x^i=x_1,x_2,...x_d)
!U(d1,d1)             ==>Normal mode eigenvectors
!==============================================================================!
!                            Global Parameters                                 !
!==============================================================================!
double precision,parameter::pi=4.*atan(1d0)
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
character(len=20)::potential
character(len=2),allocatable::atom_type(:)
integer::Natoms,NG,d,d2
double precision,allocatable,dimension(:)::x0,mass,sqrt_mass,alpha
double precision,allocatable,dimension(:,:)::r,U
!==============================================================================!
contains
!==============================================================================!
function Atom_Mass(atom)
!==============================================================================!
!Compute mass of each atom, assumes water as input (H,O only)
!==============================================================================!
implicit none
double precision::Atom_Mass
character(len=2)::atom
!==============================================================================!
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
!To call a water potential; need the number of water molecules in the system
!d=3*Natoms, Nmolecules=Natoms/3 ==> d/9=Nmolecules
!==============================================================================!
!x(d)                 ==>Cartesian coordinates
!V                    ==>Potential Energy evaluation V(x)
!forces(d)            ==>Forces (computed by the PES)
!==============================================================================!
use iso_c_binding
use TIP4P_module
implicit none
double precision::x(d),V,forces(d)
!==============================================================================!
if(potential=='tip4p'.or.potential=='TIP4P') then
  call TIP4P(d/9,x,V,forces)
!else if(potential == 'mbpol' .or. potential == 'MBPOL') then
  !call calcpotg(d/9, V, x*bohr, forces)
  !forces=-forces*bohr/autokcalmol
  !V=V/autokcalmol
else
  stop 'Cannot Identify Potential, Check "potentials" Subroutine'
endif
end subroutine water_potential
!==============================================================================!
subroutine normal_to_cartesian(r,x,flag)
!==============================================================================!
!Convert normal-mode coordinates to mass-scaled cartesian coordinates
!The PES must be called using cartesian coordinates
!if flag==true;     r(1:d2)==>x(1:d)
!if flag==false;    x(1:d)==>r(1:d2)
!==============================================================================!
!r(d2)                ==>normal mode coordinate
!x(d)                 ==>mass-scaled cartesian coordinate
!==============================================================================!
implicit none
logical::flag
integer::i
double precision::x(d),r(d2),rr(d1)
!==============================================================================!
if(flag)then
  rr=0
  do i=1,d2
    rr(1:d1)=rr(1:d1)+U(1:d1,i)*r(i)
  enddo
  x(1:d)=x0(1:d)
  x(1:d1)=x(1:d1)+rr(1:d1)/sqrt_mass(1:d1)
else
  rr=(x(1:d1)-x0(1:d1))*sqrt_mass(1:d1)     !cartesian==>mass-scaled coordinates
  r(1:d2)=0
  do i=1,d2
    r(i)=r(i)+sum(U(1:d1,i)*rr(1:d1))    !mass-scaled coordinates==>normal modes
  enddo
endif
end subroutine normal_to_cartesian
!==============================================================================!
subroutine gradient_transform(force,force1)
!==============================================================================!
!Convert force from mass-scaled cartesian to normal-mode coordinates (force1)
!==============================================================================!
!force(d)             ==>forces in mass-scaled cartesian coordinates
!force1(d2)           ==>forces in normal-mode coordinates
!==============================================================================!
implicit none
integer::i
double precision::force(d),force1(d2),rr(d1)
!==============================================================================!
rr(1:d1)=force(1:d1)/sqrt_mass(1:d1)
force1(1:d2)=0
do i=1,d2
  force1(i)=force1(i)+sum(U(1:d1,i)*rr(1:d1))        !mass-scaled==>normal modes
enddo
end subroutine gradient_transform
!==============================================================================!
subroutine matrix_elements_grad(i,j,S_ij,H_ij,GH_order,z,w)
!==============================================================================!
!Compute the i-jth matrix elements for the Overlap,Kinetic,Hamiltonian.
!GH_Order<0; use quadratic interpolation to compute potential. Else quadrature
!Uses the gradient, use "Matrix Elements" if gradient is not available/costly
!==============================================================================!
!S_ij                 ==>i-jth element of the Overlap Matrix
!H_ij                 ==>i-jth element of the Hamiltonian Matrix
!GH_order             ==>Gauss-Hermite Quadrature Order
!ss                   ==>Perturbation Parameter
!Lap                  ==>Laplace (see equations)
!z(GH_order)          ==>Quadrature points
!w(GH_order)          ==>Quadrature weights
!==============================================================================!
integer::i,j,k,ll,l(d2),GH_order
double precision,parameter::ss=1d-6            !Should test for desired accuracy
double precision::S_ij,H_ij,r2,a_ij,E,V_ij,Lap,ss1
double precision::r_ij(d2),rr(d2),w(GH_order),z(GH_order),r1(d2),x(d),force(d)
double precision::force1(d2),force0(d2)
!==============================================================================!
ss1=ss*4*(alpha(i)+alpha(j))
a_ij=alpha(i)*alpha(j)/(alpha(i)+alpha(j))
r2=sum((r(:,i)-r(:,j))**2)
!==============================================================================!
!               Compute i-jth element of the overlap matrix
!==============================================================================!
S_ij=(2.*sqrt(alpha(i)*alpha(j))/(alpha(i)+alpha(j)))**(0.5*d2)*exp(-a_ij*r2)
!==============================================================================!
!               Compute i-jth element of the kinetic matrix
!==============================================================================!
H_ij=a_ij*(d2-2.*a_ij*r2)
!==============================================================================!
!    Compute Potential Energy by expanding the potential up to second order
!==============================================================================!
r_ij(:)=(alpha(i)*r(:,i)+alpha(j)*r(:,j))/(alpha(i)+alpha(j))
!if(i.ne.j) then                               !use quadratic interpolation of V
if(GH_order<0) then                            !use quadratic interpolation of V
  r1=r_ij
  call normal_to_cartesian(r1,x,.true.)
  call water_potential(x,V_ij,force)
  call gradient_transform(force,force0)
  Lap=0d0
  do k=1,d2
    r1(k)=r_ij(k)+ss
    call normal_to_cartesian(r1,x,.true.)
    call water_potential(x,E,force)
    call gradient_transform(force,force1)
    Lap=Lap+force0(k)-force1(k)
    r1=r_ij
  enddo
  V_ij=V_ij+Lap/ss1
else                                 !use quadrature to compute potential energy
  V_ij=0d0
  l(:)=1
  do ll=1,GH_order**d2
    do k=1,d2
      rr(k)=z(l(k))
    enddo
    rr=r_ij+rr/sqrt(alpha(i)+alpha(j))
    call normal_to_cartesian(rr,x,.true.)
    call water_potential(x,E,force)
    do k=1,d2
      E=E*w(l(k))
    enddo
    V_ij=V_ij+E
    do k=1,d2
      l(k)=mod(l(k),GH_order)+1
      if(l(k).ne.1) exit
    enddo
  enddo
endif
H_ij=(H_ij+V_ij)*S_ij
end subroutine matrix_elements_grad
!==============================================================================!
subroutine matrix_elements(i,j,S_ij,H_ij,GH_order,z,w)
!==============================================================================!
!Compute the i-jth matrix elements for the Overlap,Kinetic,Hamiltonian.
!GH_Order<0; use quadratic interpolation to compute potential. Else quadrature
!Does not require the gradient to compute potential
!==============================================================================!
!S_ij                 ==>i-jth element of the Overlap Matrix
!H_ij                 ==>i-jth element of the Hamiltonian Matrix
!GH_order             ==>Gauss-Hermite Quadrature Order
!ss                   ==>Perturbation Parameter
!Lap                  ==>Laplace (see equations)
!z(gh_order)          ==>quadrature points
!w(GH_order)          ==>Quadrature weights
!==============================================================================!
integer::i,j,k,ll,l(d2),GH_order
double precision,parameter::ss=1d-6            !Should test for desired accuracy
double precision::S_ij,H_ij,r2,a_ij,E,E1,E2,V_ij,Lap,ss1
double precision::r_ij(d2),rr(d2),w(GH_order),z(GH_order),r1(d2),x(d),force(d)
!==============================================================================!
ss1=ss**2*4*(alpha(i)+alpha(j))
a_ij=alpha(i)*alpha(j)/(alpha(i)+alpha(j))
r2=sum((r(:,i)-r(:,j))**2)
!==============================================================================!
!             Compute i-jth element of the overlap matrix
!==============================================================================!
S_ij=(2.*sqrt(alpha(i)*alpha(j))/(alpha(i)+alpha(j)))**(0.5*d2)*exp(-a_ij*r2)
!==============================================================================!
!             Compute i-jth element of the kinetic matrix
!==============================================================================!
H_ij=a_ij*(d2-2.*a_ij*r2)
!==============================================================================!
!    Compute Potential Energy by expanding the potential up to second order
!============4b==================================================================!
r_ij(:)=(alpha(i)*r(:,i)+alpha(j)*r(:,j))/(alpha(i)+alpha(j))
!if(i.ne.j) then                               !use quadratic interpolation of V
if(GH_order<0)then                             !use quadratic interpolation of V
  r1=r_ij
  call normal_to_cartesian(r1,x,.true.)
  call water_potential(x,E,force)
  Lap=0d0
  do k=1,d2
    r1(k)=r_ij(k)+ss
    call normal_to_cartesian(r1,x,.true.)
    call water_potential(x,E1,force)
    r1(k)=r_ij(k)-ss
    call normal_to_cartesian(r1,x,.true.)
    call water_potential(x,E2,force)
    Lap=Lap-2*E+E1+E2
    r1(k)=r_ij(k)
  enddo
  V_ij=E+Lap/(ss**2*4*(alpha(i)+alpha(j)))
else                                                             !use quadrature
  V_ij=0d0
  l(:)=1
  do ll=1,GH_order**d2
    do k=1,d2
      rr(k)=z(l(k))
    enddo
    rr=r_ij+rr/sqrt(alpha(i)+alpha(j))
    call normal_to_cartesian(rr,x,.true.)
    call water_potential(x,E,force)
    do k=1,d2
      E=E*w(l(k))
    enddo
    V_ij=V_ij+E
    do k=1,d2
      l(k)=mod(l(k),GH_order)+1
      if(l(k).ne.1) exit
    enddo
  enddo
endif
H_ij=(H_ij+V_ij)*S_ij
end subroutine matrix_elements
!==============================================================================!
subroutine Get_Hessian(Hess_Mat)
!==============================================================================!
!Numerically evaluate the Hessian
!==============================================================================!
!potential            ==>potential name
!d                    ==>Total System Dimensionality  (d:=3*Natoms)
!x0(d)                ==>Initial Configuration (entire system)
!force(d)             ==>Forces from potential
!E0                   ==>Potential Energy of x0
!s                    ==>Perturbation Parameter
!Hess_Mat(d1,d1)      ==>Numerical Hessian
!==============================================================================!
implicit none
integer::i,j
double precision::Hess_Mat(d1,d1),x1(d),force0(d),force1(d),E0
double precision,parameter::ss=1d-6
!==============================================================================!
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
!d                    ==>Total System Dimensionality  (d:=3*Natoms)
!Hess_Mat(d1,d1)      ==>Numerical Hessian
!sqrt_mass(d)         ==>Square Root Mass
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
!Reverse the elements in array A(N)
!==============================================================================!
integer::N,i
double precision::A(N),temp
!==============================================================================!
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
!sqrt_mass(d)         ==>Square Root Mass
!Hess_Mat(d1,d1)      ==>Numerical Hessian
!omega(d)             ==>Eigenvalues
!U(d,d)               ==>Eigenvectors
!     LLAPACK(dsyev):
!v                    ==>Compute both Eigenvalues and Eigenvectors
!u                    ==>Use Upper-Triangle of matrix
!Lwork                ==>Allocation size; LLAPACK Suggetion: Lwork=max(1,3*d1-1)
!==============================================================================!
implicit none
integer::i,info,Lwork
double precision::Hess_mat(d1,d1),omega(d1),dummy(d1)
double precision,allocatable::work(:)
!==============================================================================!
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
  dummy(:)=U(i,:)
  call reverse(d1,dummy)
  U(i,:)=dummy(:)
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
!grid_in              ==>Quasi-Regular Grid file name
!coord_in             ==>Water (minimized structure) XYZ filename
!GH_order             ==>Gauss-Hermite Quadrature Order
!alpha0               ==>Flat scaling parameter for gaussian widths
!omega(d1)            ==>Mass-Scaled Hessian eigenvalues
!eigenvalues(NG)      ==>Hamiltonian eigenvalues
!z(GH_order)          ==>quadrature points
!w(GH_order)          ==>Quadrature weights
!Smat(NG,NG)          ==>Overlap Matrix
!Hmat(NG,NG)          ==>Hamiltonian Matrix
!Hess_Mat(d1,d1)      ==>Hamiltonian Matrix
!force(d)             ==>Forces (computed by water potential)
!==============================================================================!
implicit none
character(len=50)::coord_in,grid_in
integer::Lwork,info,i,j,itype,GH_order
double precision::E0,alpha0,time1,time2,r2
double precision,allocatable,dimension(:)::omega,eigenvalues,z,w,force,work
double precision,allocatable,dimension(:,:)::Hess_Mat,Smat,Hmat
call cpu_time(time1)
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
!==============================================================================!
!                                  Read xyz
!==============================================================================!
open(17,file=coord_in)
read(17,*) Natoms
read(17,*)
d=3*Natoms
!==============================================================================!
allocate(atom_type(Natoms),mass(Natoms),sqrt_mass(d),x0(d),force(d))
allocate(Hess_Mat(d1,d1),U(d1,d1),omega(d1),alpha(NG),eigenvalues(NG))
allocate(Smat(NG,NG),Hmat(NG,NG),z(GH_order),w(GH_order),r(d2,NG))
!==============================================================================!
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
call water_potential(x0,E0,force)
write(18,*) 'E0 (atomic) ==> ', E0
write(18,*) 'E0 (cm-1) ==> ', E0*autocm
write(18,*) 'E0 (kcal/mol) ==> ', E0*autokcalmol
close(18)
!==============================================================================!
call Get_Hessian(Hess_Mat)
call Mass_Scale_Hessian(Hess_Mat)
call Frequencies_Scaled_Hess(Hess_mat,omega)
!==============================================================================!
open(19,file=grid_in)
do i=1,NG
  read(19,*) r(:,i)                         !grid is generated in d2 subspace
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
if(GH_order>0) then
call cgqf(GH_order,6,0d0,0d0,0d0,1d0,z,w)                                   !GHQ
w=w/sqrt(pi)
endif
!==============================================================================!
do i=1,NG
  do j=i,NG
!==============================================================================!
!Choose when to use quadrature vs quadratic interpolation
!suggest quadrature for diagonal elements for improved accuracy
!GH_order<0;use quadratic interpolation
!==============================================================================!
    if(i==j) then
      call matrix_elements(i,j,Smat(i,j),Hmat(i,j),GH_order,z,w)
      !call matrix_elements(i,j,Smat(i,j),Hmat(i,j),-1,z,w)
    else
      !call matrix_elements(i,j,Smat(i,j),Hmat(i,j),GH_order,z,w)
      call matrix_elements(i,j,Smat(i,j),Hmat(i,j),-1,z,w)
    endif
    Hmat(j,i)=Hmat(i,j)
    Smat(j,i)=Smat(i,j)
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
call cpu_time(time2)
open(99,file='out')
write(99,*) 'Simulation Time==> ',time2-time1
write(99,*) 'd ==> ', d
write(99,*) 'd1 ==> ', d1
write(99,*) 'd2 ==> ', d2
write(99,*) 'potential ==> ', potential
write(99,*) 'Number of atoms==> ', Natoms
write(99,*) 'Number of Gaussians==> ', NG
write(99,*) 'alpha0==>', alpha0
write(99,*) 'GH Order ==>', GH_order
write(99,*) 'E0==>', E0
write(99,*) 'Initial configuration==> ', x0
close(99)
!==============================================================================!
end program main_spectra
