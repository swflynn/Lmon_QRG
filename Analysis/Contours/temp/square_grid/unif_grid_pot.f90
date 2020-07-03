!=============================================================================80
!                    QRG-Lmon Water Potential Evaluation
!==============================================================================!
!    Discussion:
!Generate a uniform grid and return the normal mode/cartesian
!coordiantes+potential (r1,r2,V) (x1,x2,V)
!==============================================================================!
!    Modified:
!26 June 2020
!    Author:
!Shane Flynn
!==============================================================================!
module QRG_Lmon_Grid
implicit none
!==============================================================================!
!potential            ==>Potential name
!d                    ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!d1                   ==>Monomer Space := 9 for warer
!d2                   ==>Monomer Subspace Lmon-d2
!r(d2,N_1D**2)        ==>All grid points  (direct product)
!V                    ==>Potential Energy evaluation V(x_i)
!N_1D                 ==>Number of points for a single dimension
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
integer::Natoms,d,d2
character(len=2),allocatable::atom_type(:)
character(len=20)::potential
double precision,allocatable::sqrt_mass(:),mass(:),x0(:)
logical::writer                       !if true write cartesian coordinates and V
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
double precision::x(d),V,forces(d)
!==============================================================================!
if(potential=='tip4p'.or.potential=='TIP4P') then
  call TIP4P(d/9,x,V,forces)
else
  stop 'Cannot Identify Potential, Check "potentials" Subroutine'
endif
if (writer .eqv. .true.) then
  write(*,*) x(1:d2), V
endif
end subroutine water_potential
!==============================================================================!
subroutine normal_cartesian_potential(r_i,V,U)
!==============================================================================!
!Evaluate potential in d-dimensionality. Lmon moves 1 monomer of the cluster
!We only have dynamics over d2 subspace of the d1 monomer
!==============================================================================!
!x0(d)        ==>initial cluster configuration
!V            ==>Potential Energy Evalation
!U            ==>Normal mode eigenvectors
!r_i(d2)      ==>coordinate we want to evaluate the potential at
!r(d1)        ==>need d1 size to evaluate potential
!x(d)         ==>scaled coordinate (cartesian space) to call potential with
!==============================================================================!
implicit none
double precision::x(d),r_i(d2),rr(d1),forces(d),V,U(d1,d1)
rr=0               !need d1-dimensional coordinate for monomer, (1:d2)=r, rest=0
rr(1:d2)=r_i(1:d2)
x=x0
!x(1:d1)=x(1:d1)+matmul(U,rr)/sqrt_mass(1:d1)
!Test different order for linear algebra, rr,U is much faster....
x(1:d1)=x(1:d1)+matmul(rr,U)/sqrt_mass(1:d1)
call water_potential(x,V,forces)
end subroutine normal_cartesian_potential
!==============================================================================!==============================================================================!
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
subroutine Frequencies_Scaled_Hess(Hess_mat,omega,U)
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
integer::i,info,Lwork
double precision::Hess_mat(d1,d1),omega(d1),U(d1,d1)
double precision,allocatable::work(:)
lwork=max(1,3*d1-1)
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
end module QRG_Lmon_Grid
!==============================================================================!
program main_grid
use QRG_Lmon_Grid
!==============================================================================!
implicit none
character(len=50)::coord_in
integer::N_1D,counter,i,j
double precision::lower,upper,V
double precision,allocatable,dimension(:)::forces,omega,points
double precision,allocatable,dimension(:,:)::Hess_Mat,U,r
!==============================================================================!
!                              Read Input File                                 !
!==============================================================================!
read(*,*) d2
read(*,*) upper
read(*,*) lower
read(*,*) N_1D
read(*,*) potential
read(*,*) coord_in
!==============================================================================!
!                                  Read xyz
!==============================================================================!
open(17,file=coord_in)
read(17,*) Natoms
read(17,*)
d=3*Natoms
!==============================================================================!
allocate(atom_type(Natoms),mass(Natoms),sqrt_mass(d),x0(d),forces(d))
allocate(Hess_Mat(d1,d1),omega(d1),U(d1,d1))
allocate(r(d2,N_1D**2),points(N_1D))
do i=1,Natoms
  read(17,*) atom_type(i),x0(3*i-2:3*i)                            !input is xyz
  mass(i)=Atom_mass(atom_type(i))
  sqrt_mass(3*i-2:3*i)=sqrt(mass(i))
enddo
close(17)
!==============================================================================!
!        Input coordinates are in angstroms, convert to atomic units
!==============================================================================!
x0=x0/bohr
!==============================================================================!
writer=.false.                                    !don't write cartesian coords
call Get_Hessian(Hess_Mat)
call Mass_Scale_Hessian(Hess_Mat)
call Frequencies_Scaled_Hess(Hess_mat,omega,U)
do i=1,N_1D
    points(i)=lower+(i-1.)*(upper-lower)/(N_1D-1.)
enddo
counter=1
do i=1,N_1D
    do j=1,N_1D
        r(1,counter)=points(i)
        r(2,counter)=points(j)
        counter=counter+1
    enddo
enddo

writer=.true.                               !write cartesian coords and potential
open(55,File='normal_pot.dat')
do i=1,N_1D**2
  call normal_cartesian_potential(r(:,i),V,U)
  write(55,*) r(:,i), V
enddo
close(55)
end program main_grid
