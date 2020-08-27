!=============================================================================80
!                       QRG-Lmon Water Monomer Implementation                  !
!==============================================================================!
!    Discussion:
!Fortran 90 QRG-Lmon implementation.
!Compute the normalization in a seperate code using MC to normalize
!    Modified:
!21 August 2020
!    Author:
!Shane Flynn
!==============================================================================!
module int_mod
implicit none
!==============================================================================!
!potential            ==>Potential name
!d                    ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!d1                   ==>Monomer Space := 9 for warer
!d2                   ==>Monomer Subspace Lmon-d2
!r(d2,Npoints)        ==>All grid points coordinates (x^i=x_1,x_2,...x_d)
!V_i                  ==>Potential Energy evaluation V(x_i)
!E_cut                ==>Distribution cutoff contour (Kcal/mol input)
!rmin(d)              ==>Minimum of normalization box size
!rmax(d)              ==>Maximum of normalization box size
!N_1D                 ==>Number of points in 1 dimension for computing integral_P
!N_MMC_box            ==>Number of MMC Iterations to determine box-size
!c_LJ                 ==>Parameter for q-LJ pseudo-potential
!integral_P           ==>Normalization constant for the distribtion P(x)
!x0(d)                ==>initial cluster configuration
!c_LJ                 ==>Parameter for q-LJ pseudo-potential
!E_cut                ==>Energy Cutoff Contour (kcal/mol input)
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
double precision::E_cut,integral_P,c_LJ
double precision,allocatable::sqrt_mass(:),mass(:),x0(:),U(:,:),r(:,:)
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
  rr=(x(1:d1)-x0(1:d1))*sqrt_mass(1:d1)  !cartesian==>mass-scaled coordinates
  r_i(1:d2)=0
  do i=1,d2
    r_i(i)=r_i(i)+sum(U(:,i)*rr(:))      !mass-scaled coordinates==>normal modes
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
double precision::Hess_mat(d1,d1),omega(d1)
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
function P_i(r_i,V,x)
!==============================================================================!
!Target Distribution Function, !defined according to a semi-classical argument:
!B. Poirier, “Algebraically self-consistent quasiclassical approximation on
!phase space,” Found. Phys. 30, 1191–1226 (2000).
!==============================================================================!
!x0(d)              ==>initial cluster configuration
!r_i(d2)            ==>i-th grid points coordinates (x^i=x_1,x_2,...x_d)
!V                  ==>Potential Energy evaluation V(x_i)
!E_cut              ==>Distribution cutoff contour
!integral_P         ==>Normalization constant for the distribtion P(x)
!P_i                ==>evaluate P(x)
!==============================================================================!
implicit none
double precision::r_i(d2),V,P_i,x(d),forces(d)
!==============================================================================!
call normal_to_cartesian(r_i,x,.true.)
call water_potential(x,V,forces)
if(V.lt.E_cut) P_i=(E_cut-V)**(d2/2.)/integral_P
if(V.ge.E_cut) P_i=1d-20                      !Define distribution=0 beyond Ecut
end function P_i
!==============================================================================!
function random_integer(Nmin,Nmax)
!==============================================================================!
!Randomly generate an integer in the range Nmin-Nmax
!==============================================================================!
!Nmin           ==>minimum index value
!Nmax           ==>maximum index value
!a              ==>uniform pseudo-random number
!==============================================================================!
implicit none
integer::Nmin,Nmax,random_integer
double precision::a
!==============================================================================!
call random_number(a)
random_integer=floor(a*(Nmax-Nmin+1))+Nmin
end function random_integer
!==============================================================================!
end module int_mod
!==============================================================================!
program main_grid
use int_mod
!==============================================================================!
implicit none
character(len=50)::coord_in
integer::N_MMC_box,N_1D,i,j,my_size
integer(kind=8)::Ntotal,ii
integer,allocatable::seed(:)
double precision::E0,V,time1,time2,dummy,moment,mv_cutoff,Pmax,V_trial
double precision,allocatable,dimension(:)::forces,omega,rmin,rmax,x1,r_i
double precision,allocatable,dimension(:)::delr,index1,r_trial,s
double precision,allocatable,dimension(:)::force0,force1,x
double precision,allocatable,dimension(:,:)::Hess_Mat
!==============================================================================!
!                              Read Input File                                 !
!==============================================================================!
call cpu_time(time1)
read(*,*) d2
read(*,*) E_cut                                        !should be kcal/mol input
read(*,*) c_LJ
read(*,*) N_MMC_box
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
allocate(atom_type(Natoms),mass(Natoms),sqrt_mass(d),x0(d),x1(d),forces(d))
allocate(rmin(d2),rmax(d2),Hess_Mat(d1,d1),omega(d1),U(d1,d1))
allocate(r_i(d2),delr(d2),index1(d2),x(d))
allocate(r_trial(d2),s(d2),force0(d),force1(d))
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
E_cut=E_cut/autokcalmol
!==============================================================================!
call water_potential(x0,E0,forces)
!==============================================================================!
call Get_Hessian(Hess_Mat)
call Mass_Scale_Hessian(Hess_Mat)
call Frequencies_Scaled_Hess(Hess_mat,omega)
!==============================================================================!
!Set seed for random number generator (make it reproducible)
!==============================================================================!
call random_seed(size=my_size)
allocate(seed(my_size))
seed=my_size+1                           !gfortran seed must be larger than size
call random_seed(put=seed)                        !seed must be an integer array
!==============================================================================!
!                     Box Size for normalizing P (MMC)
!Be careful about making mv_cutoff large with MMC, PES does not go to infinity
!==============================================================================!
integral_P=1d0                                     !Initially set to 1 to call P
r_i=0d0                                   !minimum is at 0 in scaled coordinates
rmin=r_i
rmax=r_i
mv_cutoff=0.1
!==============================================================================!
do i=1,N_MMC_box
  call random_number(s)
  r_trial=r_i+mv_cutoff*(2*s-1)                !make trial move (-1,1) not (0,1)
  call normal_to_cartesian(r_trial,x,.true.)
  call water_potential(x,V_trial,forces)
  if(V_trial<E_cut) then
     r_i=r_trial
     do j=1,d2
        if(rmin(j).gt.r_i(j)) rmin(j)=r_i(j)
        if(rmax(j).lt.r_i(j)) rmax(j)=r_i(j)
     enddo
  endif
enddo
!==============================================================================!
Ntotal=(N_1D+1)**d2 !setting Ntotal the same way as square for easy comparison
moment=0d0           !integral_P must be defined to call P_i, use seperate variable
do ii=1,Ntotal
  call random_number(r_trial(:))
  r_trial(:)=rmin(:)+r_trial(:)*(rmax(:)-rmin(:))
  moment=moment+P_i(r_trial,V,x)
enddo
dummy=1./N_1D**d2
do j=1,d2
  dummy=dummy*(rmax(j)-rmin(j))
enddo
integral_P=moment
integral_P=integral_P*dummy
open(90,file='plot.dat')
write(90,*) Ntotal, integral_P
close(90)
!==============================================================================!
!                                 Output
!==============================================================================!
call cpu_time(time2)
open(99,file='out')
write(99,*) 'd ==> ', d
write(99,*) 'd1 ==> ', d1
write(99,*) 'd2 ==> ', d2
write(99,*) 'potential ==> ', potential
do i=1,d2
  write(99,*) 'Box Dimensions==>', rmin(i),rmax(i)
enddo
write(99,*) 'integral_P ==> ',integral_P
write(99,*) 'E_cut ==> ',E_cut
write(99,*) 'N_1D ==> ',N_1D
write(99,*) 'Ntotal==> ',Ntotal
write(99,*) 'c_LJ ==> ',c_LJ
write(99,*) 'N_MMC_box ==> ',N_MMC_box
write(99,*) 'my_size==> ',my_size
write(99,*) 'seed==> ',seed(1)
write(99,*) 'Simulation Time==> ',time2-time1
close(99)
end program main_grid
