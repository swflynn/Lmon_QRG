!=============================================================================80
!                       QRG-Lmon Water Monomer Implementation                  !
!==============================================================================!
!    Discussion:
!Fortran 90 QRG-Lmon implementation.
!Generates QRG for a monomer subspace using the Lmon approximation.
!    Modified:
!30 June 2020
!    Author:
!Shane Flynn
!==============================================================================!
module QRG_Lmon_Grid
implicit none
!==============================================================================!
!potential            ==>Potential name
!Npoints              ==>Number of points to generate
!d                    ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!d1                   ==>Monomer Space := 9 for warer
!d2                   ==>Monomer Subspace Lmon-d2
!r(d2,Npoints)        ==>All grid points coordinates (x^i=x_1,x_2,...x_d)
!Uij(Npoints,Npoints) ==>Pairwise energies for all i-j points
!V_i                  ==>Potential Energy evaluation V(x_i)
!E_cut                ==>Distribution cutoff contour (Kcal/mol input)
!rmin(d)              ==>Minimum of normalization box size
!rmax(d)              ==>Maximum of normalization box size
!N_1D                 ==>Number of points in 1 dimension for computing integral_P
!N_MMC_box            ==>Number of MMC Iterations to determine box-size
!c_LJ                 ==>Parameter for q-LJ pseudo-potential
!N_MMC_grid           ==>Number of MMC Iterations to optimize QRG
!MMC_freq             ==>Frequency to update QRG MMC grid mv_cutoff
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
integer::Natoms,Npoints,d,d2
character(len=2),allocatable::atom_type(:)
character(len=20)::potential
double precision::E_cut,integral_P,c_LJ
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
function Pair_LJ_NRG(r_i,r_j,V)
!==============================================================================!
!Computes the quasi-Lennard Jones pairwise energy between grid points used in
!our QRG algorithm
!This function computes the q-LJ energy between 2 grid-points
!==============================================================================!
!potential          ==>Potential name
!Npoints            ==>Number of points to generate
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!r_i(d2)            ==>Grid points coordinates (x^i=x_1,x_2,...x_d)
!V_i                ==>Potential Energy evaluation V(x_i)
!E_cut              ==>Distribution cutoff contour
!c_LJ               ==>Parameter for q-LJ pseudo-potential
!integral_P         ==>Normalization constant for the distribtion P(x)
!Pair_LJ_NRG        ==>Energy of the i-j q-LJ potential
!==============================================================================!
implicit none
double precision::r_i(d2),r_j(d2),V,a,b,Pair_LJ_NRG,sigma1,sigma2,x(d)
!==============================================================================!
a=sum((r_i(:)-r_j(:))**2)
sigma1=c_LJ*(P_i(r_i,V,x)*Npoints)**(-1./d2)
sigma2=c_LJ*(P_i(r_j,V,x)*Npoints)**(-1./d2)
b=(sigma2**2/a)
a=(sigma1**2/a)
Pair_LJ_NRG=a**(d2+9)-a**(d2+3)+b**(d2+9)-b**(d2+3)
end function Pair_LJ_NRG
!==============================================================================!
end module QRG_Lmon_Grid
!==============================================================================!
program main_grid
use QRG_Lmon_Grid
!==============================================================================!
implicit none
character(len=50)::coord_in
integer::N_MMC_box,N_1D,N_MMC_grid,MMC_freq,accept,counter,i,j,k
integer::Ntotal
double precision::E0,V,time1,time2,Delta_E,deltae1,dummy,moment,mv_cutoff
double precision,allocatable,dimension(:)::forces,omega,rmin,rmax,x1,r_i
double precision,allocatable,dimension(:)::delr,index1,U_move,rr,r_trial,s
double precision,allocatable,dimension(:)::force0,force1,x
double precision,allocatable,dimension(:,:)::Hess_Mat,r,Uij
!==============================================================================!
!                              Read Input File                                 !
!==============================================================================!
call cpu_time(time1)
read(*,*) Npoints
read(*,*) d2
read(*,*) N_MMC_box
read(*,*) E_cut                                        !should be kcal/mol input
read(*,*) N_1D
read(*,*) c_LJ
read(*,*) N_MMC_grid
read(*,*) MMC_freq
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
allocate(rmin(d2),rmax(d2),Hess_Mat(d1,d1),omega(d1),U(d1,d1),r(d2,Npoints))
allocate(Uij(Npoints,Npoints),r_i(d2),delr(d2),index1(d2),x(d))
allocate(U_move(Npoints),rr(d2),r_trial(d2),s(d2),force0(d),force1(d))
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
open(18,File='cluster_initial.dat')
write(18,*) 'x0 in atomic units (xo/bohr)', x0
write(18,*) 'E_cut in atomic units (assumed kcal/mol input)', E_cut
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
!                     Box Size for normalizing P (MMC)
!==============================================================================!
integral_P=1d0                                     !Initially set to 1 to call P
r_i=0d0
rmin=r_i
rmax=r_i
mv_cutoff=0.1
do i=1,N_MMC_box
  call random_number(s)
  r_trial=r_i+mv_cutoff*(2*s-1)     !trial move (-1,1), random numbers are (0,1)
  call random_number(dummy)                             !MMC acceptance criteria
  if(P_i(r_trial,V,x)/P_i(r_i,V,x).ge.dummy) then
    r_i=r_trial
    do j=1,d2
      if(rmin(j).gt.r_i(j)) rmin(j)=r_i(j)
      if(rmax(j).lt.r_i(j)) rmax(j)=r_i(j)
    enddo
  endif
enddo
!==============================================================================!
!Compute Integral P with square grid         P(x)~Area_Square/N sum_n=1,N P(x_n)
!==============================================================================!
open(20,File='direct_grid.dat')
Moment=0.
Ntotal=(N_1D+1)**d2
index1=0
delr(:)=(rmax(:)-rmin(:))/N_1D
do i=1,Ntotal
  do j=1,d2
    if(index1(j).eq.N_1D) then
      index1(j)=0
    else
      index1(j)=index1(j)+1
      exit
    endif
  enddo
  r_i(:)=rmin(:)+index1(:)*delr(:)
  dummy=P_i(r_i,V,x)
  Moment=Moment+dummy
  if(V.lt.E_cut) write(20,*) r_i
enddo
dummy=1./N_1D**d2
do j=1,d2
  dummy=dummy*(rmax(j)-rmin(j))
enddo
integral_P=dummy*Moment
close(20)
!==============================================================================!
!           Generate initial distribution to then convert to a QRG
!==============================================================================!
i=1
do while(i.le.Npoints)
  call random_number(r_i)
  r_i(:)=rmin(:)+r_i(:)*(rmax(:)-rmin(:))
  call normal_to_cartesian(r_i,x,.true.)
  call water_potential(x,V,forces)
  if(V.lt.E_cut)then
    r(:,i)=r_i(:)
    i=i+1
  endif
enddo
open(21,File='grid_ini.dat')
do i=1,Npoints
  write(21,*) r(:,i)
enddo
close(21)
!==============================================================================!
!   Compute pairwise energies for all the initial Grid Points U[x_{ij}]
!==============================================================================!
do i=2,Npoints
  do j=1,i-1
    Uij(i,j)=Pair_LJ_NRG(r(:,i),r(:,j),V)
    Uij(j,i)=Uij(i,j)
  enddo
enddo
!==============================================================================!
!               Generate QRG using MMC and forced minimization
!==============================================================================!
accept=0
counter=0
mv_cutoff=0.01
deltae1=0
do i=1,N_MMC_grid
  k=random_integer(1,Npoints)                               !Select Atom to Move
  call random_number(s)
  rr=r(:,k)+mv_cutoff*(2*s-1)               !random number (0,1), make it (-1,1)
  call normal_to_cartesian(rr,x,.true.)
  call water_potential(x,V,forces)
  if(V.lt.E_cut) then                            !Only consider if V(trial)<Ecut
    counter=counter+1
    U_move(k)=P_i(rr,V,x)
    Delta_E=0d0
    do j=1,Npoints
      if(j.ne.k) then
        U_move(j)=Pair_LJ_NRG(r(:,j),rr,V)
        Delta_E=Delta_E+Uij(j,k)-U_move(j)      !Energy change due to trial move
      endif
    enddo
    if(Delta_E.ge.0d0)then      !Forced minimization: accept if energy decreases
      Uij(:,k)=U_move(:)
      Uij(k,:)=U_move(:)
      accept=accept+1
      r(:,k)=rr(:)
      deltae1=deltae1+Delta_E
    endif
  endif
!for MMC want acceptance ~30-50%, adjust trial movement displacement accordingly
  if(mod(i,MMC_freq)==0)then
    if(dble(accept)/counter.lt.0.3)then
      mv_cutoff=mv_cutoff*0.9
    else
      mv_cutoff=mv_cutoff*1.1
    endif
  accept=0
  counter=0
  deltae1=0
  endif
enddo
open(22,File='grid.dat')
do i=1,Npoints
  write(22,*) r(:,i)
enddo
close(22)
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
write(99,*) 'Npoints ==> ',Npoints
write(99,*) 'c_LJ ==> ',c_LJ
write(99,*) 'N_MMC_box ==> ',N_MMC_box
write(99,*) 'N_MMC_grid ==> ',N_MMC_grid
write(99,*) 'MMC_freq ==> ',MMC_freq
write(99,*) 'Simulation Time==> ',time2-time1
close(99)
end program main_grid
