!=================20==================40==================60==================80
!                          Grid Generation Module
!==============================================================================!
!Generate a n-dimensional QRG, optimized using a quasi-Lennard Jones Potential
!Forced minimization; accept any trial moves that reduces the system's energy
!Use Metropolis Monte Carlo to determine the box size for evaluating the
!distribution.
!A (uniform) square grid is then used to normalize the distribution
!==============================================================================!
!       Modified:
!   29 April 2020
!       Author:
!   Shane Flynn
!==============================================================================!
module qrg_grid_mod
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
open(15,File='cage_tip4p.xyz')
read(15,*) Natoms
d=3*Natoms
close(15)
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
!x0(d)              ==>Coordinates (entire system)
!atom_type(Natoms)  ==>Element abbreviation from xyz
!==============================================================================!
implicit none
integer::d,Natoms,i
character(len=2)::atom_type(Natoms)
double precision::mass(Natoms),sqrt_mass(d),x0(d)
open(16,File='cage_tip4p.xyz')
read(16,*)
read(16,*)
do i=1,Natoms
  read(16,*) atom_type(i),x0(3*i-2:3*i)                 !input is xyz therefore 3
  mass(i)=Atom_mass(atom_type(i))
  sqrt_mass(3*i-2:3*i)=sqrt(mass(i))
enddo
close(16)
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
subroutine convert_to_atomic_units(potential,d,x0,E0,forces,E_cut,bohr,autocm,&
  autokcalmol)
!==============================================================================!
!Use atomic units for calculations (Assumes input geometry is in Angstroms)
!==============================================================================!
!potential          ==>potential name
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
double precision::x0(d),E0,forces(d),E_cut,bohr,autocm,autokcalmol
x0=x0/bohr
E_cut=E_cut/autokcalmol
open(17,File='cluster_energy.dat')
write(17,*) 'x0 in atomic units (xo/bohr)', x0
write(17,*) 'E_cut in atomic units (assumed kcal/mol input)', E_cut
call water_potential(potential,d,x0,E0,forces)
write(17,*) 'E0 (atomic) ==> ', E0
write(17,*) 'E0 (cm-1) ==> ', E0*autocm
write(17,*) 'E0 (kcal/mol) ==> ', E0*autokcalmol
close(17)
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
if(potential=='tip4p'.or.potential=='TIP4P') then
  call TIP4P(d/9,x,V,forces)
else
  stop 'Cannot Identify Potential, Check "potentials" Subroutine'
endif
end subroutine water_potential
!==============================================================================!
subroutine Get_Hessian(potential,d,d1,x0,Hess_Mat)
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
character(len=20)::potential
integer::d,d1,i,j
double precision::Hess_Mat(d1,d1),x0(d),x1(d),force0(d),force1(d),E0
double precision,parameter::s=1d-6
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
double precision::Hess_mat(d1,d1),omega(d1),U(d1,d1),work(max(1,lwork))
double precision::temp(d1),temp2(d1,d1)
U=Hess_mat
call dsyev('v','u',d1,U,d1,omega,work,Lwork,info)
open(18,File='freq_scaled_hess.dat')
write(*,*) 'Frequencies from the Mass-Scaled Hessian:'
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
subroutine normal_cartesian_potential(potential,d,d1,d2,x0,r_i,V,U,sqrt_mass)
!==============================================================================!
!Evaluate potential in d-dimensionality. Lmon moves 1 monomer of the cluster
!We only have dynamics over d2 subspace of the d1 monomer
!==============================================================================!
!x0(d)        ==>initial cluster configuration
!V            ==>Potential Energy Evalation
!U            ==>Normal mode eigenvectors
!r(d2) coordinate to evaluate the potential at
!r(d1) need monomer size to scale coordinates
!x scaled coordinate in cartesian space to call potential with
!==============================================================================!
implicit none
character(len=20)::potential
integer::d,d1,d2
double precision::x0(d),x(d),r_i(d2),rr(d1),U(d1,d1),forces(d),sqrt_mass(d),V
rr=0        !need d1-dimensional coordinate for monomer, (1:d2)=r d(2+1,d1)=0
rr(1:d2)=r_i(1:d2)
x=x0
x(1:d1)=x(1:d1)+matmul(U,rr)/sqrt_mass(1:d1)
call water_potential(potential,d,x,V,forces)
end subroutine normal_cartesian_potential
!==============================================================================!
function P_i(potential,d,d1,d2,x0,r_i,U,sqrt_mass,V,E_cut,integral_P)
!==============================================================================!
!Target Distribution Function, !defined according to a semi-classical argument:
!B. Poirier, “Algebraically self-consistent quasiclassical approximation on
!phase space,” Found. Phys. 30, 1191–1226 (2000).
!==============================================================================!
!potential          ==>Potential name
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!d1                 ==>subspace dimensionality
!x(d)               ==>i-th grid points coordinates (x^i=x_1,x_2,...x_d)
!V                  ==>Potential Energy evaluation V(x_i)
!E_cut              ==>Distribution cutoff contour
!integral_P         ==>Normalization constant for the distribtion P(x)
!P_i                ==>evaluate P(x)
!==============================================================================!
implicit none
character(len=20)::potential
integer::d,d1,d2
double precision::x0(d),r_i(d2),U(d1,d1),sqrt_mass(d),E_cut,integral_P,V,P_i
!==============================================================================!
call normal_cartesian_potential(potential,d,d1,d2,x0,r_i,V,U,sqrt_mass)
if(V<E_cut) P_i=(E_cut-V)**(d2/2.)/integral_P
if(V>=E_cut) P_i=1d-20                        !Define distribution=0 beyond Ecut
end function P_i
!==============================================================================!
function random_integer(Nmin,Nmax)
!==============================================================================!
!Randomly generate an integer in the range Nmin-Nmax
!==============================================================================!
!Nmin           ==>minimum index value
!Nmax           ==>maximum index value
!a              ==>uniform pseudo-random number
!random_integer ==>integer returned
!==============================================================================!
implicit none
integer::Nmin,Nmax,random_integer
double precision::a
!==============================================================================!
call random_number(a)
random_integer=floor(a*(Nmax-Nmin+1))+Nmin
end function random_integer
!==============================================================================!
subroutine box_size_P(potential,d,d1,d2,x0,U,sqrt_mass,E_cut,rmin,rmax,&
  N_MMC_box,integral_P)
!==============================================================================!
!Determine the box size for normalizing P(x) using Metropolis Monte Carlo
!This subroutine Computes rmin(d2) and rmax(d2)
!==============================================================================!
!potential          ==>Potential name
!d2                 ==>Monomer Subspace dimensionality (x^i=x_1,x_2,...x_d2)
!V                  ==>Potential Energy evaluation V(x)
!E_cut              ==>Distribution cutoff contour
!rmin(d2)           ==>Minimum of normalization box size
!rmax(d2)           ==>Maximum of normalization box size
!N_MMC_box          ==>Number of MMC Iterations to determine box-size
!integral_P         ==>Normalization constant for the distribtion P(x)
!==============================================================================!
double precision,parameter::mv_cutoff=0.1
character(len=20)::potential
integer::d,d1,d2,N_MMC_box,i,j
double precision::x0(d),U(d1,d1),sqrt_mass(d),E_Cut,rmin(d2),rmax(d2),integral_P
double precision::r_i(d2),s(d2),r_trial(d2),dummy,V
!==============================================================================!
integral_P=1d0                                     !Initially set to 1 to call P
r_i=0d0
rmin=r_i
rmax=r_i
do i=1,N_MMC_box
  call random_number(s)
!trial move needs to be (-1,1), random numbers are (0,1) ==>s=2*s-1
  r_trial =r_i+mv_cutoff*(2*s-1)
  call random_number(dummy)                             !MMC acceptance criteria
  if(P_i(potential,d,d1,d2,x0,r_trial,U,sqrt_mass,V,E_cut,integral_P)/&
  P_i(potential,d,d1,d2,x0,r_i,U,sqrt_mass,V,E_cut,integral_P).ge.dummy) then
    r_i=r_trial
    do j=1,d2
      if(rmin(j)>r_i(j)) rmin(j)=r_i(j)
      if(rmax(j)<r_i(j)) rmax(j)=r_i(j)
    enddo
  endif
enddo
end subroutine box_size_P
!==============================================================================!
subroutine compute_integral_P(potential,d,d1,d2,x0,V,U,sqrt_mass,E_cut,&
  integral_P,rmin,rmax,N_1D)
!==============================================================================!
!Use a (uniform) square grid to integrate P(x)
!Box size for the grid is determined in the box_size subroutine
!int P(x)~Area_Square/N sum_n=1,N P(x_n)
!This subroutine Computes integral_P
!==============================================================================!
!potential          ==>Potential name
!n                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_n)
!V                  ==>Potential Energy evaluation V(x)
!E_cut              ==>Distribution cutoff contour
!xmin(n)            ==>Minimum of normalization box size
!xmax(n)            ==>Maximum of normalization box size
!N_1D               ==>Number of points in 1 dimension for computing integral_P
!integral_P         ==>Normalization constant for the distribtion P(x)
!Ntotal             ==>Total number of evaluations for all dimensions
!Moment             ==>First Moment for the distribution
!==============================================================================!
character(len=20)::potential
integer::d,d1,d2,N_1D,i,j,index1(d2),Ntotal
double precision::x0(d),r_i(d2),V,U(d1,d1),sqrt_mass(d),E_cut,integral_P,Moment
double precision::delr(d2),rmin(d2),rmax(d2),dummy
!==============================================================================!
open(19,File='direct_grid.dat')
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
  dummy=P_i(potential,d,d1,d2,x0,r_i,U,sqrt_mass,V,E_cut,integral_P)
  Moment=Moment+dummy
  if(V<E_cut) write(19,*) r_i
enddo
dummy=1./N_1D**d2
do j=1,d2
  dummy=dummy*(rmax(j)-rmin(j))
enddo
integral_P=dummy*Moment
close(19)
end subroutine compute_integral_P
!==============================================================================!
subroutine initial_distribution(potential,d,d1,d2,x0,r,V,U,sqrt_mass,Npoints,&
  E_cut,rmin,rmax)
!==============================================================================!
!Generate an initial distribution to optimize with the Quasi-Regular algorithm
!Accepts any pseudo-random number generated within the cutoff contour E_cut
!==============================================================================!
!potential          ==>Potential name
!Npoints            ==>Number of points to generate
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!r(d,Npoints)       ==>All grid points coordinates (x^i=x_1,x_2,...x_d)
!V                  ==>Potential Energy evaluation V(x_i)
!E_cut              ==>Distribution cutoff contour
!rmin(d)            ==>Minimum of normalization box size
!rmax(d)            ==>Maximum of normalization box size
!==============================================================================!
implicit none
character(len=20)::potential
integer::Npoints,d,d1,d2,i
double precision::x0(d),U(d1,d1),sqrt_mass(d),E_Cut,r(d2,Npoints),r_i(d2)
double precision::rmin(d2),rmax(d2),V
!==============================================================================!
i=1
do while(i.le.Npoints)
  call random_number(r_i)
  r_i(:)=rmin(:)+r_i(:)*(rmax(:)-rmin(:))
  call normal_cartesian_potential(potential,d,d1,d2,x0,r_i,V,U,sqrt_mass)
  if(V<E_cut)then
    r(:,i)=r_i(:)
    i=i+1
  endif
enddo
open(20,File='grid_ini.dat')
do i=1,Npoints
  write(20,*) r(:,i)
enddo
close(20)
end subroutine initial_distribution
!==============================================================================!
function Pair_LJ_NRG(potential,d,d1,d2,x0,r_i,r_j,V,U,sqrt_mass,E_cut,&
  integral_P,c_LJ,Npoints)
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
character(len=20)::potential
integer::d,d1,d2,Npoints
double precision::x0(d),r_i(d2),r_j(d2),V,U(d1,d1),sqrt_mass(d),E_cut,integral_P
double precision::c_LJ,a,b,Pair_LJ_NRG,sigma1,sigma2
!==============================================================================!
a=sum((r_i(:)-r_j(:))**2)
sigma1=c_LJ*(P_i(potential,d,d1,d2,x0,r_i,U,sqrt_mass,V,E_cut,integral_P)*&
        Npoints)**(-1./d2)
sigma2=c_LJ*(P_i(potential,d,d1,d2,x0,r_j,U,sqrt_mass,V,E_cut,integral_P)*&
        Npoints)**(-1./d2)
b=(sigma2**2/a)
a=(sigma1**2/a)
Pair_LJ_NRG=a**(d2+9)-a**(d2+3)+b**(d2+9)-b**(d2+3)
end function Pair_LJ_NRG
!==============================================================================!
subroutine initial_pair_energy(potential,d,d1,d2,x0,V,U,sqrt_mass,E_cut,&
  integral_P,c_LJ,Npoints,r,Uij)
!==============================================================================!
!Compute the pairwise energies for all the initial Grid Points U[x_{ij}]
!==============================================================================!
implicit none
character(len=20)::potential
integer::d,d1,d2,Npoints,i,j
double precision::x0(d),U(d1,d1),sqrt_mass(d),E_cut,c_LJ,Uij(Npoints,Npoints)
double precision::r(d2,Npoints),integral_P,V
!==============================================================================!
!potential            ==>Potential name
!Npoints              ==>Number of points to generate
!d                    ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!r(d2,Npoints)         ==>All grid points coordinates (x^i=x_1,x_2,...x_d)
!Uij(Npoints,Npoints) ==>Pairwise energies for all i-j points
!V_i                  ==>Potential Energy evaluation V(x_i)
!E_cut                ==>Distribution cutoff contour
!c_LJ                 ==>Parameter for q-LJ pseudo-potential
!integral_P           ==>Normalization constant for the distribtion P(x)
!==============================================================================!
do i=2,Npoints
  do j=1,i-1
    Uij(i,j)=Pair_LJ_NRG(potential,d,d1,d2,x0,r(:,i),r(:,j),V,U,sqrt_mass,&
    E_cut,integral_P,c_LJ,Npoints)
    Uij(j,i)=Uij(i,j)
  enddo
enddo
end subroutine initial_pair_energy
!==============================================================================!
subroutine Quasi_Regular(potential,d,d1,d2,x0,V,U,sqrt_mass,E_cut,integral_P,r,&
  Uij,c_LJ,N_MMC_grid,MMC_freq,Npoints)
!==============================================================================!
!Generate grid using a Quasi-Regular sequence of points
!Grid is generated using MMC and forced minimization using our pseudo-potential
!==============================================================================!
!potential          ==>Potential name
!Npoints            ==>Number of points to generate
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x(d,Npoints)       ==>All grid points coordinates (x^i=x_1,x_2,...x_d)
!U(Npoints,Npoints) ==>Pairwise energies for all i-j points
!V_i                ==>Potential Energy evaluation V(x_i)
!E_cut              ==>Distribution cutoff contour
!c_LJ               ==>Parameter for q-LJ pseudo-potential
!N_MMC_grid         ==>Number of MMC Iterations to optimize QRG
!MMC_freq           ==>Frequency to update QRG MMC grid mv_cutoff
!integral_P         ==>Normalization constant for the distribtion P(x)
!==============================================================================!
implicit none
character(len=20)::potential
integer::d,d1,d2,Npoints,N_MMC_grid,MMC_freq,accept,counter,i,j,k
double precision::x0(d),U(d1,d1),sqrt_mass(d),r(d2,Npoints),Uij(Npoints,Npoints)
double precision::c_LJ,E_cut,mv_cutoff,deltae1,Delta_E,U_move(Npoints),s(d2)
double precision::rr(d2),integral_P,V
!==============================================================================!
accept=0
counter=0
mv_cutoff=0.01
deltae1=0
do i=1,N_MMC_grid
  k=random_integer(1,Npoints)                               !Select Atom to Move
  call random_number(s)
  rr=r(:,k)+mv_cutoff*(2*s-1)   !random number (0,1), make it (-1,1) ==> s=2*s-1
!==============================================================================!
  call normal_cartesian_potential(potential,d,d1,d2,x0,rr,V,U,sqrt_mass)
  if(V.lt.E_cut) then                          !Only consider if V(trial)<Ecut
    counter=counter+1
    U_move(k)=P_i(potential,d,d1,d2,x0,rr,U,sqrt_mass,V,E_cut,integral_P)
    Delta_E=0d0
    do j=1,Npoints
      if(j.ne.k) then
        U_move(j)=Pair_LJ_NRG(potential,d,d2,d2,x0,r(:,j),rr,V,U,sqrt_mass,&
        E_cut,integral_P,c_LJ,Npoints)
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
    if(dble(accept)/counter<0.3)then
      mv_cutoff=mv_cutoff*0.9
    else
      mv_cutoff=mv_cutoff*1.1
    endif
  accept=0
  counter=0
  deltae1=0
  endif
enddo
open(21,File='grid.dat')
do i=1,Npoints
  write(21,*) r(:,i)
enddo
close(21)
!==============================================================================!
end subroutine Quasi_Regular
!==============================================================================!
subroutine write_out(potential,Npoints,d,d1,d2,E_cut,rmin,rmax,N_1D,N_MMC_box,&
  c_LJ,N_MMC_grid,MMC_freq,integral_P,time1,time2)
!==============================================================================!
!Write output file
!==============================================================================!
!potential          ==>Potential name
!Npoints            ==>Number of points to generate
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!d1                 ==>Subspace dimensionality (x^i=x_1,x_2,...x_d)
!E_cut              ==>Distribution cutoff contour
!xmin(d)            ==>Minimum of normalization box size
!xmax(d)            ==>Maximum of normalization box size
!N_1D               ==>Number of points in 1 dimension for computing integral_P
!N_MMC_box          ==>Number of MMC Iterations to determine box-size
!c_LJ               ==>Parameter for q-LJ pseudo-potential
!N_MMC_grid         ==>Number of MMC Iterations to optimize QRG
!MMC_freq           ==>Frequency to update QRG MMC grid mv_cutoff
!integral_P         ==>Normalization constant for the distribtion P(x)
!==============================================================================!
implicit none
integer::Npoints,d,d1,d2,N_1D,N_MMC_grid, MMC_freq,i,N_MMC_box
double precision::E_cut,c_LJ,integral_P,rmin(d2),rmax(d2),time1,time2
character(len=20)::potential
!==============================================================================!
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
end subroutine write_out
!==============================================================================!
end module qrg_grid_mod
