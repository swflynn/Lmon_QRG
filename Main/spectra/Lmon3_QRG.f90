!=============================================================================80
!                        Lmon-3 QRG for water clusters                         !
!==============================================================================!
!Things to do:
!starting with old code as a template, need to make this code work for our new
!project.
!==============================================================================!
!DGB-QRG analysis for water cluster using the local monomer approximation
!Normal mode coordinates for analysis, Integration with quadrature
!Requires llapack (linear algebra package) for generalized eigenvalue problem
!==============================================================================!
!       Modified:
!16 December 2019
!       Author:
!Shane Flynn
!==============================================================================!
module Lmon_QRG
implicit none
!==============================================================================!
!                            Global Paramaters
!d   ==> Monomer Dimensionality
!d1  ==> Dimensionality of Monomer Subspace
!==============================================================================!
double precision,parameter::bohr=0.52917721092
double precision,parameter::autocm=2.194746313D5
double precision,parameter::autokcalmol=627.5096
double precision,parameter::melectron=1822.88839
double precision,parameter::Hmass=1.00782503223*melectron
double precision,parameter::Omass=15.99491461957*melectron
integer,parameter::d1=9
integer,parameter::d2=3
!==============================================================================!
!                            Global Variables                                  !
!==============================================================================!
!Natoms       ==>Number of Atoms
!d            ==>Total System Dimensionality
!==============================================================================!
integer::Natoms,d
character(len=2),allocatable,dimension(:)::atom_type
double precision,allocatable,dimension(:)::mass(:),sqrt_mass
character(len=5)::potential
!==============================================================================!
contains
!==============================================================================!
function Atom_Mass(atom)
!Numerically computed Hessian
!s          ==> Perturbation Parameter
!H	    ==> (d1,d1); Symmetrized Mass-Scaled Hessian
!q          ==> (d); XYZ Configuration
!force	    ==> (d); Forces from Water Potential
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
  STOP 'Check Atom_Mass Function'
endif
end function Atom_Mass
!==============================================================================!
subroutine Hessian(q,H)
!==============================================================================!
!Numerically computed Hessian
!s          ==> Perturbation Parameter
!H	    ==> (d1,d1); Symmetrized Mass-Scaled Hessian
!q          ==> (d); XYZ Configuration
!force	    ==> (d); Forces from Water Potential
!==============================================================================!
!!!implicit none
!!!integer::i,j
!!!double precision::H(d1,d1),q(d),r(d),E,force(d),force0(d)
!!!double precision,parameter::s=1d-6
!!!r=q
!!!call water_potential(Natoms/3,r,E,force0)
!!!do i=1,d1
!!!   r(i)=q(i)+s
!!!   call water_potential(Natoms/3,r,E,force)
!!!   r(i)=q(i)
!!!   do j=1,d1
!!!      H(i,j)=(force0(j)-force(j))/s
!!!   enddo
!!!enddo
!==============================================================================!
!!!!                   Symmetrize and Mass-Scale the Hessian
!==============================================================================!
!!!do i=1,d1
!!!   do j=1,i
!!!      if(i.ne.j) H(i,j)=(H(i,j)+H(j,i))/2.
!!!      H(i,j)=H(i,j)/(sqrt_mass(i)*sqrt_mass(j))
!!!      if(i.ne.j) H(j,i)=H(i,j)
!!!   enddo
!!!enddo
!!!!write(*,*) 'Mass-Scaled Hessian', H
!!!end subroutine Hessian
!==============================================================================!
!!!subroutine Freq_Hess(d,Hess,omega,U)
!==============================================================================!
!!!!Compute Eigenvalues and Eigenvectors of Hessian
!!!!Uses the LLAPACK real symmetric eigen-solver (dsygev)
!!!!Hess   ==> (d,d); Hessian Matrix
!!!!omega  ==> (d); Eigenvalues of the Hessian
!!!!U      ==> (d,d); Eigenvectors of the Hessian
!!!!       LLAPACK (dsyev):
!!!!v      ==> Compute both Eigenvalues and Eigenvectors
!!!!u      ==> Use Upper-Triangle of matrix
!==============================================================================!
!!!implicit none
!!!integer::i,info,lwork,d1
!!!double precision::Hess(d1,d1),omega(d1),U(d1,d1)
!!!double precision,allocatable,dimension(:)::work
!!!lwork=max(1,3*d1-1)                      !suggested by LAPACK Developers
!!!allocate(work(max(1,lwork)))            !suggested by LAPACK Developers
!!!U=Hess
!!!call dsyev('v','u',d1,U,d1,omega,work,lwork,info)
!!!write(*,*) 'Frequencies from the Hessian:'
!!!do i=d1,1,-1
!!!    omega(i)=sign(sqrt(abs(omega(i))),omega(i))
!!!    write(*,*) omega(i), omega(i)*autocm, 'normalized==1?', sum(U(:,i)**2)
!!!enddo
!!!end subroutine Freq_Hess
!==============================================================================!
!!!subroutine water_potential(NO,q,energy,force)
!!!use iso_c_binding
!!!use TIP4P_module
!==============================================================================!
!!!!Potential Energy module for water
!!!!       Variables:
!!!!NO         ==> number of water molecules
!!!!q          ==> (9*NO); coordinates
!!!!force      ==> (9*NO) Forces computed in external water potential
!!!!energy     ==> Potential Energy
!==============================================================================!
!!!implicit none
!!!integer,intent(in)::NO
!!!double precision,dimension(9*NO),intent(in)::q
!!!double precision,dimension(9*NO),intent(inout)::force
!!!double precision,intent(inout)::energy
!!!if(potential=='tip4p') then
!!!   call TIP4P(NO,q,energy,force)
!!!elseif(potential=='mbpol') then
!!!    call calcpotg(NO,energy,q*bohr,force)
!!!    force=-force*bohr/autokcalmol
!!!    energy=energy/autokcalmol
!!!else
!!!   stop 'Cannot Identify Potential, Check "Water_Potential" Subroutine'
!!!endif
!!!end subroutine water_potential
!==============================================================================!
end module Lmon_QRG
!==============================================================================!
program main
use Lmon_QRG
!==============================================================================!
!==============================================================================!
!coord_in       ==> Input Geometry File-Name (xyz file) !!!!NG             ==> Number of Gaussians
!Nsobol         ==> Number of Sobol Points for numerical integration
!ii,skip        ==> Integer(kind=8): necessary for using Sobol Module
!alpha          ==> (NG): Width of Gaussian
!q0             ==> (d): Input Geometry coordinates
!r              ==> (d,NG): Gaussian Centers
!r2             ==> (d): ith Gaussian coordinate for integration
!force          ==> (d): Forces. See "Toy_Force"
!Hess           ==> (d,d): Mass-Scaled Hessian. See "Toy_Hessian"
!omega          ==> (d): Frequencies (Eigenvalues). See "Freq_Hess"
!U              ==> (d,d): Hess Eigenvectors. See "Freq_Hess"
!z              ==> (d,NG): Quasi-Random Sequence for Gaussian Centers
!z2             ==> (d,Nsobol): Quasi-Random Sequence for Integration
!Smat           ==> (NG,NG): Overlap Matrix
!Vmat           ==> (NG,NG): Potential Energy Matrix
!Hmat           ==> (NG,NG): Hamiltonian Matrix (V+T) for eigenvalue problem
!eigenvalues    ==> (NG): Eigenvalues of the Hamiltonian Matrix
!lambda         ==> (NG): Eigenvalues of the Overlap Matrix
!r_ij           ==> center of the i,j matrix element
!E0             ==> Energy Evaluation at the minimum configuration
!pot_ene        ==> Potential Energy evaluated in q space
!d2              ==> subspace dimensionality
!mon_num        ==> monomer dynamics are being run on
!==============================================================================!
implicit none
character(len=50)::coord_in
integer::NG,Nsobol,mon_num
integer::i,j,l,m,n,counter,data_freq
integer*8::ii,skip,skip2
double precision::E0,pot_ene,s_sum,kinetic,alpha0,t1,t2
!==============================================================================!
!!!double precision,allocatable,dimension(:)::q0,q1,force,r2,eigenvalues,temp,r_ij
!!!double precision,allocatable,dimension(:)::omega,alpha,lambda
!!!double precision,allocatable,dimension(:,:)::r,Hess,U,z,z2,Hmat,Smat,Vmat,temp2
!==============================================================================!
!!!!                       LLAPACK dsygv variables                                !
!==============================================================================!
!!!integer::itype,info,lwork
!!!double precision,allocatable,dimension(:)::work
!==============================================================================!
!!!!                           Read Input Data File                               !
!==============================================================================!
!!!call cpu_time(t1)
!!!read(*,*) potential
!!!read(*,*) alpha0
!!!read(*,*) NG
!!!read(*,*) Nsobol
!!!read(*,*) data_freq
!!!read(*,*) coord_in
!!!read(*,*) mon_num
!!!skip=NG
!!!skip2=Nsobol
!!!write(*,*) 'Test 1; Successfully Read Input Data File'
!==============================================================================!
!!!!                         Set Input Water Geometry
!==============================================================================!
!!!open(16,File=coord_in)
!!!read(16,*) Natoms
!!!read(16,*)
!!!d=3*Natoms
!==============================================================================!
!!!allocate(atom_type(Natoms),mass(Natoms),sqrt_mass(d),q0(d),q1(d))
!!!allocate(force(d),Hess(d1,d1),omega(d1),U(d1,d1),z(d2,NG),r_ij(d1))
!!!allocate(z2(d1,Nsobol),r(d1,NG),alpha(NG),Smat(NG,NG),eigenvalues(NG))
!!!allocate(Vmat(NG,NG),Hmat(NG,NG),r2(d1),lambda(NG),temp(d1),temp2(d1,d1))
!==============================================================================!
!!!!                         Input Configuration Energy
!==============================================================================!
!!!do i=1,Natoms
!!!    read(16,*) atom_type(i),q0(3*i-2:3*i)
!!!    mass(i)=Atom_mass(atom_type(i))
!!!    sqrt_mass(3*i-2:3*i)=sqrt(mass(i))
!!!enddo
!!!close(16)
!!!write(*,*) 'Dimen ==> ', Dimen, 'd ==> ', d, 'd1 ==> ', d1, 'd2 ==> ', d2
!==============================================================================!
!!!! 		Convert to Atomic Units for Calculations
!==============================================================================!
!!!q0=q0/bohr
!!!call water_potential(Natoms/3,q0,E0,force)
!!!write(*,*) 'E0 ==> ',E0*autocm, 'cm-1', E0*autokcalmol, 'kcal/mol'
!==============================================================================!
!!!! 			Compute Hessian and Frequencies
!==============================================================================!
!!!call Hessian(q0,Hess)
!!!call Freq_Hess(d,Hess,omega,U)
!!!write(*,*) 'Test 2; Successfully Computed Hessian and Frequencies'
!==============================================================================!
!!!!                       Subspace Needs Largest Eigenvalues
!!!!llapack outputs smallest to largest, need to reorder for subspace
!==============================================================================!
!!!temp=omega
!!!temp2=U
!!!do i=1,d1
!!!    omega(i)=temp(d1+1-i)
!!!    U(:,i)=temp2(:,d1+1-i)
!!!enddo
!==============================================================================!
!!!!           Generate Gaussian Centers with a quasi-random grid
!!!!Assume centers are in normal-mode space r (not coordinate space q)
!==============================================================================!
!!!r=0d0
!!!do ii=1,NG
!!!    call sobol_stdnormal(d2,skip,z(:,ii))
!!!    do i=1,d2
!!!        r(i,ii)=z(i,ii)/sqrt(abs(omega(i)))
!!!    enddo
!!!enddo
!!!write(*,*) 'Test 3; Successfully Generated Gaussian Centers'
!==============================================================================!
!!!!                       Generate Alpha Scaling
!==============================================================================!
!!!do i=1,NG
!!!    alpha(i)=alpha0*NG**(2./d2)*exp(-1./d2*sum(r(:,i)**2*omega(:)))
!!!enddo
!!!open(unit=17,file='centers.dat')
!!!do i=1,NG
!!!    write(17,*) alpha(i), r(:,i)
!!!enddo
!!!close(17)
!!!write(*,*) 'Test 4; Successfully Generated Gaussian Widths'
!==============================================================================!
!!!!                           Overlap Matrix (S)
!==============================================================================!
!!!do i=1,NG
!!!    do j=i,NG
!!!        s_sum=sum(omega(:)*(r(:,i)-r(:,j))**2)
!!!        Smat(i,j)=sqrt(alpha(i)+alpha(j))**(-d2)&
!!!                 *exp(-0.5*alpha(i)*alpha(j)/(alpha(i)+alpha(j))*s_sum)
!!!        Smat(j,i)=Smat(i,j)
!!!    enddo
!!!enddo
!==============================================================================!
!!!!                   Check to see if S is positive definite
!!!! If this is removed, you need to allocate llapack arrays before Hamiltonian
!==============================================================================!
!!!lwork=max(1,3*NG-1)
!!!allocate(work(max(1,lwork)))
!!!call dsyev('v','u',NG,Smat,NG,lambda,work,Lwork,info)
!!!write(*,*) 'Info (Overlap Matrix) ==> ', info
!!!open(unit=18,file='overlap_eigenvalues.dat')
!!!do i=1,NG
!!!    write(18,*) lambda(i)
!!!enddo
!!!close(18)
!!!write(*,*) 'Test 5; Successfully Computed Overlap Matrix'
!==============================================================================!
!!!!                   Generate Sequence For Evaluating Potential
!==============================================================================!
!!!do ii=1,Nsobol
!!!    call sobol_stdnormal(d1,skip2,z2(:,ii))
!!!enddo
!!!write(*,*) 'Test 6; Successfully Generated Integration Sequence'
!==============================================================================!
!!!!                              Evaluate Potential
!==============================================================================!
!!!Vmat=0d0
!!!counter=0
!!!open(unit=19,file='eigenvalues.dat')
!!!open(unit=20,file='fundamentals.dat')
!!!do while(counter.lt.Nsobol/data_freq)
!!!    do i=1,NG
!!!       do j=i,NG
!!!          r_ij(:)=(alpha(i)*r(:,i)+alpha(j)*r(:,j))/(alpha(i)+alpha(j))
!!!            do l=1+(counter*data_freq),(counter+1)*data_freq
!!!                r2(:)=r_ij(:)+z2(:,l)/sqrt(omega(:)*(alpha(i)+alpha(j)))
!==============================================================================!
!!!!                   Evaluate Potential
!!!!Use initial configuration q0 to keep dimensions correct
!!!!Only modify subspace coordinates undergoing dynamics, to cartesian space (q1)
!==============================================================================!
!!!                q1=q0
!!!                q1(1:d1)=q1(1:d1)+matmul(U,r2)/sqrt_mass(1:d1)
!!!                call water_potential(Natoms/3,q1,pot_ene,force)
!!!                Vmat(i,j)=Vmat(i,j)+pot_ene
!!!            enddo
!!!        enddo
!!!    enddo
!==============================================================================!
!!!!                     Solve Generalized Eigenvalue Problem
!==============================================================================!
!!!    write(*,*) 'Iteration ==> ', (counter+1)*data_freq
!!!    Hmat=Vmat
!!!    do m=1,NG
!!!        do n=m,NG
!!!            s_sum=sum(omega(:)*(r(:,m)-r(:,n))**2)
!!!            Smat(m,n)=sqrt(alpha(m)+alpha(n))**(-d2)&
!!!                *exp(-0.5*alpha(m)*alpha(n)/(alpha(m)+alpha(n))*s_sum)
!!!            Smat(n,m)=Smat(m,n)
!!!            kinetic=Smat(m,n)*0.5*alpha(m)*alpha(n)/(alpha(m)+alpha(n))&
!!!           *sum(omega(:)-(alpha(m)*alpha(n)*(omega(:)**2*(r(:,m)-r(:,n))**2)&
!!!                /(alpha(m)+alpha(n))))
!!!            Hmat(m,n)=Hmat(m,n)*(Smat(m,n)/((counter+1)*data_freq))+kinetic
!!!            Hmat(n,m)=Hmat(m,n)
!!!        enddo
!!!    enddo
!!!    itype=1
!!!    eigenvalues=0d0
!!!    call dsygv(itype,'n','u',NG,Hmat,NG,Smat,NG,eigenvalues,work,Lwork,info)
!!!    write(19,*) (counter+1)*data_freq, eigenvalues(:)*autocm
!!!    write(20,*) (counter+1)*data_freq, eigenvalues(1)*autocm, &
!!!        eigenvalues(2)*autocm, (eigenvalues(2)-eigenvalues(1))*autocm, &
!!!        eigenvalues(3)*autocm, (eigenvalues(3)-eigenvalues(1))*autocm, &
!!!        eigenvalues(4)*autocm, (eigenvalues(4)-eigenvalues(1))*autocm
!!!    write(*,*) 'info ==> ', info
!!!    counter=counter+1
!!!enddo
!!!close(19)
!!!close(20)
!!!open(unit=21,file='frequency.dat')
!!!write(21,*) mon_num, (eigenvalues(2)-eigenvalues(1))*autocm
!!!write(21,*) mon_num, (eigenvalues(3)-eigenvalues(1))*autocm
!!!write(21,*) mon_num, (eigenvalues(4)-eigenvalues(1))*autocm
!!!close(21)
!!!open(unit=22,file='harmonic.dat')
!!!write(22,*) mon_num, omega(1)*autocm
!!!write(22,*) mon_num, omega(2)*autocm
!!!write(22,*) mon_num, omega(3)*autocm
!!!close(22)
!!!open(unit=23,file='alpha_conv.dat')
!!!write(23,*) alpha0, eigenvalues(:)*autocm
!!!close(23)
!!!write(*,*) 'Hello Universe!'
!!!call cpu_time(t2)
!==============================================================================!
!!!!                               output file                                    !
!==============================================================================!
!!!open(90,file='simulation.dat')
!!!write(90,*) 'Natoms ==> ', Natoms
!!!write(90,*) 'Dimensionality Input System ==> ', d
!!!write(90,*) 'Monomer Dimensionality ==> ', d1
!!!write(90,*) 'Subspace Dimensionality ==> ', d2
!!!write(90,*) 'potential ==> ', potential
!!!write(90,*) 'alpha0 ==> ', alpha0
!!!write(90,*) 'N_gauss ==> ', NG
!!!write(90,*) 'N_Sobol ==> ', Nsobol
!!!write(90,*) 'omega ==> ', omega
!!!write(90,*) 'Total Time ==> ', t2-t1
!!!close(90)
!!!end program main
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!morse 3D code for computing eigenspectra is pasted below, will need to merge
!!!!these two codes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=============================================================================80
!!!!                      3D Morse Code with QRG and GHQ
!==============================================================================!
!!!!Quasi-Regular Grid for computing vibrational spectra
!!!!generate alpha based on nearest neghbor
!!!!see distribtuion code for alpha as a funciton of the target distribution
!!!!This code needs the gen_hermite_rule.f90 code (Gauss-Hermite-Quadriture)
!==============================================================================!
!!!!       Modified:
!!!!   15 May 2019
!!!!       Author:
!!!!   Shane Flynn
!==============================================================================!
!!!module QRGB_mod
!!!implicit none
!==============================================================================!
!!!!                            Global Variables
!==============================================================================!
!!!!d              ==>i-th gaussian dsionality (x^i=x^i_1,x^i_2,..,x^i_d)
!==============================================================================!
!!!integer::d
!==============================================================================!
!!!contains
!==============================================================================!
!!!function V(x)
!==============================================================================!
!!!!Hard-coded Morse Potential Energy
!==============================================================================!
!!!!x              ==>(d) ith particles coordinate x^i_1,..,x^i_d
!!!!V              ==>evaluate V(x)
!!!!D_morse        ==>Parameter for Morse Potential
!!!!omega          ==>(d) Parameter for Morse Potential
!==============================================================================!
!!!implicit none
!!!double precision::x(d),V
!!!double precision,parameter::omega(3)=(/0.2041241,0.18371169,0.16329928/)
!!!double precision,parameter::D_morse=12.
!!!V=D_morse*sum((exp(-omega(:)*x(:))-1.)**2 )
!!!end function V
!==============================================================================!
!==============================================================================!
!!!end module QRGB_mod
!==============================================================================!
!==============================================================================!
!!!program main
!!!use QRGB_mod
!==============================================================================!
!!!!unif_grid      ==>If true alpha=constant, else alpha=alpha(distribution)
!!!!grid_in        ==>Filename Containing Gridpoints, see qlj_morse2d_grid.f90
!!!!theory_in      ==>Filename Containing Analytic Eigenvalues, see theory.f90
!!!!theory         ==>(NG) Analytic Eigenvalues
!!!!NG             ==>Number of Gaussian Basis Functions (gridpoints)
!!!!GH_order       ==>Number of Points for evaluating the potential (Gauss-Hermite)
!!!!alpha0         ==>Flat Scaling Parameter for Gaussian Widths
!!!!alpha          ==>(d) Gaussian Widths
!!!!RCN            ==>Recriprical Convergence Number, stability of Overlap Matrix
!!!!x              ==>(d) ith atoms coordinates
!!!!x_ij           ==>i-jth gaussian center (product of Gaussians is a Gaussian)
!!!!Smat           ==>(NG,NG) Overlap Matrix
!!!!Hmat           ==>(NG,NG) Hamiltonian Matrix
!==============================================================================!
!!!implicit none
!!!logical::unif_grid
!!!character(len=50)::grid_in,theory_in
!!!integer::NG,GH_order,i,j,l1,l2,l3
!!!double precision,parameter::pi=4.*atan(1d0)
!!!double precision::aij,r2,alpha0,Vij,RCN
!!!double precision,allocatable,dimension(:)::alpha,eigenvalues,x_ij,z,w,rr,theory
!!!double precision,allocatable,dimension(:,:)::x,Smat,Hmat
!==============================================================================!
!!!!                       LLAPACK dsygv variables                                !
!==============================================================================!
!!!integer::itype,info,lwork
!!!double precision,allocatable,dimension(:)::work
!==============================================================================!
!!!!                           Read Input Data File                               !
!==============================================================================!
!!!read(*,*) d
!!!read(*,*) NG
!!!read(*,*) GH_order
!!!read(*,*) grid_in
!!!read(*,*) theory_in
!!!read(*,*) alpha0
!==============================================================================!
!!!!                               Allocations
!==============================================================================!
!!!allocate(x(d,NG),x_ij(d),rr(d),alpha(NG),eigenvalues(NG),Smat(NG,NG))
!!!allocate(Hmat(NG,NG),theory(NG),z(GH_order),w(GH_order))
!!!write(*,*) 'Test 0; Successfully Read Input File'
!==============================================================================!
!!!!                           Read GridPoints x(d,NG)
!==============================================================================!
!!!open(17,File=grid_in)
!!!do i=1,NG
!!!    read(17,*) x(:,i)
!!!enddo
!!!close(17)
!==============================================================================!
!!!!                       Generate Alphas alpha(NG)
!!!!           Use nearest neighbor to define width alpha(NG)
!==============================================================================!
!!!do i=1,NG
!!!    alpha(i)=1d20                            !large distance for placeholder
!!!    do j=1,NG
!!!        if(j.ne.i) then
!!!            r2=sum((x(:,i)-x(:,j))**2)          !distance between gridpoints
!!!            if(r2<alpha(i)) alpha(i)=r2
!!!        endif
!!!    enddo
!!!    alpha(i)=alpha0/alpha(i)
!!!enddo
!!!write(*,*) 'Test 1; Successfully Generated Alphas from nearest neighbor'
!==============================================================================!
!!!!                          Write Alphas to File
!==============================================================================!
!!!open(unit=18,file='alphas.dat')
!!!do i=1,NG
!!!    write(18,*) alpha(i)
!!!enddo
!!!close(18)
!==============================================================================!
!!!!                           Overlap Matrix (S)
!==============================================================================!
!!!do i=1,NG
!!!    do j=i,NG
!!!         aij=alpha(i)*alpha(j)/(alpha(i)+alpha(j))
!!!         r2=sum((x(:,i)-x(:,j))**2)
!!!!different from 2D
!!!         Smat(i,j)=(2*sqrt(alpha(i)*alpha(j))/(alpha(i)+alpha(j)))**(0.5*d)&
!!!             *exp(-aij*r2)
!!!       Smat(j,i)=Smat(i,j)
!!!    enddo
!!!enddo
!==============================================================================!
!!!!                   Check to see if S is positive definite
!!!!If this is removed, you need to allocate llapack arrays before Hamiltonian
!==============================================================================!
!!!lwork=max(1,3*NG-1)
!!!allocate(work(max(1,lwork)))
!!!call dsyev('v','u',NG,Smat,NG,eigenvalues,work,Lwork,info)
!!!write(*,*) 'Info (Initial Overlap Matrix) ==> ', info
!!!RCN=eigenvalues(1)/eigenvalues(NG)
!!!write(*,*) 'RCN =', RCN
!!!open(unit=19,file='overlap_eigenvalues.dat')
!!!do i=1,NG
!!!    write(19,*) eigenvalues(i)
!!!enddo
!!!close(19)
!!!write(*,*) 'Test 2; Overlap Matrix is Positive Definite'
!==============================================================================!
!!!!             Use Gauss Hermit quadrature to evaluate the potential matrix
!!!! z(GH-order) w(GH-order) --- quadrature points and weights
!==============================================================================!
!!!call cgqf(GH_order,6,0d0,0d0,0d0,1d0,z,w)
!!!!different from 2D
!!!w=w/sqrt(pi)
!==============================================================================!
!!!!                   Solve Generalized Eigenvalue Problem
!==============================================================================!
!!!do i=1,NG
!!!  do j=i,NG
!!!     aij=alpha(i)*alpha(j)/(alpha(i)+alpha(j))
!!!     r2=sum((x(:,i)-x(:,j))**2)
!!!     Smat(i,j)=(2*sqrt(alpha(i)*alpha(j))/(alpha(i)+alpha(j)))**(0.5*d)&
!!!         *exp(-aij*r2)
!!!     Smat(j,i)=Smat(i,j)
!==============================================================================!
!!!!                          Kinetic Energy Matrix
!==============================================================================!
!!!     Hmat(i,j)=aij*(d-2*aij*r2)
!==============================================================================!
!!!!                         Potential Energy Matrix
!==============================================================================!
!!!     x_ij(:)=(alpha(i)*x(:,i)+alpha(j)*x(:,j))/(alpha(i)+alpha(j))
!!!     Vij=0d0
!!!     do l1=1,GH_order
!!!        do l2=1,GH_order
!!!           do l3=1,GH_order
!!!              rr(1)=z(l1)
!!!              rr(2)=z(l2)
!!!              rr(3)=z(l3)
!!!              rr=x_ij+rr/sqrt(alpha(i)+alpha(j))
!!!              Vij=Vij+w(l1)*w(l2)*w(l3)*V(rr)
!!!           enddo
!!!        enddo
!!!     enddo
!==============================================================================!
!!!!                      Hamiltonian = Kinetic + Potential
!==============================================================================!
!!!     Hmat(i,j)=(Hmat(i,j)+Vij)*Smat(i,j)
!!!     Hmat(j,i)=Hmat(i,j)
!!!  enddo
!!!enddo
!!!itype=1
!!!call dsygv(itype,'n','u',NG,Hmat,NG,Smat,NG,eigenvalues,work,Lwork,info)
!!!write(*,*) 'info ==> ', info
!!!open(unit=20,file='eigenvalues.dat')
!!!write(20,*) alpha0, eigenvalues(:)
!!!close(20)
!==============================================================================!
!!!!                              Exact Eigenvalues
!==============================================================================!
!!!open(21,File=theory_in)
!!!do i=1,NG
!!!    read(21,*) theory(i)
!!!enddo
!!!close(21)
!!!open(unit=22,file='abs_error.dat')
!!!open(unit=23,file='rel_error.dat')
!!!open(unit=24,file='alpha_abs_error.dat')
!!!open(unit=25,file='alpha_rel_error.dat')
!!!do i=1,NG
!!!        write(22,*) i, abs(theory(i)-eigenvalues(i))
!!!        write(23,*) i, (eigenvalues(i)-theory(i))/theory(i)
!!!        write(24,*) alpha0, abs(theory(i)-eigenvalues(i))
!!!        write(25,*) alpha0, (eigenvalues(i)-theory(i))/theory(i)
!!!enddo
!!!close(22)
!!!close(23)
!!!close(24)
!!!close(25)
!!!open(unit=26,file='alpha_rel_150.dat')
!!!do i=1,150
!!!        write(26,*) alpha0, (eigenvalues(i)-theory(i))/theory(i)
!!!enddo
!!!close(26)
!==============================================================================!
!                               Output File                                    !
!==============================================================================!
open(99,file='simulation.dat')
!write(99,*) 'dimensionality ==> ', d
close(99)
write(*,*) 'Hello Universe!'
end program main
