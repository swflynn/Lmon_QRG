!=============================================================================80
!                     Metropolis Monte Carlo Code Grid
!==============================================================================!
!       Discussion:
!using optimized grid solve generalized eigenvalue problem
!omegas can all be the same value (set=1) since our basis is symmetric 
!gaussians now due to spacing between gridpoints
!integrals computed iwht quasi Monte Carlo
!==============================================================================!
!This code is intended to use real coordinates vi mass scaling and a Hessian
!This is just a template code, not tested, use at your own risk!!!
!==============================================================================!
!       Modified:
!   1 May 2019
!       Author:
!   Shane Flynn 
!==============================================================================!
module dgb_mod
implicit none
!==============================================================================!
!                            Global Variables 
!==============================================================================!
!d              ==>Gaussian dimensionality
!==============================================================================!
integer::d,NG,Natoms
character(len=2),allocatable::atom_type(:)
double precision,allocatable::mass(:),sqrt_mass(:)
double precision,parameter::Hmass=1d0
double precision,parameter::D_morse=12.
double precision,parameter::c_x=0.2041241
double precision,parameter::c_y=0.18371169
!==============================================================================!
contains
!==============================================================================!
function Atom_Mass(atom)
!==============================================================================!
implicit none
double precision::Atom_Mass
character(len=2)::atom
if(atom=='H'.or.atom=='h')then 
    Atom_mass=Hmass
else 
    write(*,*) 'atom ', atom, ' is not recognized'
    stop 'Check Atom_Mass Function'
endif
end function Atom_Mass
!==============================================================================!
subroutine Toy_Potential(x_i,energies)
!==============================================================================!
!       Discussion:
!Hard-coded Morse Potential Energy 
!==============================================================================!
implicit none
double precision::x_i(d),energies
double precision,parameter::c_x=0.2041241
double precision,parameter::c_y=0.18371169
energies=D_morse*((exp(-c_x*x_i(1))-1)**2+(exp(-c_y*x_i(2))-1)**2)
end subroutine Toy_Potential
!==============================================================================!
subroutine Toy_Force(x,forces)
!==============================================================================!
!       Discussion:
!Returns the Forces associated with Toy_Potential Subroutine
!Forces are hard-coded based on Toy_Potential 
!==============================================================================!
implicit none
integer::i
double precision::x(d),forces(d)
forces(1)=-2.*D_morse*c_x*exp(-2.*c_x*x(1))*(exp(c_x*x(1))-1)
forces(2)=-2.*D_morse*c_y*exp(-2.*c_y*x(2))*(exp(c_y*x(2))-1)
!write(*,*) 'Forces from Toy_Force Subroutine ==> ', forces
end subroutine Toy_Force
!==============================================================================!
subroutine Toy_Hessian(x,Hess_Mat)
!==============================================================================!
!       Discussion:
!Numerically computed Hessian using forces from "Toy_Force Subroutine"
!       Variables:
!s          ==> Perturbation Parameter 
!Hess_Mat   ==> (d,d); Symmetrized Mass-Scaled Hessian
!x          ==> (d); XYZ Configuration 
!==============================================================================!
implicit none 
integer::i,j
double precision::Hess_Mat(d,d),x(d),r(d),force(d)
double precision::force0(d)
double precision,parameter::s=1d-6
r=x
call Toy_Force(r,force0)
do i=1,d
    r(i)=x(i)+s
    call Toy_Force(r, force)
    r(i)=x(i)
    do j=1,d
        Hess_Mat(i,j)=(force0(j)-force(j))/s
    enddo
enddo
!==============================================================================!
!                   Symmetrize and Mass-Scale the Hessian
!==============================================================================!
do i=1,d
    do j=1,i
        if(i.ne.j) Hess_Mat(i,j)=(Hess_Mat(i,j)+Hess_Mat(j,i))/2
        Hess_Mat(i,j)=Hess_Mat(i,j)/(sqrt_mass(i)*sqrt_mass(j))
        if(i.ne.j) Hess_Mat(j,i)=Hess_Mat(i,j)
    enddo
enddo
!write(*,*) 'Hessian from Toy_Hessian Subroutine ==> ', Hess_Mat
end subroutine Toy_Hessian
!==============================================================================!
subroutine Freq_Hess(d,Hess,omega,U)
!==============================================================================!
!       Discussion:
!Compute Eigenvalues and Eigenvectors of Hessian
!Uses the LLAPACK real symmetric eigen-solver (dsygev)
!       Variables:
!Hess   ==> (d,d); Hessian Matrix
!omega  ==> (d); Hessian Eigenvalues
!U      ==> (d,d); Hessian Eigenvectors
!       LLAPACK (dsyev):
!v      ==> Compute both Eigenvalues and Eigenvectors
!u      ==> Use Upper-Triangle of matrix
!==============================================================================!
implicit none
integer::i,info,lwork,d
double precision::Hess(d,d),omega(d),U(d,d)
double precision,allocatable::work(:) 
lwork=max(1,3*d-1)          !suggested by LAPACK Developers
allocate(work(max(1,lwork)))    !suggested by LAPACK Developers 
U=Hess  
call dsyev('v','u',d,U,d,omega,work,lwork,info) 
write(*,*) 'Frequencies from the Hessian:'
do i=d,1,-1
    omega(i)=sign(sqrt(abs(omega(i))),omega(i))
    write(*,*) omega(i), 'normalized = 1?', sum(U(:,i)**2)
enddo
end subroutine Freq_Hess
!==============================================================================!
end module dgb_mod
!==============================================================================!
!==============================================================================!
program main
use dgb_mod
!==============================================================================!
!       Discussion:
!==============================================================================!
!d              ==> i-th gaussian dsionality (x^i=x^i_1,x^i_2,..,x^i_d)
!NG             ==> Number of basis functions
!x              ==>(d) all atom coordinates
!U              ==>(NG,NG) All i,j pair-wise energies (LJ-12-6)
!x0             ==>(d) store previous coordinate before trial move
!accept         ==> number of accepted trial moves, for acceptance~50%  
!counter        ==> total number of moves, for acceptance~50%
!x2 ith gaussians coordinates for matrix elements
!z integration sequence form matrix elements
!x_ij= i-jth gaussian center
!U_Hess(d,d) eigenvectors hessian, for now set equal to 1
!z=sequence for integration
!t_i,t_f        ==> cpu time to ~ simulation time
!==============================================================================!
implicit none
character(len=50)::coord_in,grid_in,alpha_in
integer::Nsobol,data_freq,n,i,j
integer*8::ii,skip
double precision::E0,s_sum,time1,time2
double precision,allocatable,dimension(:)::q0,omega,alpha,eigenvalues
double precision,allocatable,dimension(:,:)::x,Hess,U,Smat,z
!==============================================================================!
!                       LLAPACK dsygv variables                                !
!==============================================================================!
integer::itype,info,lwork
double precision,allocatable,dimension(:)::work
!==============================================================================!
!                           Read Input Data File                               !
!==============================================================================!
call cpu_time(time1)
read(*,*) d
read(*,*) NG
read(*,*) coord_in
read(*,*) grid_in
read(*,*) alpha_in
read(*,*) Nsobol
read(*,*) data_freq
open(15,File=coord_in)
read(15,*) Natoms
read(15,*)
!==============================================================================!
!                               Allocations
!==============================================================================!
allocate(atom_type(Natoms),mass(Natoms),sqrt_mass(d),q0(d))
allocate(x(d,NG),Hess(d,d),U(d,d),omega(d),alpha(NG),eigenvalues(NG),Smat(NG,NG))
allocate(z(d,Nsobol))
write(*,*) 'Test 0; Successfully Read Input File'
!==============================================================================!
!                         Input Configuration Energy
!==============================================================================!
do i=1,Natoms
    read(15,*) atom_type(i),q0(2*i-1:2*i)   
    mass(i)=Atom_mass(atom_type(i))
    sqrt_mass(2*i-1:2*i)=sqrt(mass(i))
enddo
close(15)
write(*,*) 'd ==> ', d
call toy_potential(q0,E0)       
write(*,*) 'E0 ==> ', E0
!==============================================================================!
! 			Compute Hessian and Frequencies
!==============================================================================!
call Toy_Hessian(q0,Hess)
call Freq_Hess(d,Hess,omega,U)
write(*,*) 'Test 2; Successfully Computed Hessian'
!==============================================================================!
!                           Input GridPoints x(d,NG)
!==============================================================================!
open(16,File=grid_in)
do n=1,NG
    read(16,*) x(:,n)
enddo 
close(16)
!==============================================================================!
!                               Read Alphas
!==============================================================================!
open(17,File=alpha_in)
do n=1,NG
    read(17,*) alpha(n)
enddo 
close(17)
!==============================================================================!
!                           Overlap Matrix (S)
!==============================================================================!
do i=1,NG
    do j=i,NG
        s_sum=sum(omega(:)*(x(:,i)-x(:,j))**2)
        Smat(i,j)=sqrt(alpha(i)+alpha(j))**(-d)&
                 *exp(-0.5*alpha(i)*alpha(j)/(alpha(i)+alpha(j))*s_sum)
        Smat(j,i)=Smat(i,j)
    enddo
enddo
!==============================================================================!
!                   Check to see if S is positive definite
!If this is removed, you need to allocate llapack arrays before Hamiltonian 
!==============================================================================!
lwork=max(1,3*NG-1)
allocate(work(max(1,lwork)))
call dsyev('v','u',NG,Smat,NG,eigenvalues,work,Lwork,info)
write(*,*) 'Info (Initial Overlap Matrix) ==> ', info
open(unit=19,file='overlap_eigenvalues.dat')
do i=1,NG
    write(19,*) eigenvalues(i)
enddo
close(19)
write(*,*) 'Test 4; Overlap Matrix is Positive Definite'
!==============================================================================!
!                   Generate Sequence For Evaluating Potential
!==============================================================================!
skip=Nsobol
do ii=1,Nsobol
    call sobol_stdnormal(d,skip,z(:,ii))
enddo
write(*,*) 'Test 5; Successfully Generated Integration Sequence'
!==============================================================================! 
!                             Evaluate Potential
!==============================================================================!
Vmat=0d0
counter=0
open(unit=21,file='eigenvalues.dat')
do while(counter.lt.Nsobol/data_freq)
    do i=1,NG
        do j=i,NG
            x_ij(:)=(alpha(i)*x(:,i)+alpha(j)*x(:,j))/(alpha(i)+alpha(j))
            do l=1+(counter*data_freq),(counter+1)*data_freq
                x2(:)=x_ij(:)+z(:,l)/sqrt(omega(:)*(alpha(i)+alpha(j)))
!==============================================================================!
!                           Transformation to Cartesian
!==============================================================================!
                Vmat(i,j)=Vmat(i,j)+V_x(matmul(U_Hess,x2)) 
            enddo
        enddo
    enddo
!==============================================================================!
!                   Solve Generalized Eigenvalue Problem
!==============================================================================!
    write(*,*) 'Generalized Eigenvalue Iteration', (counter+1)*data_freq
    Hmat=Vmat
    do m=1,NG
        do n=m,NG
            s_sum=sum(omega(:)*(x(:,m)-x(:,n))**2)
            Smat(m,n)=sqrt(alpha(m)+alpha(n))**(-d)&           
                     *exp(-0.5*alpha(m)*alpha(n)/(alpha(m)+alpha(n))*s_sum)
            Smat(n,m)=Smat(m,n)
            kinetic=Smat(m,n)*0.5*alpha(m)*alpha(n)/(alpha(m)+alpha(n))&
            *sum(omega(:)-(alpha(m)*alpha(n)*(omega(:)**2*(x(:,m)-x(:,n))**2)&
            /(alpha(m)+alpha(n))))
            Hmat(m,n)=Hmat(m,n)*(Smat(m,n)/((counter+1)*data_freq))+kinetic
            Hmat(n,m)=Hmat(m,n)
        enddo
    enddo
        itype=1
        eigenvalues=0d0
        call dsygv(itype,'n','u',NG,Hmat,NG,Smat,NG,eigenvalues,work,Lwork,info)
        write(21,*) eigenvalues(:)
        write(*,*) 'info ==> ', info
        counter=counter + 1
enddo
close(21)
!==============================================================================!
!                               output file                                    !
!==============================================================================!
call cpu_time(time2)
open(90,file='simulation.dat')
write(90,*) 'particle dimensionality ==> ', d
write(90,*) 'NG ==> ', NG
write(90,*) 'Nsobol==>', Nsobol
write(90,*) 'Total Time ==> ', time2-time1
close(90)
write(*,*) 'Hello Universe!'
end program main
