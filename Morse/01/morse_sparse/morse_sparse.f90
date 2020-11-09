!=============================================================================80
!                  ND-Morse EigenSpectra (Sparse Quadrature)
!==============================================================================!
!Compute the EigenSpectra for an nD-Morse Potential
!Quasi-Regular Grid as Input; generate (Gaussian Widths) from nearest neghbor
!This code needs the gen_hermite_rule.f90 code (Gauss-Hermite-Quadriture)
!This code uses LLAPACK to solve the generalized eigenvalue problem
!==============================================================================!
!       Modified:
!   2 November 2020
!       Author:
!   Shane Flynn
!==============================================================================!
!start writing morse test case, use main_spare.f90 as template
!==============================================================================!
module QRGB_mod
implicit none
!==============================================================================!
!                            Global Variables
!==============================================================================!
!d              ==>i-th gaussian dimensionality (x^i=x^i_1,x^i_2,..,x^i_d)
!NG             ==>Number of Gaussian to generate
!omega          ==>(d) Parameter for Morse Potential
!==============================================================================!
integer::d,NG
double precision,allocatable,dimension(:)::omega,alpha
double precision,allocatable,dimension(:,:)::x
!==============================================================================!
contains
!==============================================================================!
function V(x)
!==============================================================================!
!Sum of 1D-Morse Potential
!==============================================================================!
!x              ==>(d) ith atoms coordinates
!V              ==>evaluate V(x_i)
!D_morse        ==>Parameter for Morse Potential
!==============================================================================!
implicit none
double precision::x(d),V
double precision,parameter::D_morse=12.
V=D_morse*sum((exp(-omega(:)*x(:))-1.)**2 )
end function V
!==============================================================================!
subroutine matrix_elements_quad(i,j,S_ij,H_ij,quad_size,z,w)
!==============================================================================!
integer::i,j,k,quad_size
double precision,parameter::ss=1d-5                               !need to check
double precision::S_ij,H_ij,x2,a_ij,E,E1,E2,V_ij,Lap
double precision::x_ij(d),xx(d),z(d,quad_size),w(quad_size),x1(d)
!==============================================================================!
a_ij=alpha(i)*alpha(j)/(alpha(i)+alpha(j))
x2=sum((x(:,i)-x(:,j))**2)
!==============================================================================!
!               Compute i-jth element of the overlap matrix
!==============================================================================!
S_ij=(2.*sqrt(alpha(i)*alpha(j))/(alpha(i)+alpha(j)))**(0.5*d)*exp(-a_ij*x2)
!==============================================================================!
!               Compute i-jth element of the kinetic matrix
!==============================================================================!
H_ij=a_ij*(d-2.*a_ij*x2)
!==============================================================================!
!    Compute Potential Energy by expanding the potential to second order
!==============================================================================!
x_ij(:)=(alpha(i)*x(:,i)+alpha(j)*x(:,j))/(alpha(i)+alpha(j))
!==============================================================================!
if(quad_size==0)then                              !quadratic interpolation of V
  x1=x_ij
  E=V(x1)
  Lap=0d0
  do k=1,d
    x1(k)=x_ij(k)+ss
    E1=V(x1)
    x1(k)=x_ij(k)-ss
    E2=V(x1)
    Lap=Lap-2*E+E1+E2
    x1(k)=x_ij(k)
  enddo
  V_ij=E+Lap/(ss**2*4*(alpha(i)+alpha(j)))
else                                                             !Use quadrature
  V_ij=0d0
  do k=1,quad_size
    xx=x_ij+z(:,k)/sqrt(alpha(i)+alpha(j))
    E=V(xx)
    V_ij=V_ij+E*w(k)
  enddo
  H_ij=(H_ij+V_ij)*S_ij
end if
end subroutine matrix_elements_quad
!==============================================================================!
end module QRGB_mod
!==============================================================================!
program QRG_Spec
use QRGB_mod
!==============================================================================!
!unif_grid      ==>If True set alpha to be a constant, else use nearest neighbor
!grid_in        ==>Filename Containing Gridpoints
!theory_in      ==>Filename Containing Analytic Eigenvalues
!NG             ==>Number of Gaussian Basis Functions (gridpoints)
!GH_order       ==>Number of Points for evaluating the potential (Gauss-Hermite)
!alpha0         ==>Flat Scaling Parameter for Gaussian Widths
!alpha          ==>(d) Gaussian Widths
!RCN            ==>Recriprical Convergence Number, stability of Overlap Matrix
!x              ==>(d) ith atoms coordinates
!x_ij           ==>i-jth gaussian center (product of Gaussians is a Gaussian)
!Smat           ==>(NG,NG) Overlap Matrix
!Hmat           ==>(NG,NG) Hamiltonian Matrix
!==============================================================================!
implicit none
character(len=50)::grid_in,theory_in
integer::GH_order,sg_level,quad_size,i,j,k,tsize
external gqn,gqn_order                                        !sparse_quadrature
integer,allocatable,dimension(:)::l
double precision::x2,alpha0,time1,time2
double precision,parameter::pi=4.*atan(1d0)
double precision,allocatable,dimension(:)::eigenvalues,z1,w,w1,theory
double precision,allocatable,dimension(:,:)::Smat,Hmat,z
!==============================================================================!
!                            LLAPACK dsygv variables
!==============================================================================!
integer::itype,info,lwork
double precision,allocatable,dimension(:)::work
!==============================================================================!
!                             Read Input Data File
!==============================================================================!
call cpu_time(time1)
read(*,*) d
allocate(omega(d))
read(*,*) omega
read(*,*) NG
read(*,*) GH_order
read(*,*) grid_in
read(*,*) theory_in
read(*,*) tsize     !how many theory eigenvalues to read in
read(*,*) alpha0
!==============================================================================!
!                               Allocations
!==============================================================================!
allocate(x(d,NG),alpha(NG),eigenvalues(NG),Smat(NG,NG),Hmat(NG,NG),theory(tsize))
!==============================================================================!
!                           Read GridPoints x(d,NG)
!==============================================================================!
open(17,File=grid_in)
do i=1,NG
  read(17,*) x(:,i)
enddo
close(17)
!==============================================================================!
!         Symmetric gaussians. Use nearest neighbor to determine alpha
!==============================================================================!
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
open(unit=18,file='alphas.dat')
do i=1,NG
  write(18,*) alpha(i)
enddo
close(18)
!==============================================================================!
!                   Generate quadrature points and weights
!==============================================================================!
if(GH_order>0)then                     !Direct product Gauss-Hermite quadrature
  quad_size=GH_order**d
  allocate(z(d,quad_size),w(quad_size),z1(GH_order),w1(GH_order),l(d))
  call cgqf(GH_order,6,0d0,0d0,0d0,1d0,z1,w1)                               !GHQ
  w1=w1/sqrt(pi)
  l(:)=1
  do i=1,quad_size
    do k=1,d
      z(k,i)=z1(l(k))
    enddo
    w(i)=product(w1(l(:)))
    !w(i)=w1(l(1))
    !do k=2,d2
      !w(i)=w(i)*w1(l(k))
    !enddo
    do k=1,d
      l(k)=mod(l(k),GH_order)+1
      if(l(k).ne.1) exit
    enddo
  enddo
elseif(GH_order<0)then                                              !Sparse grid
!==============================================================================!
!For given growth order GQN_ORDER, dimension d, and level sg_l
!get the number of sample points and weights.
!==============================================================================!
  sg_level=-GH_order
!==============================================================================!
!For given growth order GQN_ORDER, dimension d, and level sg_l
!get the number of sample points and weights.
!==============================================================================!
  call nwspgr_size(gqn_order,d,sg_level,quad_size)
  allocate(z(d,quad_size),w(quad_size))
!==============================================================================!
!                     Compute sample points and weights
!==============================================================================!
  i=quad_size
  call nwspgr(gqn,gqn_order,d,sg_level,i,quad_size,z,w)
!==============================================================================!
!(2 pi)^(-D/2) int( -oo, +oo )^D F(X) exp(-(X1^2+X2^2+..+XD^2)/2) dX1 dX2...dXD
!=\sum_{n=1^quand_size} F(z_n) w_n
!Therefore, we scale z and the scaling for w is canceled:
!==============================================================================!
  z=z/sqrt(2.)
!==============================================================================!
else
  quad_size=0                                       !Use quadratic interpolation
endif
!==============================================================================!
do i=1,NG
  do j=i,NG
    if(i==j) then
      call  matrix_elements_quad(i,j,Smat(i,j),Hmat(i,j),quad_size,z,w)
    else
      call  matrix_elements_quad(i,j,Smat(i,j),Hmat(i,j),quad_size,z,w)
    endif
    Hmat(j,i)=Hmat(i,j)
    Smat(j,i)=Smat(i,j)
  enddo
enddo
!==============================================================================!
!               Eigenvalues of the Hamiltonian matrix
!==============================================================================!
lwork=max(1,3*NG-1)
allocate(work(max(1,lwork)))
itype=1
call dsygv(itype,'n','u',NG,Hmat,NG,Smat,NG,eigenvalues,work,Lwork,info)
write(*,*) 'info ==> ', info
open(unit=20,file='eigenvalues.dat')
write(20,*) alpha0, eigenvalues(:)
close(20)
!==============================================================================!
!                              Exact Eigenvalues
!==============================================================================!
open(21,File=theory_in)
do i=1,tsize
  read(21,*) theory(i)
enddo
close(21)
open(unit=22,file='abs_error.dat')
open(unit=23,file='rel_error.dat')
open(unit=24,file='alpha_abs_error.dat')
open(unit=25,file='alpha_rel_error.dat')
do i=1,tsize
  write(22,*) i, abs(theory(i)-eigenvalues(i))
  write(23,*) i, (eigenvalues(i)-theory(i))/theory(i)
  write(24,*) alpha0, abs(theory(i)-eigenvalues(i))
  write(25,*) alpha0, (eigenvalues(i)-theory(i))/theory(i)
enddo
close(22)
close(23)
close(24)
close(25)
open(unit=26,file='alpha_rel_150.dat')
do i=1,150
  write(26,*) alpha0, (eigenvalues(i)-theory(i))/theory(i)
enddo
close(26)
!==============================================================================!
!                                Output file                                   !
!==============================================================================!
call cpu_time(time2)
open(99,file='simulation.dat')
write(99,*) 'dimensionality ==> ', d
write(99,*) 'omega ==> ', omega
write(99,*) 'NG ==> ', NG
write(99,*) 'alpha0==>', alpha0
write(99,*) 'GH Order==>', GH_order
write(99,*) 'Quad size==>', quad_size
write(99,*) 'theory size ==>', tsize
write(99,*) 'time ==>', time2-time1
close(99)
write(*,*) 'Hello Universe!'
end program QRG_Spec
