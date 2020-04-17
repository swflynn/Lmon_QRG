!=================20==================40==================60==================80
!                              Potentials Module
!==============================================================================!
!Any user defined potential is valid assuming a single point x_i = (x_1,..x_d)
!results in a single energy evaluation [V_i(x_i) = #]
!After defining a potential add the appropriate call statement to the
!"potentials" subroutine in grids_mod.f90
!All equations in grid_generation and vibrational_spectra are derived
!and implemented as d-dimensional.
!==============================================================================!
!       Modified:
!   24 March 2020
!       Author:
!   Shane Flynn
!==============================================================================!
module potential_mod
implicit none
!==============================================================================!
contains
!==============================================================================!
subroutine morse(d,x_i,V_i)
!==============================================================================!
!Potential Energy: Sum of 1D Morse potentials
!==============================================================================!
!d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x_i(d)           ==>i-th grid points coordinates (x^i=x_1,x_2,...x_d)
!D_morse          ==>Parameter for Morse Potential
!omega(d)         ==>Parameter for Morse Potential
!V_i              ==>Potential Energy evaluation V(x_i)
!==============================================================================!
implicit none
integer::d
double precision::x_i(d),V_i,omega(d)
double precision,parameter::omega2(2)=(/0.2041241,0.18371169/)
double precision,parameter::omega3(3)=(/0.2041241,0.18371169,0.16329928/)
double precision,parameter::D_morse=12.
if(d==2)then
  omega=omega2
elseif(d==3)then
  omega=omega3
else
  stop 'Subroutine morse omega not defined, check potentials_mod.f90'
endif
V_i=D_morse*sum((exp(-omega(:)*x_i(:))-1.)**2 )
end subroutine morse
!==============================================================================!
subroutine henon(d,x_i,V_i)
!==============================================================================!
!Potential Energy: 2d Henon-Heiles
!==============================================================================!
!d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x_i(d)           ==>i-th grid points coordinates (x^i=x_1,x_2,...x_d)
!V_i              ==>Potential Energy evaluation V(x_i)
!lambda           ==>Parameter for Henon-Heiles Potential
!==============================================================================!
implicit none
integer::d
double precision::x_i(d),V_i
double precision,parameter::lambda=sqrt(0.0125)
V_i=0.5*(x_i(1)**2+x_i(2)**2)+lambda*(x_i(1)**2*x_i(2)-x_i(2)**3/3.)
end subroutine henon
!==============================================================================!
subroutine toy_3d(d,x,V_i)
!==============================================================================!
!Potential Energy (hard-coded 3D test potential)
!V:= c1 x^4 + c2 x^2 + c3 y^4 + c4 y^2 + c5 z^4 + c6 z^2 + c7 xy + c8 yz + c9 xz
!==============================================================================!
!d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x_i(d)           ==>i-th grid points coordinates (x^i=x_1,x_2,...x_d)
!V_i              ==>Potential Energy evaluation V(x_i)
!c_i              ==>Parameters for Potential Scaling
!==============================================================================!
implicit none
integer::d
double precision::x(d),V_i
double precision,parameter::c1=0.1
double precision,parameter::c2=0.5
double precision,parameter::c3=0.2
double precision,parameter::c4=0.6
double precision,parameter::c5=0.3
double precision,parameter::c6=0.7
double precision,parameter::c7=1.
double precision,parameter::c8=2.
double precision,parameter::c9=3.
V_i=c1*x(1)**4+c2*x(1)**2+c3*x(2)**4+c4*x(2)**2+c5*x(3)**4+c6*x(3)**2+c7*x(1)*&
x(2)+c8*x(2)*x(3)+c9*x(1)*x(3)
end subroutine toy_3d
!==============================================================================!
end module potential_mod
