!=================20==================40==================60==================80
!                     Quasi-Regular Grid Generation Main
!==============================================================================!
!Generate Quasi-Regular Grid for any d-dimensional system of interest.
!==============================================================================!
!       Modified:
!   29 March 2020
!       Author:
!   Shane Flynn
!==============================================================================!
program main
use qrg_grid_mod
!==============================================================================!
implicit none
integer,parameter::d1=9
double precision,parameter::bohr=0.52917721092
double precision,parameter::autocm=2.194746313D5
double precision,parameter::autokcalmol=627.5096
!==============================================================================!
character(len=20)::potential
character(len=2),allocatable,dimension(:)::atom_type
integer::d,d2,Natoms,Npoints,N_MMC_box,N_1D,N_MMC_grid,MMC_freq,Lwork0
double precision::E0,E_cut,c_LJ,integral_P,V,time1,time2
double precision,allocatable,dimension(:)::mass,sqrt_mass,x0,forces,omega,rmin
double precision,allocatable,dimension(:)::rmax
double precision,allocatable,dimension(:,:)::Hess_Mat,U,r,Uij
!==============================================================================!
!potential            ==>Potential name
!Npoints              ==>Number of points to generate
!d                    ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!r(d2,Npoints)        ==>All grid points coordinates (x^i=x_1,x_2,...x_d)
!Uij(Npoints,Npoints) ==>Pairwise energies for all i-j points
!V_i                  ==>Potential Energy evaluation V(x_i)
!E_cut                ==>Distribution cutoff contour
!rmin(d)              ==>Minimum of normalization box size
!rmax(d)              ==>Maximum of normalization box size
!N_1D                 ==>Number of points in 1 dimension for computing integral_P
!N_MMC_box            ==>Number of MMC Iterations to determine box-size
!c_LJ                 ==>Parameter for q-LJ pseudo-potential
!N_MMC_grid           ==>Number of MMC Iterations to optimize QRG
!MMC_freq             ==>Frequency to update QRG MMC grid mv_cutoff
!integral_P           ==>Normalization constant for the distribtion P(x)
!==============================================================================!
call cpu_time(time1)
read(*,*) Npoints
read(*,*) d2
read(*,*) N_MMC_box
read(*,*) E_cut
read(*,*) N_1D
read(*,*) c_LJ
read(*,*) N_MMC_grid
read(*,*) MMC_freq
read(*,*) potential
!==============================================================================!
call get_xyz_parameters(d,Natoms)
!==============================================================================!
allocate(atom_type(Natoms),mass(Natoms),sqrt_mass(d),x0(d),forces(d),rmin(d2))
allocate(rmax(d2),Hess_Mat(d1,d1),omega(d1),U(d1,d1),r(d2,Npoints))
allocate(Uij(Npoints,Npoints))
Lwork0=max(1,3*d1-1)
!==============================================================================!
call read_input_geometry(d,Natoms,atom_type,mass,sqrt_mass,x0)
call convert_to_atomic_units(potential,d,x0,E0,forces,E_cut,bohr,autocm,&
  autokcalmol)
call Get_Hessian(potential,d,d1,x0,Hess_Mat)
call Mass_Scale_Hessian(d,d1,Hess_Mat,sqrt_mass)
call Frequencies_Scaled_Hess(d1,Hess_mat,omega,U,Lwork0)
!==============================================================================!
call box_size_P(potential,d,d1,d2,x0,U,sqrt_mass,E_cut,rmin,rmax,N_MMC_box,&
  integral_P)
call compute_integral_P(potential,d,d1,d2,x0,V,U,sqrt_mass,E_cut,integral_P,&
  rmin,rmax,N_1D)
call initial_distribution(potential,d,d1,d2,x0,r,V,U,sqrt_mass,Npoints,E_cut,&
  rmin,rmax)
call initial_pair_energy(potential,d,d1,d2,x0,V,U,sqrt_mass,E_cut,integral_P,&
  c_LJ,Npoints,r,Uij)
call Quasi_Regular(potential,d,d1,d2,x0,V,U,sqrt_mass,E_cut,integral_P,r,Uij,&
  c_LJ,N_MMC_grid,MMC_freq,Npoints)
call write_out(potential,Npoints,d,d1,d2,E_cut,rmin,rmax,N_1D,N_MMC_box,c_LJ,&
  N_MMC_grid,MMC_freq,integral_P)
call cpu_time(time2)
write(*,*) 'CPU time ==> ', time1-time2
end program main
