!=============================================================================80
!                Gaussian Convolution For EigenSpectra Smoothing
!==============================================================================!
!Given vibrational eigenspectra and dipole (amplitude) compute a 
!gaussian-smoothed vibrational spectra.
!If you do not have a dipole surface (assume ampliture=1)==>Density of States
!==============================================================================!
Program Gauss_conv
implicit none
!==============================================================================!
!data_in        ==>Data file containing: eigenvalue dipole
!Nmax           ==>Number of data points in data_in
!K              ==>Number of points to plot the spectrum
!sigma          ==>Gaussian Width
!Spectra        ==>Gaussian-Convolution Spectrum (K evaluations from Fmin-Fmax)
!Fmin/Fmax      ==>Minimum/Maximum Frequency for spectra.dat
!E              ==>Eigenvalue
!d              ==>dipole (amplitude)
!==============================================================================!
character(len=50)::data_in
double precision,parameter::autocm=2.194746313D5
integer::Nmax,K,i
double precision::sigma,F,Fmin,Fmax,Spectrum
double precision,allocatable,dimension(:)::E,d
!==============================================================================!
!                             Read Input Data File
!==============================================================================!
read(*,*) data_in
read(*,*) Nmax
read(*,*) K
read(*,*) sigma                                 !adjust to your data accordingly
read(*,*) Fmin
read(*,*) Fmax
!==============================================================================!
allocate(E(Nmax),d(Nmax))
!==============================================================================!
open(21,File=data_in)
do i=1,Nmax
  read(21,*) E(i),d(i)
enddo
close(21)
!E=E*autocm
E(2:Nmax)=E(2:Nmax)-E(1)
Spectrum=0d0
open(22,File='spectra.dat')
do i=0,K
  F=Fmin+dble(i)*((Fmax-Fmin)/K)
  Spectrum=sum(d(2:Nmax)*exp(-((E(2:Nmax)-F)/sigma)**2))
  write(22,*) F,Spectrum
enddo
close(22)
end Program Gauss_conv
