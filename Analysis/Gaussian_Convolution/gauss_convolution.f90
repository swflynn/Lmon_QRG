!=============================================================================80
!                Gaussian Convolution For EigenSpectra Smoothing
!==============================================================================!
!Given vibrational eigenspectra and amplitudes compute a gaussian-smoothed
!vibrational spectra.
!==============================================================================!
Program Gauss_conv
implicit none
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
read(*,*) sigma
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
