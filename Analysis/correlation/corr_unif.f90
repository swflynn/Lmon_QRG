!==============================================================================!
!                            Correlation Function 
!==============================================================================!
!       Discussion:
!Pair Correlation Function for a uniform distribution
!Checking to see if these distributions have a clear "solvation shell" type 
!distribution...Conclusion, they DO NOT. 
!==============================================================================!
!       Modified:
!   10 April 2019
!       Author:
!   Shane Flynn 
!==============================================================================!
program paircorr
!==============================================================================!
implicit none
integer::Dimen,m,xmin,hissize
double precision::delx,dist
integer::Nsobol,i,j
integer*8::ii,skip
double precision,allocatable::z(:,:)
double precision,allocatable::hist(:)
!==============================================================================!
!                          Read Input Data File
!==============================================================================!
read(*,*) xmin
read(*,*) delx
read(*,*) Dimen
read(*,*) Nsobol
skip=Nsobol
hissize=100000
allocate(z(Dimen,Nsobol),hist(hissize))
!==============================================================================!
open(unit=66,file='points.dat')
!==============================================================================!
do i=1,Nsobol
    call sobol_stdnormal(Dimen,skip,z(:,i))
    write(66,*) z(:,i)
enddo
close(66)
hist=0
do i=1,Nsobol-1
    do j=i+1,Nsobol
        dist=sqrt(sum((z(:,i)-z(:,j))**2))
        m=(dist-xmin)/delx
        hist(m) = hist(m) + 1
    enddo
enddo
hist=hist / (sum(hist)*delx)

open(unit=67,file='hist.dat')
do i=1,hissize
    write(67,*) i, hist(i)
enddo

open(unit=68,file='sim.dat')
write(68,*) 'Nsobol ', Nsobol
write(68,*) 'Dimensionality', dimen
close(68)
end program paircorr
