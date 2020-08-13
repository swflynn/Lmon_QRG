!=============================================================================80
!                       QRG Morse Implementation                               !
!==============================================================================!
!ndmorse implementation for using greedy minimization
!Points only repel no attractions
!==============================================================================!
!    Modified:
!12 August 2020
!    Author:
!Shane Flynn
!==============================================================================!
module morse_mod
implicit none
!==============================================================================!
!Npoints              ==>Number of points to generate
!d                    ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
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
!c_LJ                 ==>Parameter for q-LJ pseudo-potential
!E_cut                ==>Energy Cutoff Contour (kcal/mol input)
!omega                ==>(d) Parameter for Morse Potential
!==============================================================================!
!                            Global Variables                                  !
!==============================================================================!
integer::d,Npoints
double precision,allocatable,dimension(:)::omega
double precision,allocatable,dimension(:,:)::r
double precision::E_cut,integral_P,c_LJ
!==============================================================================!
contains
!==============================================================================!
function P(x)
!==============================================================================!
!Target Distribution Function
!==============================================================================!
!P              ==>evaluate P(x)
!x              ==>(d) ith particles coordinate x^i_1,..,x^i_d
!==============================================================================!
implicit none
double precision::x(d),P
if(V(x)<E_cut) P=(E_cut-V(x))**(d/2.)/integral_P
if(V(x)>=E_cut) P=1d-20
end function P
!==============================================================================!
function V(x)
!==============================================================================!
!Potential Energy (sum of 1D morse potentials)
!==============================================================================!
!x              ==>(d) ith particles coordinate x^i_1,..,x^i_d
!V              ==>evaluate V(x)
!D_morse        ==>Parameter for Morse Potential
!omega          ==>(d) Parameter for Morse Potential
!==============================================================================!
implicit none
double precision::x(d),V
double precision,parameter::D_morse=12.
V=D_morse*sum((exp(-omega(:)*x(:))-1.)**2)
end function V
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
function Pair_rpl(x1,x2,sigma1,sigma2)
!==============================================================================!
!Pairwise repulsive energy between grid points (replaced quas-lennard jones)
!==============================================================================!
!x1             ==>(d) ith atoms coordinates
!x2             ==>(d) jth atoms coordinates
!==============================================================================!
implicit none
double precision::x1(d),x2(d),Pair_rpl,sigma1,sigma2,dist
dist=sqrt(sum((x1(:)-x2(:))**2))
Pair_rpl=(sigma1/dist)**(d+9)+(sigma2/dist)**(d+9)
end function Pair_rpl
!==============================================================================!
end module morse_mod
!==============================================================================!
program main_grid
use morse_mod
!==============================================================================!
implicit none
integer::N_MMC_box,N_1D,N_MMC_grid,MMC_freq,accept,counter,i,j,k,Ntotal
double precision::time1,time2,Delta_E,dummy,moment,mv_cutoff,Pmax,V1,P_trial
double precision::sigma_trial,E_total,V_trial
double precision,allocatable,dimension(:)::delr,index1,U_move,r_trial,s,r_i
double precision,allocatable,dimension(:)::sigma,r0,rmin,rmax
!==============================================================================!
!                              Read Input File                                 !
!==============================================================================!
call cpu_time(time1)
read(*,*) d
allocate(omega(d))
read(*,*) omega
read(*,*) Npoints
read(*,*) E_cut
read(*,*) c_LJ
read(*,*) N_1D
read(*,*) N_MMC_box
read(*,*) N_MMC_grid
read(*,*) MMC_freq
!==============================================================================!
allocate(r(d,Npoints),sigma(Npoints),r_i(d),delr(d),index1(d))
allocate(r_trial(d),s(d),r0(d),rmin(d),rmax(d))
!==============================================================================!
!                     Box Size for normalizing P (MMC)
!==============================================================================!
integral_P=1d0                                     !Initially set to 1 to call P
r_i=0d0
rmin=r_i
rmax=r_i
accept=0
mv_cutoff=1
do i=1,N_MMC_box
  call random_number(s)
  r_trial=r_i+mv_cutoff*(2*s-1)     !trial move (-1,1), random numbers are (0,1)
  if(V(r_trial)<E_cut) then
     r_i=r_trial
     accept=accept+1
     do j=1,d
        if(rmin(j).gt.r_i(j)) rmin(j)=r_i(j)
        if(rmax(j).lt.r_i(j)) rmax(j)=r_i(j)
     enddo
  endif
  if(mod(i,MMC_freq)==0)then
    if(dble(accept)/MMC_freq.lt.0.5)then
      mv_cutoff=mv_cutoff*0.9
    else
      mv_cutoff=mv_cutoff*1.1
    endif
    accept=0
  endif
enddo
!==============================================================================!
!Compute Integral P with square grid         P(x)~Area_Square/N sum_n=1,N P(x_n)
!direct grids get huge, only write if you are dubugging.
!==============================================================================!
!open(20,File='direct_grid.dat')
Moment=0.
Ntotal=(N_1D+1)**d
index1=0
delr(:)=(rmax(:)-rmin(:))/N_1D
do i=1,Ntotal
  do j=1,d
    if(index1(j).eq.N_1D) then
      index1(j)=0
    else
      index1(j)=index1(j)+1
      exit
    endif
  enddo
  r_i(:)=rmin(:)+index1(:)*delr(:)
  Moment=Moment+P(r_i)
enddo
dummy=1./N_1D**d
do j=1,d
  dummy=dummy*(rmax(j)-rmin(j))
enddo
integral_P=dummy*Moment
!close(20)
!==============================================================================!
!           Generate initial distribution to then convert to a QRG
!==============================================================================!
i=1
r0=0d0
Pmax=P(r0)
do while(i.le.Npoints)
  call random_number(r_i)
  r_i(:)=rmin(:)+r_i(:)*(rmax(:)-rmin(:))
  call random_number(dummy)
  if(P(r_i)/Pmax.gt.dummy)then
     r(:,i)=r_i(:)
     sigma(i)=c_LJ*(P(r(:,i))*Npoints)**(-1./d)
     i=i+1
  endif
enddo
E_total=0d0
do i=2,Npoints
   do j=1,i-1
      E_total=E_total+Pair_rpl(r(:,i),r(:,j),sigma(i),sigma(j))
   enddo
enddo
write(*,*) 'E_total=',E_total
open(21,File='grid_ini.dat')
do i=1,Npoints
  write(21,*) r(:,i), P(r(:,i)),sigma(i)
enddo
close(21)
!==============================================================================!
!                             Generate QRG (greedy search)
!==============================================================================!
accept=0
mv_cutoff=1.
do i=1,N_MMC_grid
  k=random_integer(1,Npoints)                               !Select Atom to Move
  call random_number(s)
  r_trial=r(:,k)+mv_cutoff*sigma(k)*(2*s-1)      !random # (0,1), make it (-1,1)
  P_trial=P(r_trial)
  if(P_trial>1d-20) then                                !Only consider if V<Ecut
    sigma_trial=c_LJ*(P_trial*Npoints)**(-1./d)
    Delta_E=0d0
    do j=1,Npoints
      if(j.ne.k) then
        Delta_E=Delta_E+Pair_rpl(r(:,j),r_trial,sigma(j),sigma_trial) - &
        Pair_rpl(r(:,j),r(:,k),sigma(j),sigma(k)) !Energy change from trial move
      endif
    enddo
    if(Delta_E.le.0d0)then
       r(:,k)=r_trial(:)
       sigma(k)=sigma_trial
       accept=accept+1
       E_total=E_total+Delta_E
    endif
  endif
!for MMC want acceptance ~30-50%, adjust trial movement displacement accordingly
  if(mod(i,MMC_freq)==0)then
    if(dble(accept)/MMC_freq.lt.0.4)then
      mv_cutoff=mv_cutoff*0.9
    else
      mv_cutoff=mv_cutoff*1.1
    endif
    accept=0
    write(72,*) i,mv_cutoff,E_total
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
open(99,file='out.txt')
write(99,*) 'd ==> ', d
do i=1,d
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
