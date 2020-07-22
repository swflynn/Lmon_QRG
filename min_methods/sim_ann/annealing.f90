!=============================================================================80
!                       QRG Morse Implementation                  !
!==============================================================================!
!simple ndmorse implementation for easy testing
!==============================================================================!
!    Discussion:
!    Modified:
!30 June 2020
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
!x0(d)                ==>initial cluster configuration
!c_LJ                 ==>Parameter for q-LJ pseudo-potential
!E_cut                ==>Energy Cutoff Contour (kcal/mol input)
!==============================================================================!
!                            Global Variables                                  !
!==============================================================================!
integer::d,Npoints
double precision,allocatable,dimension(:)::omega,rmin,rmax
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
function Pair_LJ_NRG(x1,x2)
!==============================================================================!
!quasi-Lennard Jones pairwise energy between grid points
!==============================================================================!
!x1             ==>(d) ith atoms coordinates
!x2             ==>(d) jth atoms coordinates
!a/b            ==>evaluate q-LJ
!sigma          ==>Gaussian widths (c*sigma(P))
!Pair_LJ_NRG    ==>Energy of the i-j q-LJ potential
!==============================================================================!
implicit none
double precision::x1(d),x2(d),a,b,Pair_LJ_NRG,sigma1,sigma2
a=sum((x1(:)-x2(:))**2)
sigma1=c_LJ*(P(x1)*Npoints)**(-1./d)
sigma2=c_LJ*(P(x2)*Npoints)**(-1./d)
b=(sigma2**2/a)
a=(sigma1**2/a)
Pair_LJ_NRG=a**(d+9)-a**(d+3)+b**(d+9)-b**(d+3)
end function Pair_LJ_NRG
!==============================================================================!
subroutine T_schedule(cooling,T,cc)
!==============================================================================!
!Select a cooling schedule for the simulated annealing
!==============================================================================!
implicit none
double precision::T,cc
character(len=50)::cooling
!==============================================================================!
if(cooling=='linear'.or.cooling=='LINEAR') then
  call linear_cooling(T,cc)
elseif(cooling=='geometric'.or.cooling=='GEOMETRIC') then
  call geometric_cooling(T,cc)
elseif(cooling=='log'.or.cooling=='LOG') then
  call log_cooling(T,cc)
else
  stop 'Cannot Identify Cooling Schedule, Check "T_schedule" Subroutine'
endif
end subroutine T_schedule
!==============================================================================!
subroutine linear_cooling(T,cc)
!==============================================================================!
!assumes cc>0
!==============================================================================!
implicit none
double precision::T,cc
!==============================================================================!
T=T-cc
end subroutine linear_cooling
!==============================================================================!
subroutine geometric_cooling(T,cc)
!==============================================================================!
!assumes 0<cc<1
!==============================================================================!
implicit none
double precision::T,cc
!==============================================================================!
T=T*cc
end subroutine geometric_cooling
!==============================================================================!
subroutine log_cooling(T,cc)
!==============================================================================!
!log cooling to vary T
!==============================================================================!
implicit none
double precision::T,cc
!==============================================================================!
T=T/(1+cc*T)
end subroutine log_cooling
!==============================================================================!
end module morse_mod
!==============================================================================!
program main_grid
use morse_mod
!==============================================================================!
implicit none
character(len=50)::cooling
integer::N_MMC_box,N_1D,i,j,k,Ntotal,T_iter,counter,accept,reject,out_cut
double precision::time1,time2,Delta_E,deltae1,dummy,moment,mv_cutoff,t1,T,T_i
double precision::T_f,cc
double precision,allocatable,dimension(:)::delr,index1,U_move,rr,r_trial,s,r_i
double precision,allocatable,dimension(:,:)::Uij
!==============================================================================!
!                              Read Input File                                 !
!==============================================================================!
call cpu_time(time1)
read(*,*) d
allocate(omega(d))
read(*,*) omega
read(*,*) Npoints
read(*,*) N_1D
read(*,*) N_MMC_box
read(*,*) E_cut
read(*,*) c_LJ
read(*,*) T_i
read(*,*) T_f
read(*,*) T_iter
read(*,*) cc
read(*,*) cooling
!==============================================================================!
allocate(r(d,Npoints),Uij(Npoints,Npoints),r_i(d),delr(d),index1(d))
allocate(U_move(Npoints),rr(d),r_trial(d),s(d))
!==============================================================================!
!                     Box Size for normalizing P (MMC)
!==============================================================================!
integral_P=1d0                                     !Initially set to 1 to call P
r_i=0d0
rmin=r_i
rmax=r_i
mv_cutoff=0.1
do i=1,N_MMC_box
  call random_number(s)
  r_trial=r_i+mv_cutoff*(2*s-1)     !trial move (-1,1), random numbers are (0,1)
  call random_number(dummy)                             !MMC acceptance criteria
  if(P(r_trial)/P(r_i).ge.dummy) then
    r_i=r_trial
    do j=1,d
      if(rmin(j).gt.r_i(j)) rmin(j)=r_i(j)
      if(rmax(j).lt.r_i(j)) rmax(j)=r_i(j)
    enddo
  endif
enddo
!write(*,*) 'test 1'
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
  dummy=P(r_i)
  Moment=Moment+dummy
!  if(V.lt.E_cut) write(20,*) r_i
enddo
dummy=1./N_1D**d
do j=1,d
  dummy=dummy*(rmax(j)-rmin(j))
enddo
integral_P=dummy*Moment
!close(20)

!write(*,*) 'test 2'
!write(*,*) 'integral p', integral_P
!==============================================================================!
!           Generate initial distribution to then convert to a QRG
!==============================================================================!
i=1
do while(i.le.Npoints)
  call random_number(r_i)
  r_i(:)=rmin(:)+r_i(:)*(rmax(:)-rmin(:))
  if(V(r_i).lt.E_cut)then
    r(:,i)=r_i(:)
    i=i+1
  endif
enddo
open(21,File='grid_ini.dat')
do i=1,Npoints
  write(21,*) r(:,i)
enddo
close(21)
!==============================================================================!
!   Compute pairwise energies for all the initial Grid Points U[x_{ij}]
!==============================================================================!
do i=2,Npoints
  do j=1,i-1
    Uij(i,j)=Pair_LJ_NRG(r(:,i),r(:,j))
    Uij(j,i)=Uij(i,j)
  enddo
enddo
!write(*,*) 'test 3'
!==============================================================================!
!                           Generate QRG (simulated annealing)
!==============================================================================!
deltae1=0
T=T_i
counter=0
accept=0
reject=0
!open(25,File='rejects.dat')
do while (T > T_f)
  do i=1,T_iter
    counter=counter+1
    k=random_integer(1,Npoints)                               !Select Atom to Move
    call random_number(s)
    rr=r(:,k)+(2*s-1)                       !random number (0,1), make it (-1,1)
    if(V(rr).lt.E_cut) then                            !Only consider if V(trial)<Ecut
      U_move(k)=P(rr)
      Delta_E=0d0
      do j=1,Npoints
        if(j.ne.k) then
          U_move(j)=Pair_LJ_NRG(r(:,j),rr)
          Delta_E=Delta_E+Uij(j,k)-U_move(j)      !Energy change due to trial move
        endif
      enddo
      call random_number(t1)
      if(exp(Delta_E/T).ge.t1)then
        accept=accept+1
        Uij(:,k)=U_move(:)
        Uij(k,:)=U_move(:)
        r(:,k)=rr(:)
        deltae1=deltae1+Delta_E
      else
        reject=reject+1
!        write(25,*) T, reject
      endif
    endif
  enddo
  call T_schedule(cooling,T,cc)
enddo
!close(25)
out_cut = counter-(accept+reject)
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
write(99,*) 'Initial Temperature ==> ', T_i
write(99,*) 'Final Temperature ==> ', T_f
write(99,*) '# moves outside Ecut==> ', out_cut, float(out_cut)/float(counter)*100.
write(99,*) 'Iterations per Temperature ==> ', T_iter
write(99,*) 'Total number of iterations ==> ', counter
write(99,*) 'Total number of accepted moves ==> ', accept, float(accept)/float(counter-out_cut)*100.
write(99,*) 'Total number of rejected moves ==> ', reject, float(reject)/float(counter-out_cut)*100.
write(99,*) 'Simulation Time==> ',time2-time1
close(99)
end program main_grid
