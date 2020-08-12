!=============================================================================80
!                       QRG Morse Implementation                  !
!==============================================================================!
!implement move cutoff that scales wrt sigma.
!writes grid to screen at each move cutoff update
!initial distribution made using rejection.
!==============================================================================!
!    Modified:
!20 August 2020
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
double precision,allocatable,dimension(:)::omega,rmin,rmax,x0
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
a=sum(abs(x1(:)-x2(:)))
sigma1=c_lj*(P(x1)*Npoints)**(-1./d)
sigma2=c_LJ*(P(x2)*Npoints)**(-1./d)
Pair_LJ_NRG=(sigma1/a)**(9+d)-(sigma1/a)**(3+d) + &
            (sigma2/a)**(9+d)-(sigma2/a)**(3+d)
end function Pair_LJ_NRG
!==============================================================================!
function sigma_local(x1)
!==============================================================================!
!sigma for a given point
!==============================================================================!
!x1             ==>(d) ith atoms coordinates
!a/b            ==>evaluate q-LJ
!sigma          ==>Gaussian widths (c*sigma(P))
!Pair_LJ_NRG    ==>Energy of the i-j q-LJ potential
!==============================================================================!
implicit none
double precision::x1(d),sigma_local
sigma_local=c_LJ*(P(x1)*Npoints)**(-1./d)
end function sigma_local
!==============================================================================!
end module morse_mod
!==============================================================================!
program main_grid
use morse_mod
!==============================================================================!
implicit none
integer::N_MMC_box,N_1D,N_MMC_grid,MMC_freq,accept,counter,i,j,k,l,Ntotal
integer::total_accept,total_reject,out_cut,my_size
integer,allocatable::seed(:)
double precision::time1,time2,Delta_E,dummy,moment,mv_cutoff,reject,P_max
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
read(*,*) E_cut
read(*,*) c_LJ
read(*,*) N_1D
read(*,*) N_MMC_box
read(*,*) N_MMC_grid
read(*,*) MMC_freq
!==============================================================================!
allocate(r(d,Npoints),Uij(Npoints,Npoints),r_i(d),delr(d),index1(d),x0(d))
allocate(U_move(Npoints),rr(d),r_trial(d),s(d))
!==============================================================================!
!Set seed for random number generator (make it reproducible)
!==============================================================================!
call random_seed(size=my_size)
allocate(seed(my_size))
seed=my_size+1          !for gfortran seed must be larger than size, simple fix
call random_seed(put=seed)
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
  r_trial=r_i+mv_cutoff*(2*s-1) !Make trial move (-1,1), not (0,1)
  call random_number(dummy)                             !MMC acceptance criteria
  if(P(r_trial)/P(r_i).ge.dummy) then
    r_i=r_trial
    do j=1,d
      if(rmin(j).gt.r_i(j)) rmin(j)=r_i(j)
      if(rmax(j).lt.r_i(j)) rmax(j)=r_i(j)
    enddo
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
!==============================================================================!
!           Generate initial distribution to then convert to a QRG
!Use pseudo-random number generator and rejection to make the initial dist.
!==============================================================================!
x0=0                                                  !for morse minimum is at 0
P_max=P(x0)
i=1
do while(i.le.Npoints)
  call random_number(r_i)
  r_i(:)=rmin(:)+r_i(:)*(rmax(:)-rmin(:)) !scalue unif. dist. to span domain
  call random_number(reject)
  if(P(r_i)/P_max.gt.reject)then
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
!==============================================================================!
!                             Generate QRG (greedy search)
!==============================================================================!
accept=0
counter=0
mv_cutoff=1.
out_cut=0
total_accept=0
total_reject=0
open(77,File='trajectory.dat')
do i=1,N_MMC_grid
  k=random_integer(1,Npoints)                               !Select Atom to Move
  call random_number(s)
  rr=r(:,k)+mv_cutoff*sigma_local(r(:,k))*(2*s-1)               !random number (0,1), make it (-1,1)
  if(V(rr).lt.E_cut) then                            !Only consider if V(trial)<Ecut
    counter=counter+1
    U_move(k)=P(rr)
    Delta_E=0d0
    do j=1,Npoints
      if(j.ne.k) then
        U_move(j)=Pair_LJ_NRG(r(:,j),rr)
        Delta_E=Delta_E+Uij(j,k)-U_move(j)      !Energy change due to trial move
      endif
    enddo
    if(Delta_E.ge.0d0)then
      total_accept=total_accept+1
      Uij(:,k)=U_move(:)
      Uij(k,:)=U_move(:)
      accept=accept+1
      r(:,k)=rr(:)
    else
      total_reject=total_reject+1
    endif
  else
    out_cut=out_cut+1
  endif
!for MMC want acceptance ~30-50%, adjust trial movement displacement accordingly
  if(mod(i,MMC_freq)==0)then
    if(dble(accept)/counter.lt.0.4)then
      mv_cutoff=mv_cutoff*0.9
    else
      mv_cutoff=mv_cutoff*1.1
    endif
  accept=0
  counter=0
  write(77,*)
  do l=1,Npoints
    write(77,*) r(:,l)
  enddo
  endif
enddo
close(77)
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
write(99,*) '# moves outside Ecut==> ', out_cut, float(out_cut)/float(N_MMC_grid)*100.
write(99,*) 'Total number of accepted moves ==> ', total_accept, &
  float(total_accept)/float(N_MMC_grid-out_cut)*100.
write(99,*) 'Total number of rejected moves ==> ', total_reject, &
float(total_reject)/float(N_MMC_grid-out_cut)*100.
write(99,*) 'Simulation Time==> ',time2-time1
write(99,*) 'my_size==> ',my_size
write(99,*) 'seed==> ',seed(1)
write(99,*) 'x0==> ',x0
write(99,*) 'P_max==> ',P_max
close(99)
end program main_grid
