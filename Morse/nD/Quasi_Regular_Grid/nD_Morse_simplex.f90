!=============================================================================80
!           nD-Morse Potential Quasi-Regular Grid (Simplex Minimization)
!==============================================================================!
!       Modified:
!   31 Janary 2020
!       Author:
!   Shane W. Flynn, Yang Liu, Vladimir A. Mandelshtam
!==============================================================================!
module QRG_mod
implicit none
!==============================================================================!
!                            Global Variables
!==============================================================================!
!d              ==>Particle Dimensionality
!Npoints        ==>Number of grid points
!c_LJ           ==>parameter for LJ
!E_cut          ==>Energy Cutoff Contour
!integral_P     ==>Normalization for P_x
!omega          ==>(d) Parameter for Morse Potential
!k              ==>Needed global for calling amoeba
!==============================================================================!
integer::d,Npoints,k
double precision,allocatable,dimension(:,:)::x
double precision,allocatable,dimension(:)::omega,xmin,xmax
double precision::integral_P,c_LJ,E_cut
!==============================================================================!
contains
!==============================================================================!
function random_integer(Nmin,Nmax)
!==============================================================================!
!Randomly generate an integer in the range Nmin-Nmax
!==============================================================================!
!Nmin           ==>minimum index value
!Nmax           ==>maximum index value
!random_integer ==>integer returned
!a              ==>random number (0,1)
!==============================================================================!
implicit none
integer::Nmin,Nmax,random_integer
double precision::a
call random_number(a)
random_integer=floor(a*(Nmax-Nmin+1))+Nmin
end function random_integer
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
!Target Distribution Function
!==============================================================================!
!P              ==>evaluate P(x)
!x              ==>(d) ith particles coordinate x^i_1,..,x^i_d
!==============================================================================!
function P(x)
implicit none
double precision::x(d),P
if(V(x)<E_cut) P=(E_cut-V(x))**(d/2.)/integral_P
if(V(x)>=E_cut) P=1d-20
end function P
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
!c_LJ           ==>parameter for LJ
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
function f(PP)
!==============================================================================!
!Computes the total potential energy (for simplex minimizaiton)
!==============================================================================!
implicit none
integer::i
double precision::PP(d),f,vtot
vtot=0.
do i=1,Npoints
  if(i.ne.k) vtot=vtot+Pair_LJ_NRG(x(:,i),PP(:))
enddo
f=vtot
end function f
!==============================================================================!
subroutine box_size(N_MMC_box)
!==============================================================================!
!Determine the box size for normalizing P; (xmin,xmax) using MMC
!==============================================================================!
!N_MMC_box      ==>Number of MMC Iterations to determine box-size
!mv_cutoff      ==>trial displacement move cutoff
!r              ==>(d) coordinates
!r_trial        ==>(d) trial coordinates
!s              ==>(d) trail displacement; random number for coordinates
!xmin           ==>(d) minimum of normalization box
!xmax           ==>(d) maximum of normalization box
!==============================================================================!
integer::N_MMC_box,i,j
double precision::dummy,r_trial(d),r(d),s(d)
double precision,parameter::mv_cutoff=0.1
r=0d0
xmin=r
xmax=r
do i=1,N_MMC_box
!==============================================================================!
!                   Generate coordinates for Random move
!           random numbers generated (0,1), make it (-1,1) ==> s=2*s-1
!==============================================================================!
  call random_number(s)
  r_trial=r+mv_cutoff*(2*s-1)
!==============================================================================!
!                             Test Acceptance
!==============================================================================!
  call random_number(dummy)
  if(P(r_trial)/P(r).ge.dummy) then
    r=r_trial
    do j=1,d
      if(xmin(j)>r(j)) xmin(j)=r(j)
      if(xmax(j)<r(j)) xmax(j)=r(j)
    enddo
  endif
enddo
end subroutine box_size
!==============================================================================!
subroutine compute_integral_P(N_1D)
!==============================================================================!
!Use a (uniform) square grid to integrate P(r)
!Box size for the grid is determined in the box_size subroutine
!int P(r)~Area_Square/N sum_n=1,N P(r_n)
!==============================================================================!
!N_1D           ==> Number of points in a single dimension for the integral
!Ntotal         ==>Total number of evaluations for all dimensions
!xmin           ==>(d) minimum of normalization box
!xmax           ==>(d) maximum of normalization box
!Moment         ==>First Moments for the distribution
!r(d)           ==>Coordinates
!==============================================================================!
integer(kind=8)::N_1D
integer::i,j,index1(d),Ntotal
double precision::r(d),Moment,dummy,delx(d)
Moment=0.
Ntotal=(N_1D+1)**d
index1=0
delx(:)=(xmax(:)-xmin(:))/N_1D
do i=1,Ntotal
  do j=1,d
    if(index1(j).eq.N_1D) then
      index1(j)=0
    else
      index1(j)=index1(j)+1
      exit
    endif
  enddo
  r(:)=xmin(:)+index1(:)*delx(:)
  dummy=P(r)
  Moment=Moment+dummy
enddo
dummy=1./N_1D**d
do j=1,d
  dummy=dummy*(xmax(j)-xmin(j))
enddo
integral_P=dummy*Moment
write(*,*) 'integral_P ==> ', integral_P
end subroutine compute_integral_P
!==============================================================================!
!We have adapted the original amoeba subroutine made available by J-P Moreau
!Special thanks for making this code available!
!   http://jean-pierre.moreau.pagesperso-orange.fr/Fortran/tamoeba_f90.txt
!See Numerical Recipies for more details on downhill simplex minimization
!==============================================================================!
  SUBROUTINE AMOEBA(P,Y,Npar,FTOL,ITER)
  !-------------------------------------------------------------------
  ! Multidimensional minimization of the function f(X) where X is
  ! an Npar-dimensional vector, by the downhill simplex method of
  ! Nelder and Mead. Input is a matrix P whose Npar+1 rows are Npar-
  ! dimensional vectors which are the vertices of the starting simplex
  ! (Logical dimensions of P are P(Npar+1,Npar); physical dimensions
  ! are input as P(Npar+1,Npar)). Also input is the vector Y of length Npar
  ! +1, whose components must be pre-initialized to the values of f
  ! evaluated at the Npar+1 vertices (rows) of P; and FTOL the fractio-
  ! nal convergence tolerance to be achieved in the function value. On
  ! output, P and Y will have been reset to Npar+1 new points all within
  ! FTOL of a minimum function value, and ITER gives the number of ite-
  ! rations taken.
!==============================================================================!
!       on input ITER is the maximum number of iterations to be perform
!           on output it is the actual number of iterations
!==============================================================================!
    implicit none
    double precision, PARAMETER :: ALPHA=1.d0,BETA=0.5d0,GAMMA=2.d0
    ! Expected maximum number of dimensions, three parameters which define
    ! the expansions and contractions, and maximum allowed number of
    ! iterations.
    integer Npar, ITER, MPTS, ILO, IHI, INHI, I, J, ITMAX
    double precision ::  P(Npar+1,Npar), Y(Npar+1), PR(Npar), PRR(Npar), PBAR(Npar)
    double precision :: FTOL, RTOL, YPR, YPRR

    ITMAX=ITER
    MPTS=Npar+1
    ITER=0
  1 ILO=1
    IF(Y(1).GT.Y(2)) THEN
       IHI=1
       INHI=2
    ELSE
       IHI=2
       INHI=1
    ENDIF
    DO I=1, MPTS
       IF(Y(I).LT.Y(ILO)) ILO=I
       IF(Y(I).GT.Y(IHI)) THEN
          INHI=IHI
          IHI=I
       ELSE IF (Y(I).GT.Y(INHI)) THEN
          IF(I.NE.IHI) INHI=I
       END IF
    END DO
    ! Compute the fractional range from highest to lowest and return if
    ! satisfactory.
    RTOL=2.d0*ABS(Y(IHI)-Y(ILO))/(ABS(Y(IHI))+ABS(Y(ILO)))
    IF(RTOL.LT.FTOL) RETURN
    IF(ITER.EQ.ITMAX) then
       write(*,*) ' Amoeba exceeding maximum iterations.'
       return
    endif
    ITER=ITER+1
    DO J=1, Npar
       PBAR(J)=0.d0
    END DO
    DO I=1, MPTS
       IF(I.NE.IHI) THEN
          DO J=1,Npar
             PBAR(J)=PBAR(J) + P(I,J)
          END DO
       END IF
    END DO
    DO J=1, Npar
       PBAR(J)=PBAR(J)/Npar
       PR(J)=(1.d0+ALPHA)*PBAR(J) - ALPHA*P(IHI,J)
    END DO
    YPR=f(PR)
    IF(YPR.LE.Y(ILO)) THEN
       DO J=1,Npar
          PRR(J)=GAMMA*PR(J) + (1.d0-GAMMA)*PBAR(J)
       END DO
       YPRR=f(PRR)
       IF(YPRR.LT.Y(ILO)) THEN
          DO J=1, Npar
             P(IHI,J)=PRR(J)
          END DO
          Y(IHI)=YPRR
       ELSE
          DO J=1, Npar
             P(IHI,J)=PR(J)
          END DO
          Y(IHI)=YPR
       END IF
    ELSE IF(YPR.GE.Y(INHI)) THEN
       IF(YPR.LT.Y(IHI)) THEN
          DO J=1, Npar
             P(IHI,J)=PR(J)
          END DO
          Y(IHI)=YPR
       END IF
       DO J=1, Npar
          PRR(J)=BETA*P(IHI,J) + (1.d0-BETA)*PBAR(J)
       END DO
       YPRR=f(PRR)
       IF(YPRR.LT.Y(IHI)) THEN
          DO J=1, Npar
             P(IHI,J)=PRR(J)
          END DO
          Y(IHI)=YPRR
       ELSE
          DO I=1, MPTS
             IF(I.NE.ILO) THEN
                DO J=1,Npar
                   PR(J)=0.5d0*(P(I,J) + P(ILO,J))
                   P(I,J)=PR(J)
                END DO
                Y(I)=f(PR)
             END IF
          END DO
       END IF
    ELSE
       DO J=1, Npar
          P(IHI,J)=PR(J)
       END DO
       Y(IHI)=YPR
    END IF
    GO TO 1
  END SUBROUTINE AMOEBA
!==============================================================================!
end module QRG_mod
!==============================================================================!
program main
use QRG_mod
!==============================================================================!
!update names here
!==============================================================================!
implicit none
integer::N_MMC_box,N_MMC_grid,ITER,Npar,m,i,j
integer(kind=8)::N_1D
double precision::FTOL
double precision,allocatable,dimension(:)::x0,s,U_move,Y
double precision,allocatable,dimension(:,:)::U,PP
!==============================================================================!
!                               Allocations
!==============================================================================!
read(*,*) d
allocate(omega(d))
read(*,*) omega
read(*,*) Npoints
read(*,*) N_MMC_box
read(*,*) N_1D
read(*,*) N_MMC_grid
read(*,*) E_cut
read(*,*) c_LJ
read(*,*) ITER
read(*,*) FTOL
!==============================================================================!
!                               Allocations
!==============================================================================!
Npar=d                                                   !For simplex subroutine
allocate(x(d,Npoints),x0(d),s(d),U(Npoints,Npoints),U_move(Npoints),xmin(d))
allocate(xmax(d),PP(Npar+1,Npar),Y(Npar+1))
write(*,*) 'Test 0; Successfully Read Input File'
!==============================================================================!
!                Run Monte Carlo to Normalize P/Get Moments
!==============================================================================!
integral_P=1d0                       !set equal to 1 so you can initially call P
call box_size(N_MMC_box)
do i=1,d
  write(*,*) i,' xmin xmax ==>', xmin(i),xmax(i)
enddo
xmin=xmin-(xmax-xmin)*0.01
xmax=xmax+(xmax-xmin)*0.01
call compute_integral_P(N_1D)
write(*,*) 'Test 1; Successfully normalized P(x)'
!==============================================================================!
!                       Generate Initial Distribution
!               Initally accept any point where Potential<Ecut
!==============================================================================!
i=1
do while(i.le.Npoints)
    call random_number(s)
    s(:)=xmin(:)+s(:)*(xmax(:)-xmin(:))
    if(V(s)<E_cut)then
        x(:,i)=s(:)
        i=i+1
    endif
enddo
!==============================================================================!
!                          Write Initial Coordinates
!==============================================================================!
open(unit=17,file='coor_ini.dat')
do i=1,Npoints
    write(17,*) x(:,i)
enddo
close(17)
!==============================================================================!
!                       Begin Simplex to Optimize GridPoints
!==============================================================================!
!                           Select Atom to Move
!==============================================================================!
do i=1,N_MMC_grid
  k=random_integer(1,Npoints)
  x0=x(:,k)
  PP(1,:)=x0(:)
  do j=2,Npar+1
    PP(j,j-1)=PP(1,j-1)+0.05
  enddo
  do m=1,Npar+1
    Y(m)=f(PP(m,:))
  enddo
  CALL AMOEBA(PP,Y,Npar,FTOL,ITER)
!==============================================================================!
!                             select the minimum point
!==============================================================================!
  if(Y(1).eq.min(Y(1),Y(2),Y(3))) then
    x(:,k)=PP(1,:)
  else if(Y(2).eq.min(Y(1),Y(2),Y(3))) then
    x(:,k)=PP(2,:)
  else
    x(:,k)=PP(3,:)
  endif
enddo
!==============================================================================!
!                               Write Optimized Grid
!==============================================================================!
open(unit=20,file='grid.dat')
do i=1,Npoints
  write(20,*) x(:,i)
enddo
close(20)
write(*,*) 'Successfully Generated Quasi-Regular Grid'
end program main
