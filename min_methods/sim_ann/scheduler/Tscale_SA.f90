!=============================================================================80
!                   Simulated Annealing Temperature Scales
!==============================================================================!
!simple implementation for testing different Temperature scales for applying
!simulated annealing. Output is each Temperature cycle for a given set of values
!==============================================================================!
!    Modified:
!3 August 2020
!    Author:
!Shane Flynn
!==============================================================================!
module SA_mod
implicit none
!==============================================================================!
!cooling      ==>Method for scheduling cooling cycles
!T            ==>Temperature
!cc           ==>cooling constant (Assumed to be > 0)
!T_i          ==>Initial Temperature
!T_f          ==>Final Temperature
!T_iters      ==>Iterations per temperature cycle
!counter      ==>Total number of iterations in simulation
!==============================================================================!
contains
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
character(len=50)::cooling
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
character(len=50)::cooling
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
character(len=50)::cooling
!==============================================================================!
T=T/(1+cc*T)
end subroutine log_cooling
!==============================================================================!
end module SA_mod
!==============================================================================!
program main_grid
use SA_mod
!==============================================================================!
implicit none
character(len=50)::cooling
double precision::T,T_i,T_f,cc
integer::T_iters,i,counter
!==============================================================================!
!                              Read Input File                                 !
!==============================================================================!
read(*,*) T_i
read(*,*) T_f
read(*,*) T_iters
read(*,*) cc
read(*,*) cooling
!==============================================================================!
T=T_i
counter=0
open(27,file='temperature.dat')
do while (T>T_f)
  do i=1,T_iters
    if(i==1) write(27,*) T
    counter=counter+1
  enddo
  call T_schedule(cooling,T,cc)
enddo
close(27)
!==============================================================================!
!                                 Output
!==============================================================================!
open(99,file='out.txt')
write(99,*) 'Initial Temperature ==> ', T_i
write(99,*) 'Final Temperature ==> ', T_f
write(99,*) '# of Iterations per Temperature ==> ', T_iters
write(99,*) 'Total # of iterations ==> ', counter
close(99)
end program main_grid
