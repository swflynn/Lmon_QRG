!=================20==================40==================60==================80
!Test orthonormality of Smat and Hmat
!==============================================================================!
program main
!==============================================================================!
implicit none
integer::NG,Lwork,info,i,j,itype
double precision::dummy1,dummy2,dummy3,dummy4
double precision,allocatable,dimension(:)::work
double precision,allocatable,dimension(:,:)::Smat,Sev,Hmat,Hev
!==============================================================================!
NG=100
allocate(Smat(NG,NG),Sev(NG,NG),Hmat(NG,NG),Hev(NG,NG))
!==============================================================================!
open(unit=56,file='overlap_matrix.dat')
open(unit=57,file='hamiltonian_matrix.dat')
do i=1,NG
  do j=i,NG
    read(56,*) Smat(i,j)
    read(57,*) Hmat(i,j)
    Smat(j,i)=Smat(i,j)
    Hmat(j,i)=Hmat(i,j)
  enddo
enddo
!==============================================================================!
dummy1=0
dummy2=0
dummy3=0
dummy4=0
do i=1,NG
  do j=1,NG
    dummy1=dummy1+(Hmat(1,i)*Smat(i,j)*Hmat(1,j))
    dummy2=dummy2+(Hmat(i,1)*Smat(i,j)*Hmat(j,1))
    dummy3=dummy3+(Hmat(2,i)*Smat(i,j)*Hmat(1,j))
    dummy4=dummy4+(Hmat(i,2)*Smat(i,j)*Hmat(j,1))
  enddo
enddo
write(*,*) 'Test 1 ==>'
write(*,*) 'dummy1 test: ', dummy1
write(*,*) 'dummy2 test: ', dummy2
write(*,*) 'dummy3 test: ', dummy3
write(*,*) 'dummy4 test: ', dummy4
!==============================================================================!
!check eigenvectors
open(unit=58,file='overlap_eigenvectors.dat')
open(unit=59,file='hamiltonian_eigenvectors.dat')
do i=1,NG
  do j=1,NG
    read(58,*) Sev(i,j)
    read(59,*) Hev(i,j)
  enddo
enddo
close(58)
close(59)
dummy1=0
dummy2=0
dummy3=0
dummy4=0
do i=1,NG
  do j=1,NG
    dummy1=dummy1+(Hev(1,i)*Sev(i,j)*Hev(1,j))
    dummy2=dummy2+(Hev(i,1)*Sev(i,j)*Hev(j,1))
    dummy3=dummy3+(Hev(2,i)*Sev(i,j)*Hev(1,j))
    dummy4=dummy4+(Hev(i,2)*Sev(i,j)*Hev(j,1))
  enddo
enddo
write(*,*) 'Test 2 ==>'
write(*,*) 'dummy1 test: ', dummy1
write(*,*) 'dummy2 test: ', dummy2
write(*,*) 'dummy3 test: ', dummy3
write(*,*) 'dummy4 test: ', dummy4
!==============================================================================!
dummy1=0
dummy2=0
dummy3=0
dummy4=0
do i=1,NG
  do j=1,NG
    dummy1=dummy1+(Hev(1,i)*Smat(i,j)*Hev(1,j))
    dummy2=dummy2+(Hev(2,i)*Smat(i,j)*Hev(1,j))
    dummy3=dummy3+(Hev(i,1)*Smat(i,j)*Hev(j,1))
    dummy4=dummy4+(Hev(i,2)*Smat(i,j)*Hev(j,1))
  enddo
enddo
write(*,*) 'Test 3 ==>'
write(*,*) 'dummy1 test: ', dummy1
write(*,*) 'dummy2 test: ', dummy2
write(*,*) 'dummy3 test: ', dummy3
write(*,*) 'dummy4 test: ', dummy4
!==============================================================================!
end program main
