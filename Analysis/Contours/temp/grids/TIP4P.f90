MODULE TIP4P_module

  !DOUBLE PRECISION, PARAMETER :: epsilon = 0.1852D0            ! kcal/mol
  DOUBLE PRECISION, PARAMETER :: bohr = 0.52917721092           ! atomic units
  DOUBLE PRECISION, PARAMETER :: autokJmol =2625.5002
  DOUBLE PRECISION, PARAMETER :: epsilon = 0.0002951349           ! atomic units
  DOUBLE PRECISION, PARAMETER :: sigma = 3.1589D0/bohr             ! bohr
  DOUBLE PRECISION, PARAMETER :: qM = 1.1128D0                 ! |e|
  DOUBLE PRECISION, PARAMETER :: gamma = 0.73612D0             ! (unit-less)
  !DOUBLE PRECISION, PARAMETER :: Dr = 116.09D0                 ! kcal/mol
  DOUBLE PRECISION, PARAMETER :: Dr = 0.1850012                 ! atomic units
  DOUBLE PRECISION, PARAMETER :: alphar = 2.287*bohr             ! bohr**(-1)
  !DOUBLE PRECISION, PARAMETER :: alphar = 2.287D0             ! angstrom**(-1)
  DOUBLE PRECISION, PARAMETER :: req =0.9419D0 /bohr                ! bohr
  !DOUBLE PRECISION, PARAMETER :: ktheta = 87.85D0              ! kcal/(mol*rad**2)
  DOUBLE PRECISION, PARAMETER :: ktheta = 0.139998              ! atomic units/(rad**2)
  !DOUBLE PRECISION, PARAMETER :: thetaeq = 107.4D0            ! degrees
  DOUBLE PRECISION, PARAMETER :: thetaeq = 1.87448361664D0      ! radians; 107.4D0*pi/180



CONTAINS


  !==============================================================================!


  SUBROUTINE TIP4P(NO, q, energy, force)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: NO                                      ! number of water molecules
    DOUBLE PRECISION, DIMENSION(3, 3*NO), INTENT(IN) :: q          ! coordinates; (Ox, Oy, Oz, H1x, H1y, H1z, H2x, H2y, H2z)
    DOUBLE PRECISION, DIMENSION(3, 3*NO), INTENT(INOUT) :: force
    DOUBLE PRECISION, INTENT(INOUT) :: energy

    DOUBLE PRECISION, DIMENSION(3*NO, 3) :: q_block
    DOUBLE PRECISION :: energyinter, energyintra
    DOUBLE PRECISION, DIMENSION(3*NO, 3) :: forceinter, forceintra


    q_block = TRANSPOSE(q)

    energyinter = 0.0d0
    energyintra = 0.0d0
    forceinter = 0.0d0
    forceintra = 0.0d0

!    WRITE(*,*) 'SIZE(force): ', SIZE(force)

    CALL TIP4P_inter(NO, q_block, energyinter, forceinter)
    CALL TIP4P_intra(NO, q_block, energyintra, forceintra)


       !WRITE(*,*) 'forceinter: '
       !print '(3F15.8)', (forceinter(i,:), i=1,3*NO)

       !WRITE(*,*) 'forceintra: '
       !print '(3F15.8)', (forceintra(i,:), i=1,3*NO)


    energy = energyinter + energyintra
    force = TRANSPOSE(forceinter + forceintra)


!    WRITE(*,*) 'SIZE(force): ', SIZE(force)   ! why does this throw a floating point exception; not anymore

!    WRITE(*,*) 'force: '
!    print '(3F15.8)', (force(:,i), i=1,3*NO)

  END SUBROUTINE TIP4P

  !==============================================================================!

  SUBROUTINE TIP4P_inter(NO, r, U_inter, F)

    IMPLICIT NONE

    INTEGER :: NO
    DOUBLE PRECISION, INTENT(INOUT) :: U_inter                    ! intermolecular potential
    DOUBLE PRECISION, DIMENSION(3*NO, 3), INTENT(INOUT) :: F      ! intermolecular forces
    DOUBLE PRECISION, DIMENSION(3*NO, 3), INTENT(INOUT) :: r      ! coordinates
    DOUBLE PRECISION, ALLOCATABLE :: q(:), rO(:,:), FLJ(:,:)      ! charges, O coordinates, LJ forces (on O)




    ALLOCATE(q(9*NO), rO(NO,3), FLJ(NO,3))

    q(1::3) = -qM
    q(2::3) = 0.5d0*qM
    q(3::3) = 0.5d0*qM

    rO = r(1::3,:)
    r(1::3,:) = gamma*r(1::3,:) + 0.5*(1.0-gamma)*(r(2::3,:) + r(3::3,:))

    F = 0d0
    U_inter = 0d0

    CALL Coulomb(NO, q, r, F, U_inter)

    F(2::3,:) = F(2::3,:) + 0.5*(1.0-gamma) * F(1::3,:)
    F(3::3,:) = F(3::3,:) + 0.5*(1.0-gamma) * F(1::3,:)
    F(1::3,:) = gamma * F(1::3,:)
    r(1::3,:) = rO

    CALL LJ(NO, r, F, U_inter)


  END SUBROUTINE TIP4P_inter



  SUBROUTINE LJ(NO, r, F, Utot)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NO 
    DOUBLE PRECISION, INTENT(in) :: r(:,:)
    DOUBLE PRECISION, INTENT(inout) :: F(:,:), Utot

    DOUBLE PRECISION, DIMENSION(NO) :: rsq, ir2, sir6, sir12, Cij
    DOUBLE PRECISION, DIMENSION(NO, 3) ::  dr, Fij
    INTEGER :: i, k

    DO i=1,NO-1
       DO k=1,3
          dr(i+1:NO,k) = r(3*i+1::3,k) - r(3*i-2,k)
       END DO
       rsq(i+1:NO) = dr(i+1:NO,1)**2 + dr(i+1:NO,2)**2 + dr(i+1:NO,3)**2
       ir2(i+1:NO) = 1.0/rsq(i+1:NO)


       sir6(i+1 : NO) = (sigma**2 * ir2(i+1:NO))**3
       sir12(i+1:NO) = sir6(i+1:NO)**2
       Utot  = Utot + 4.0 *epsilon* SUM(sir12(i+1:NO) - sir6(i+1:NO))

       Cij(i+1:NO) = 4.0*epsilon*(12.0 * sir12(i+1:NO)  - 6.0 * sir6(i+1:NO)) * ir2(i+1:NO)
       DO k=1,3
          Fij(i+1:NO,k) = dr(i+1:NO,k) * Cij(i+1:NO)
       END DO


       F(3*i-2,:) = F(3*i-2,:) - SUM(Fij(i+1:NO,:), 1)
       F(3*i+1::3,:) = F(3*i+1::3,:) +  Fij(i+1:NO,:)
    END DO
  END SUBROUTINE LJ


  SUBROUTINE Coulomb(NO, q, r, F, Utot)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NO 
    DOUBLE PRECISION, INTENT(in) :: q(:), r(:,:)
    DOUBLE PRECISION, INTENT(out) :: F(:,:), Utot

    INTEGER :: Natoms
    DOUBLE PRECISION, DIMENSION(3*NO) :: rsq, ir2, Cij
    DOUBLE PRECISION, DIMENSION(3*NO, 3) ::  dr, Fij
    INTEGER :: i, k


    Natoms = 3*NO
    DO i=1,Natoms-3
       DO k=1,3
          dr(i+1:Natoms,k) = r(i+1:Natoms,k) - r(i,k)
       END DO

       rsq(i+1:Natoms) = dr(i+1:Natoms,1)**2 + dr(i+1:Natoms,2)**2 + dr(i+1:Natoms,3)**2
       ir2(i+1:Natoms) = 1.0/rsq(i+1:Natoms)

       Cij(i+1 : Natoms) = q(i) * q(i+1:Natoms) * SQRT(ir2(i+1:Natoms))

       Cij(i : 3*((i + 2)/3)) = 0.0

       Utot  = Utot + SUM(Cij(i+1:Natoms))

       Cij(i+1:Natoms) = Cij(i+1:Natoms) * ir2(i+1:Natoms)
       DO k=1,3
          Fij(i+1:Natoms,k) = dr(i+1:Natoms,k) *Cij(i+1:Natoms)
       END DO

       F(i,:) = F(i,:) - SUM(Fij(i+1:Natoms,:), 1)
       F(i+1:Natoms,:) = F(i+1:Natoms,:) +  Fij(i+1:Natoms,:)
    END DO
  END SUBROUTINE Coulomb


  !==============================================================================!


  SUBROUTINE TIP4P_intra(NO, coords, U_intra, F_intra)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: NO                                            ! number of water molecules
    DOUBLE PRECISION, DIMENSION(3*NO, 3), INTENT(IN) :: coords           ! coordinates
    DOUBLE PRECISION, DIMENSION(3*NO, 3), INTENT(INOUT) :: F_intra       ! intramolecular forces
    DOUBLE PRECISION, INTENT(INOUT) :: U_intra              

    INTEGER :: i, k
    DOUBLE PRECISION :: U1, U2, Uth
    DOUBLE PRECISION :: r1, r2, r1dotr2, rtilde, theta, denom
    DOUBLE PRECISION, DIMENSION(3) :: r1vec, r2vec
    DOUBLE PRECISION, DIMENSION(3,3) :: q, F1, F2, Fth


    U1 = 0.0d0
    U2 = 0.0d0
    Uth = 0.0d0
    U_intra = 0.0d0
    F1 = 0.0d0
    F2 = 0.0d0
    Fth = 0.0d0

    DO i = 1, 3*(NO-1)+1, 3


       q = (coords(i:i+2, :))
       
       r1vec = q(2,:) - q(1,:)
       r2vec = q(3,:) - q(1,:)

       r1dotr2= r1vec(1)*r2vec(1) + r1vec(2)*r2vec(2) + r1vec(3)*r2vec(3)
       r1 = SQRT(r1vec(1)**2 + r1vec(2)**2 + r1vec(3)**2) 
       r2 = SQRT(r2vec(1)**2 + r2vec(2)**2 + r2vec(3)**2)

       rtilde = r1dotr2/(r1*r2)
       theta = ACOS(rtilde)
       denom = SQRT(1.0d0-rtilde**2)


       ! intramolecular potential energy
       U1 = (Dr * alphar**2 * (r1 - req)**2) * (1.0d0 - alphar*(r1-req)*(1.0d0 - (7.0d0/12.0d0)*alphar*(r1-req)))
       U2 = (Dr * alphar**2 * (r2 - req)**2) * (1.0d0 - alphar*(r2-req)*(1.0d0 - (7.0d0/12.0d0)*alphar*(r2-req)))
       Uth = (ktheta/2.0D0) * (theta - thetaeq)**2

       U_intra = U_intra + (U1 + U2 + Uth)


       !Gradient associated with O-H1 stretch
       F1(1,:) = (Dr*alphar**2*(r1-req)  *  (2.0d0 - 3.0d0*alphar*(r1-req) + (7.0d0/3.0d0) * alphar**2 * (r1-req)**2) )/r1
       F1(2,:) = F1(1,:)

       DO k = 1, 3
          F1(1, k) = F1(1, k) * (-(q(2,k)-q(1,k)))
          F1(2, k) = F1(2, k) * (q(2,k)-q(1,k))
       END DO


       !Gradient associated with O-H2 stretch
       F2(1,:) = (Dr*alphar**2*(r2-req)  *  (2.0d0 - 3.0d0*alphar*(r2-req) + (7.0d0/3.0d0) * alphar**2 * (r2-req)**2) )/r2
       F2(3,:) = F2(1,:)      

       DO k = 1, 3
          F2(1, k) = F2(1, k) * (-(q(3,k)-q(1,k)))
          F2(3, k) = F2(3, k) * (q(3,k)-q(1,k))
       END DO


       ! Gradient associated with H1-O-H2 bond angle
       DO k = 1, 3
          Fth(1, k) = (2.0d0*q(1,k)-q(2,k)-q(3,k))/(r1 * r2) &
               + (q(3,k) - q(1,k))*r1dotr2/(r1 * r2**3) &
               + (q(2,k) - q(1,k))*r1dotr2/(r1**3 * r2)

          Fth(2, k) = (q(3,k)-q(1,k))/(r1*r2) &
               - (q(2,k) - q(1,k))*r1dotr2/(r1**3 * r2)

          Fth(3, k) = (q(2,k)-q(1,k))/(r1*r2) &
               - (q(3,k) - q(1,k))*r1dotr2/(r1 * r2**3)
              
       END DO

       Fth = (-ktheta*(theta - thetaeq)/denom) * Fth

       ! intramolecular force
       F_intra(i:i+2, :) = -(F1 + F2 + Fth)

    END DO


  END SUBROUTINE TIP4P_intra


  !==============================================================================!


  SUBROUTINE test_grad(NO,q)
    ! gives the force using finite differences

    IMPLICIT NONE
    INTEGER :: i, NO
    REAL(8) :: s,Ener,Ener0
    REAL(8), DIMENSION(9*NO) :: F, F0, r, q

    s=1d-4

    CALL TIP4P(NO, q, Ener, F0)
    r=q
    DO i=1,9*NO
       r(i)=q(i)-s
       !   write(*,*) i
       CALL TIP4P(NO, r, Ener0, F)
       r(i)=q(i)+s
       !   write(*,*) i,'yo'
       CALL TIP4P(NO, r, Ener, F)
       !   write(*,*) i,'yo yo'
       F(i)=(Ener0-Ener)/(2*s)
       !   write(*,*) i,'yay'
       r(i)=q(i)
       WRITE(*,*) F(i), F0(i)
    ENDDO


  END SUBROUTINE test_grad


  !==============================================================================!


END MODULE TIP4P_module 
