!  Filename : realtime3d.f90
!
!  Title : Fortran program for time-dependent Gross-Pitaevskii equation in 
!          a fully anisotropic trap
!  Authors : P. Muruganandam and S. K. Adhikari
! 
!  Fortran program for Gross-Pitaevskii equation in three-dimensional 
!  anisotropic trap by real time propagation (Fortran 90/95 Version)
!
MODULE COMM_DATA
! NX, NY, NZ : Number of space mesh points (X, Y and Z)  
  INTEGER, PARAMETER :: NX = 200, NXX = NX-1, NX2 = NX/2
  INTEGER, PARAMETER :: NY = 160, NYY = NY-1, NY2 = NY/2
  INTEGER, PARAMETER :: NZ = 120, NZZ = NZ-1, NZ2 = NZ/2
! NSTP : Number of iterations to introduce the nonlinearity
! NPAS : Number of subsequent iterations with fixed nonlinearity, and 
! NRUN : Number of final time iterations with fixed nonlinearity
  INTEGER, PARAMETER ::  NSTP = 60000, NPAS = 1000, NRUN = 4000
  REAL (8), PARAMETER :: PI = 3.14159265358979D0
END MODULE COMM_DATA

MODULE GPE_DATA
  USE COMM_DATA, ONLY : NX, NY, NZ
! DX, DY, DZ : Space step and DT : Time step 
  REAL (8), PARAMETER :: DX = 0.1D0, DY = 0.1D0, DZ = 0.1D0
  REAL (8), PARAMETER :: DT = 0.002D0 
! G0 : Final nonlinearity
  REAL (8), PARAMETER :: G0 = 22.454D0 ! Final nonlinearity
  COMPLEX (8), PARAMETER :: CI = (0.0D0,1.0D0)  ! Complex i
! OPTION and XOP decides which equation to be solved.
  ! OPTION=1 Solves -psi_xx-psi_yy-psi_zz+V(x,y,z)psi+G0|psi|^2 psi=i psi_t
  ! OPTION=2 Solves [-psi_xx-psi_yy-psi_zz+V(x,y,z)psi]/2+G0|psi|^2 psi=i psi_t
  ! OPTION=3 Solves -psi_xx-psi_yy-psi_zz+V(x,y,z)psi/4+G0|psi|^2 psi=i psi_t
  INTEGER, PARAMETER :: OPTION = 2
! AL, BL : Anisotropy coefficients
  REAL (8) :: AL, BL, G, XOP
! X(0:NX), Y(0:NY), Z(0:NZ) : Space mesh, V(0:NX,0:NY,0:NZ) : Potential 
! CP(0:NX,0:NY,0:NZ) : Wave function (COMPLEX)
  REAL (8), DIMENSION(0:NX,0:NY,0:NZ) :: V, R2
  REAL (8), DIMENSION(0:NX) :: X, X2
  REAL (8), DIMENSION(0:NY) :: Y, Y2
  REAL (8), DIMENSION(0:NZ) :: Z, Z2
  COMPLEX (8), DIMENSION(0:NX,0:NY,0:NZ) :: CP
END MODULE GPE_DATA

MODULE CN_DATA
  USE COMM_DATA, ONLY : NX, NY, NZ
  COMPLEX (8), DIMENSION(0:NX) :: CALA, CGAA
  COMPLEX (8), DIMENSION(0:NY) :: CALB, CGAB
  COMPLEX (8), DIMENSION(0:NZ) :: CALC, CGAC
  COMPLEX (8) :: CT0X, CT0Y, CT0Z
  COMPLEX (8) :: CA0, CB0, CC0, CA0R, CB0R, CC0R
END MODULE CN_DATA 

MODULE TEMP_DATA
  USE COMM_DATA, ONLY : NX, NY, NZ
  COMPLEX (8), DIMENSION(0:NX,0:NY,0:NZ) :: CBE
  REAL (8), DIMENSION(0:NX,0:NY,0:NZ) :: P, TMP3
  REAL (8), DIMENSION(0:NX,0:NY) :: TMP2
  REAL (8), DIMENSION(0:NX) :: TMP1
END MODULE TEMP_DATA

PROGRAM GROSS_PITAEVSKII_SSCN_3D
  USE COMM_DATA
  USE GPE_DATA
  IMPLICIT NONE
!------------------------ INTERFACE BLOCKS -----------------------
  INTERFACE 
    SUBROUTINE INITIALIZE()
      IMPLICIT NONE
    END SUBROUTINE INITIALIZE
  END INTERFACE

  INTERFACE 
    SUBROUTINE COEF()
      IMPLICIT NONE
    END SUBROUTINE COEF
  END INTERFACE

  INTERFACE 
    SUBROUTINE NU(CP, DT)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: CP
      REAL (8), INTENT(IN) :: DT
    END SUBROUTINE NU
  END INTERFACE

  INTERFACE 
    SUBROUTINE LUX(CP)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: CP
    END SUBROUTINE LUX
  END INTERFACE

  INTERFACE 
    SUBROUTINE LUY(CP)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: CP
    END SUBROUTINE LUY
  END INTERFACE

  INTERFACE 
    SUBROUTINE LUZ(CP)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: CP
    END SUBROUTINE LUZ
  END INTERFACE

  INTERFACE
    SUBROUTINE NORM(CP, ZNORM)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: CP
      REAL (8), INTENT(OUT) :: ZNORM
    END  SUBROUTINE NORM
  END INTERFACE

  INTERFACE
    SUBROUTINE RAD(CP, RMS)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:,0:,0:), INTENT(IN) :: CP
      REAL (8), INTENT(OUT) :: RMS
    END  SUBROUTINE RAD
  END INTERFACE

  INTERFACE
    SUBROUTINE CHEM(CP, MU, EN)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:,0:,0:), INTENT(IN) :: CP
      REAL (8), INTENT(OUT) :: MU, EN
    END  SUBROUTINE CHEM
  END INTERFACE
!----------------------- END INTERFACE BLOCKS --------------------
  INTEGER :: K, I, J
  REAL (8) :: ZNORM, MU, EN, RMS, GSTP
  SELECT CASE (OPTION)
    CASE (1,3)
      XOP = 1.0D0
    CASE (2)
      XOP = 2.0D0
    CASE (:0,4:)
      PRINT *, 'ERROR: Wrong option', OPTION
      STOP
  END SELECT
  AL = SQRT(2.0D0)      !  Anisotropy 'LAMBDA' of the trap
  BL = 2.0D0            !  Anisotropy 'KAPPA' of the trap
  WRITE(7,900) OPTION
  WRITE(7,901) AL, BL
  WRITE(7,*)  
  WRITE(7,902) NX, NY, NZ
  WRITE(7,903) NSTP, NPAS, NRUN  
  WRITE(7,904) G0 
  WRITE(7,905) DX, DY, DZ
  WRITE(7,906) DT 
  WRITE(7,*)
  900 FORMAT('  OPTION = ',I3)
  901 FORMAT('  Anisotropy AL = ',F12.6,', BL = ',F12.6)
  902 FORMAT('# Space Stp NX = ',I8,', NY = ',I8,', NZ = ',I8)
  903 FORMAT('# Time Stp : NSTP = ',I9,', NPAS = ',I9,', NRUN = ',I9)
  904 FORMAT('  Nonlinearity G = ', F16.8)
  905 FORMAT('  Space Step DX = ',F10.6,', DY = ',F10.6,', DZ = ',F10.6)
  906 FORMAT('  Time Step  DT = ', F10.6)
  G = 0.0D0
! INITIALIZE() initializes the starting normalized wave function 'CP' and
! harmonic oscillator potential V
  CALL INITIALIZE()
! COEF() defines the coefficients of Crank-Nicholson Scheme.
  CALL COEF()
! NORM() calculates norm and restores normalization
  CALL NORM(CP, ZNORM)
! RAD() calculates the r.m.s radius RMS
  CALL RAD(CP, RMS)
! CHEM() calculates the chemical potential MU and energy EN.
  CALL CHEM(CP, MU, EN)
  999 FORMAT (2F12.6, F16.8)
  WRITE (7, 1001)
  WRITE (7, 1002)
  WRITE (7, 1001)
  WRITE (7, 1003) ZNORM, MU/XOP, EN/XOP, RMS, ABS(CP(NX2, NY2, NZ2))
  1001 FORMAT (18X, 54('-'))
  1002 FORMAT (20X, 'Norm', 6X, 'Chem', 8X, 'Ener', 8X, '<r>', 5X, 'Psi(0,0,0)')
  1003 FORMAT ('Initial :  ', 3X, F11.3, 2F12.3, 3F11.3)


            DO I = 0, NX
                WRITE(1,1000)X(I),ABS(CP(I,NY/2,NZ/2))
             END DO
             DO J = 0, NY
                WRITE(2,1000)Y(J),ABS(CP(NX/2,J,NZ/2))
             END DO
             DO K = 0, NZ
                WRITE(3,1000)Z(K),ABS(CP(NX/2,NY/2,K))
             END DO


  GSTP = XOP*G0/NSTP
  DO K = 1, NSTP ! Introduces nonlinearity in 'NSTP' iterations
     G = G + GSTP
     CALL NU(CP, DT)
     CALL LUX(CP)
     CALL LUY(CP)
     CALL LUZ(CP)
  END DO
  CALL NORM(CP, ZNORM)
  CALL CHEM(CP, MU, EN)
  CALL RAD(CP, RMS)
  WRITE (7, 1004) ZNORM, MU/XOP, EN/XOP, RMS, ABS(CP(NX2, NY2, NZ2))
  1004 FORMAT('After NSTP iter.:', F8.3, 2F12.3, 3F11.3)
!
  DO K = 1, NPAS 
     CALL NU(CP, DT)
     CALL LUX(CP)
     CALL LUY(CP)
     CALL LUZ(CP)
  END DO
  CALL NORM(CP, ZNORM)
  CALL CHEM(CP, MU, EN)
  CALL RAD(CP, RMS)
  WRITE (7, 1005) ZNORM, MU/XOP, EN/XOP, RMS, ABS(CP(NX2, NY2, NZ2))
  WRITE (7, 1001)

 1000 FORMAT(F12.4,F16.6)

  1005 FORMAT('After NPAS iter.:', F8.3, 2F12.3, 3F11.3)

            DO I = 0, NX
                WRITE(11,1000)X(I),ABS(CP(I,NY/2,NZ/2))
             END DO
             DO J = 0, NY
                WRITE(12,1000)Y(J),ABS(CP(NX/2,J,NZ/2))
             END DO
             DO K = 0, NZ
                WRITE(13,1000)Z(K),ABS(CP(NX/2,NY/2,K))
             END DO


!
!  The following line defines a nonstationary problem which is studied
!  below and the time evolution written on file 8 via WRITE(8,*)
  G = 0.75*G
  DO K = 1, NRUN 
     CALL NU(CP, DT)
     CALL LUX(CP)
     CALL LUY(CP)
     CALL LUZ(CP)
     IF (MOD(K,4).EQ.0) THEN
        CALL RAD(CP, RMS)
        WRITE(8,999) K*DT*XOP, RMS 
     END IF
  END DO
END PROGRAM GROSS_PITAEVSKII_SSCN_3D

SUBROUTINE INITIALIZE()
!   Routine that initizlizes the constant and variables.
!   Calculates the potential term V and the initial wave function CP
  USE COMM_DATA, ONLY : NX, NY, NZ, NX2, NY2, NZ2, PI
  USE GPE_DATA
  USE TEMP_DATA, ONLY : TMP3
  IMPLICIT NONE
  REAL (8) :: SPI, SP, AL2, BL2
  INTEGER :: I, J, K
  SP = SQRT(PI*SQRT(PI/(AL*BL)))
  SPI = SP*SQRT(2.0D0*SQRT(2.0D0))
  AL2 = AL*AL
  BL2 = BL*BL
  FORALL (I=0:NX) X(I) = (I-NX2)*DX
  FORALL (J=0:NY) Y(J) = (J-NY2)*DY
  FORALL (K=0:NZ) Z(K) = (K-NZ2)*DZ
  X2 = X*X
  Y2 = Y*Y
  Z2 = Z*Z
  FORALL (J=0:NY, K=0:NZ)
    R2(:,J,K) = X2 + Y2(J) + Z2(K)
    V(:,J,K) = X2 + AL2*Y2(J) + BL2*Z2(K)
    TMP3(:,J,K) = X2 + AL*Y2(J) + BL*Z2(K)
  END FORALL
  SELECT CASE (OPTION)
    CASE (1,2)
       CP = EXP(-TMP3/2.0D0)/SP
    CASE (3)
       V = V/2.0D0
       CP = EXP(-TMP3/4.0D0)/SPI
  END SELECT
END SUBROUTINE INITIALIZE
!
SUBROUTINE COEF()
!   Calculates the coefficients needed in subroutines LUX, LUY, LUZ      
  USE COMM_DATA, ONLY : NX, NY, NZ, NXX, NYY, NZZ
  USE GPE_DATA, ONLY : DX, DY, DZ, DT, CI
  USE CN_DATA
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL (8) :: DX2, DY2, DZ2
  REAL (8) :: DXX, DYY, DZZ
  COMPLEX (8) :: CDT
!
  DX2 = DX*DX  ! Generating the Coefficients for the
  DY2 = DY*DY  ! C-N method (to solve the spatial part)
  DZ2 = DZ*DZ
!
  DXX = 1.0D0/DX2
  DYY = 1.0D0/DY2
  DZZ = 1.0D0/DZ2
  CDT = CI*DT
!
  CA0 = 1.0D0+CDT*DXX
  CA0R = 1.0D0-CDT*DXX
  CB0 = 1.0D0+CDT*DYY
  CB0R = 1.0D0-CDT*DYY
  CC0 = 1.0D0+CDT*DZZ
  CC0R = 1.0D0-CDT*DZZ
!
  CT0X = -CDT*DXX/2.0D0
  CALA(NXX) = 0.0D0
  CGAA(NXX) = -1.0D0/CA0
  DO I = NXX, 1, -1
     CALA(I-1) = CT0X*CGAA(I)
     CGAA(I-1) = -1.0D0/(CA0+CT0X*CALA(I-1))
  END DO
!
  CT0Y = -CDT*DYY/2.0D0
  CALB(NYY) = 0.0D0
  CGAB(NYY) = -1.0D0/CB0
  DO J = NYY, 1, -1
     CALB(J-1) = CT0Y*CGAB(J)
     CGAB(J-1) = -1.0D0/(CB0+CT0Y*CALB(J-1))
  END DO
!
  CT0Z = -CDT*DZZ/2.0D0
  CALC(NZZ) = 0.0D0
  CGAC(NZZ) = -1.0D0/CC0
  DO K = NZZ, 1, -1
     CALC(K-1) = CT0Z*CGAC(K)
     CGAC(K-1) = -1.0D0/(CC0+CT0Z*CALC(K-1))
  END DO
END SUBROUTINE COEF

SUBROUTINE NU(CP, DT) ! Exact solution
!   Solves the partial differential equation with the potential and the
!   nonlinear term.
  USE GPE_DATA, ONLY : V, G, CI
  USE TEMP_DATA, ONLY : TMP3, P
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(IN) :: DT
  P = ABS(CP)
  TMP3 = P*P
  TMP3 = DT*(V + G*TMP3)
  CP = CP*EXP(-CI*TMP3)
END SUBROUTINE NU

SUBROUTINE LUX(CP)
!  Solves the partial differential equation only with the X space
!  derivative term using the Crank-Nicholson method
  USE COMM_DATA, ONLY : NX, NY, NZ, NXX
  USE CN_DATA, ONLY : CA0, CA0R, CT0X, CALA, CGAA
  USE TEMP_DATA, ONLY : CBE
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: CP
  INTEGER :: I
  COMPLEX (8), DIMENSION(0:NY,0:NZ) :: CXX
  CBE(NXX,:,:) = CP(NX,:,:)
  DO I = NXX, 1, -1
     CXX = -CT0X*CP(I+1,:,:)+CA0R*CP(I,:,:)-CT0X*CP(I-1,:,:)
     CBE(I-1,:,:) = CGAA(I)*(CT0X*CBE(I,:,:)-CXX)
  END DO
  CP(0,:,:) = 0.0D0
  DO I = 0, NXX
     CP(I+1,:,:) = CALA(I)*CP(I,:,:) + CBE(I,:,:)
  END DO
  CP(NX,:,:) = 0.0D0
END SUBROUTINE LUX

SUBROUTINE LUY(CP)
!  Solves the partial differential equation only with the Y space
!  derivative term using the Crank-Nicholson method
  USE COMM_DATA, ONLY : NX, NY, NZ, NYY
  USE CN_DATA, ONLY : CB0, CB0R, CT0Y, CALB, CGAB
  USE TEMP_DATA, ONLY : CBE
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: CP
  INTEGER :: J
  COMPLEX (8), DIMENSION(0:NX,0:NZ) :: CYY
  CBE(:,NYY,:) = CP(:,NY,:)
  DO J = NYY, 1, -1
     CYY = -CT0Y*CP(:,J+1,:)+CB0R*CP(:,J,:)-CT0Y*CP(:,J-1,:)
     CBE(:,J-1,:) = CGAB(J)*(CT0Y*CBE(:,J,:)-CYY)
  END DO
  CP(:,0,:) = 0.0D0
  DO J = 0, NYY
     CP(:,J+1,:) = CALB(J)*CP(:,J,:) + CBE(:,J,:)
  END DO
  CP(:,NY,:) = 0.0D0
END SUBROUTINE LUY

SUBROUTINE LUZ(CP)
!  Solves the partial differential equation only with the Z space
!  derivative term using the Crank-Nicholson method
  USE COMM_DATA, ONLY : NX, NY, NZ, NZZ
  USE CN_DATA, ONLY : CC0, CC0R, CT0Z, CALC, CGAC
  USE TEMP_DATA, ONLY : CBE
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: CP
  INTEGER :: K
  COMPLEX (8), DIMENSION(0:NX,0:NY) :: CZZ
  CBE(:,:,NZZ) = CP(:,:,NZ)
  DO K = NZZ, 1, -1
     CZZ = -CT0Z*CP(:,:,K+1)+CC0R*CP(:,:,K)-CT0Z*CP(:,:,K-1)
     CBE(:,:,K-1) = CGAC(K)*(CT0Z*CBE(:,:,K)-CZZ)
  END DO
  CP(:,:,0) = 0.0D0
  DO K = 0, NZZ
     CP(:,:,K+1) = CALC(K)*CP(:,:,K) + CBE(:,:,K)
  END DO
  CP(:,:,NZ) = 0.0D0
END SUBROUTINE LUZ

SUBROUTINE NORM(CP, ZNORM)
!  Calculates the normalization of the wave function and sets it to
!  unity.
  USE COMM_DATA, ONLY : NX, NY, NZ
  USE GPE_DATA, ONLY : DX, DY, DZ
  USE TEMP_DATA
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(OUT) :: ZNORM
!-----------------------------------------------------------  
  INTERFACE
    PURE FUNCTION SIMP(F, DX) RESULT (VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
!-----------------------------------------------------------
  INTEGER :: I, J
  P = ABS(CP)
  TMP3 = P*P
  FORALL (I = 0:NX)
     FORALL (J = 0:NY) TMP2(I,J) = SIMP(TMP3(I,J,:), DZ)
     TMP1(I) = SIMP(TMP2(I,:), DY)
  END FORALL
  ZNORM = SQRT(SIMP(TMP1, DX))
  CP = CP/ZNORM
END SUBROUTINE NORM

SUBROUTINE RAD(CP, RMS)
  USE COMM_DATA, ONLY : NX, NY, NZ
  USE GPE_DATA, ONLY : DX, DY, DZ, R2
  USE TEMP_DATA
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:,0:,0:), INTENT(IN) :: CP
  REAL (8), INTENT(OUT) :: RMS   
!-----------------------------------------------------------
  INTERFACE
    PURE FUNCTION SIMP(F, DX) RESULT (VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
!-----------------------------------------------------------
  INTEGER :: I, J
  P = ABS(CP)
  TMP3 = R2*P*P
  FORALL (I = 0:NX)
     FORALL (J = 0:NY) TMP2(I,J) = SIMP(TMP3(I,J,:), DZ)
     TMP1(I) = SIMP(TMP2(I,:), DY)
  END FORALL
  RMS = SQRT(SIMP(TMP1, DX))
END SUBROUTINE RAD

SUBROUTINE CHEM(CP, MU, EN)
  USE COMM_DATA, ONLY : NX, NY, NZ, NXX, NYY, NZZ
  USE GPE_DATA, ONLY : DX, DY, DZ, X2, V, G
  USE TEMP_DATA
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:,0:,0:), INTENT(IN) :: CP
  REAL (8), INTENT(OUT) :: MU, EN
!-----------------------------------------------------------
  INTERFACE
    PURE FUNCTION DIFF(P, DX) RESULT (DP)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: P
      REAL (8), INTENT(IN) :: DX
      REAL (8), DIMENSION(0:SIZE(P)-1) :: DP
    END FUNCTION DIFF
  END INTERFACE
!-----------------------------------------------------------
  INTERFACE
    PURE FUNCTION SIMP(F, DX) RESULT (VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
!-----------------------------------------------------------
  INTEGER :: I, J, K
  REAL (8), DIMENSION(0:NX,0:NY,0:NZ) :: DPX, DPY, DPZ
  REAL (8), DIMENSION(0:NX,0:NY,0:NZ) :: P2, GP2, DP2
!
  P = ABS(CP)
  FORALL (I=0:NX, J=0:NY) DPZ(I,J,:) = DIFF(P(I,J,:), DZ)
  FORALL (I=0:NX, K=0:NZ) DPY(I,:,K) = DIFF(P(I,:,K), DY)
  FORALL (J=0:NY, K=0:NZ) DPX(:,J,K) = DIFF(P(:,J,K), DX)
!
  P2 = P*P
  GP2 = G*P2
  DP2 = DPX*DPX + DPY*DPY + DPZ*DPZ
  TMP3 = (V + GP2)*P2 + DP2
  FORALL (I = 0:NX, J = 0:NY) TMP2(I,J) = SIMP(TMP3(I,J,:), DZ)
  FORALL (I = 0:NX) TMP1(I) = SIMP(TMP2(I,:), DY)
  MU = SIMP(TMP1, DX)
  TMP3 = (V + GP2/2.0D0)*P2 + DP2
  FORALL (I = 0:NX, J = 0:NY) TMP2(I,J) = SIMP(TMP3(I,J,:), DZ)
  FORALL (I = 0:NX) TMP1(I) = SIMP(TMP2(I,:), DY)
  EN = SIMP(TMP1, DX)
END SUBROUTINE CHEM

PURE FUNCTION SIMP(F, DX) RESULT (VALUE)
! Does the spatial integration with Simpson's rule.
! N refer to the number of integration points, DX space step, and
! F is the function to be integrated.
  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(IN) :: F
  REAL (8), INTENT(IN) :: DX
  REAL (8) :: VALUE
  REAL (8) :: F1, F2
  INTEGER :: I, N
  N = SIZE(F) - 1
  F1 = F(1) + F(N-1) ! N EVEN
  F2 = F(2) 
  DO I = 3, N-3, 2
     F1 = F1 + F(I)
     F2 = F2 + F(I+1)
  END DO
  VALUE = DX*(F(0) + 4.0D0*F1 + 2.0D0*F2 + F(N))/3.0D0
END FUNCTION SIMP

PURE FUNCTION DIFF(P,DX) RESULT (DP)
  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(IN) :: P
  REAL (8), INTENT(IN) :: DX
  REAL (8), DIMENSION(0:SIZE(P)-1) :: DP
  INTEGER :: I, N
  N = SIZE(P) - 1
  DP(0) = 0.0D0
  DP(1) = (P(2) - P(0))/(2.0D0*DX)
  FORALL(I=2:N-2)
    DP(I) = (P(I-2)-8.0D0*P(I-1)+8.0D0*P(I+1)-P(I+2))/(12.0D0*DX)
  END FORALL
  DP(N-1) = (P(N) - P(N-2))/(2.0D0*DX)
  DP(N) = 0.0D0
END FUNCTION DIFF
