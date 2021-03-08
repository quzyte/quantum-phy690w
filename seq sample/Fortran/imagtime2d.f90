!  Filename : imagtime2d.f90
!  Title : Fortran program for time-dependent Gross-Pitaevskii equation in 
!          a fully anisotropic trap
!  Authors : P. Muruganandam and S. K. Adhikari
! 
!  Fortran program for Gross-Pitaevskii equation in two-dimensional 
!  anisotropic trap by imaginary time propagation (Fortran 90/95 version)
!
MODULE COMM_DATA
! NX, NY : Number of space mesh points (X and Y)  
  INTEGER, PARAMETER :: NX = 800, NXX = NX-1, NX2 = NX/2 ! No. of X steps
  INTEGER, PARAMETER :: NY = 800, NYY = NY-1, NY2 = NY/2 ! No. of Y steps
! NPAS : Number of subsequent iterations with fixed nonlinearity, and 
! NRUN : Number of final time iterations with fixed nonlinearity
  INTEGER, PARAMETER :: NPAS = 30000, NRUN = 5000
  REAL (8), PARAMETER :: PI = 3.14159265358979D0
END MODULE COMM_DATA

MODULE GPE_DATA
  USE COMM_DATA, ONLY : NX, NY
! DX, DY : Space step and DT : Time step
  REAL (8), PARAMETER :: DX = 0.02D0, DY = 0.02D0, DT = 0.0001D0
! G0 : Final nonlinearity, AL : Anisotrophy coefficient
  REAL (8), PARAMETER :: G0 = 12.5484D0 ! Final nonlinearity
  REAL (8), PARAMETER :: AL = 2.0D0 ! Anisotropy 'LAMBDA' of the trap
  ! OPTION = 1 Solves -psi_xx-psi_yy+V(x,y)psi+G0|psi|^2 psi=i psi_t
  ! OPTION = 2 Solves [-psi_xx-psi_yy+V(x,y)psi]/2+G0|psi|^2 psi=i psi_t
  ! OPTION = 3 Solves -psi_xx-psi_yy+V(x,y)psi/4+G0|psi|^2 psi=i psi_t
  INTEGER, PARAMETER :: OPTION = 2
! X(0:N), Y(0:N) : Space mesh, V(0:N) : Potential, and 
! CP(0:N) : Wave function 
  REAL (8), DIMENSION(0:NX, 0:NY) :: V, R2, CP
  REAL (8), DIMENSION(0:NX) :: X, X2
  REAL (8), DIMENSION(0:NY) :: Y, Y2
  REAL (8) :: G, XOP
END MODULE GPE_DATA

MODULE CN_DATA
  USE COMM_DATA, ONLY : NX, NY
  REAL (8), DIMENSION(0:NX) :: CALA, CGAA
  REAL (8), DIMENSION(0:NY) :: CALB, CGAB
  REAL (8) :: CT0X, CT0Y
  REAL (8) :: CA0, CB0, CA0R, CB0R
END MODULE CN_DATA 

 PROGRAM GROSS_PITAEVSKII_SSCN_2D
   USE COMM_DATA, ONLY : NX, NY, NX2, NY2, NPAS, NRUN
   USE GPE_DATA, ONLY : DT, DX, DY, OPTION, XOP, AL, G, G0, X, Y, CP
   IMPLICIT NONE
! Subroutine INTIALIZE() used to initialize the space mesh X(I), 
! potential V(I) and the initial wave function. Subroutine COEF() used to
! generate the coefficients for the Crank-Nicholson Scheme. The routine
! NU() performs time progation for the non-derivative part and LU() performs
! time propagation of derivative part. NORM() calculates the norm and 
! normalizes the wave function, CHEM() and RAD() are used to calculate the 
! chemical potential, energy and the rms radius, respectively. The functions 
! DIFF(), used to calculate the space derivatives of the wave function used 
! in CHEM(). The function SIMP() does the integration by Simpson's rule.
!------------------------ INTERFACE BLOCKS -----------------------
   INTERFACE 
     SUBROUTINE INITIALIZE(CP)
       IMPLICIT NONE
       REAL (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
     END SUBROUTINE INITIALIZE
   END INTERFACE
!
   INTERFACE 
     SUBROUTINE COEF()
       IMPLICIT NONE
     END SUBROUTINE COEF
   END INTERFACE
!
   INTERFACE 
     SUBROUTINE NU(CP, DT)
       IMPLICIT NONE
       REAL (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
       REAL (8), INTENT(IN) :: DT
     END SUBROUTINE NU
   END INTERFACE
!
   INTERFACE 
     SUBROUTINE LUX(CP)
       IMPLICIT NONE
       REAL (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
     END SUBROUTINE LUX
   END INTERFACE
!
   INTERFACE 
     SUBROUTINE LUY(CP)
       IMPLICIT NONE
       REAL (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
     END SUBROUTINE LUY
   END INTERFACE
!
   INTERFACE
     SUBROUTINE NORM(CP, ZNORM)
       IMPLICIT NONE
       REAL (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
       REAL (8), INTENT(OUT) :: ZNORM
     END  SUBROUTINE NORM
   END INTERFACE
!
   INTERFACE
     SUBROUTINE RAD(CP, RMS)
       IMPLICIT NONE
       REAL (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
       REAL (8), INTENT(OUT) :: RMS
     END  SUBROUTINE RAD
   END INTERFACE
!
   INTERFACE
     SUBROUTINE CHEM(CP, MU, EN)
       IMPLICIT NONE
       REAL (8), DIMENSION(0:,0:), INTENT(IN) :: CP
       REAL (8), INTENT(OUT) :: MU, EN
     END  SUBROUTINE CHEM
   END INTERFACE
!------------------------ END INTERFACE BLOCKS -------------------
  INTEGER :: I, J, K
  REAL (8) :: ZNORM, RMS, MU, EN

  SELECT CASE (OPTION)
    CASE (1,3)
      XOP = 1.0D0
    CASE (2)
      XOP = 2.0D0
    CASE (:0,4:)
      PRINT *, 'ERROR: Wrong option', OPTION
      STOP
  END SELECT
  G = 0.0D0
  WRITE(7,900) OPTION
  WRITE(7,901) AL
  WRITE(7,*)  
  WRITE(7,902) NX, NY
  WRITE(7,903) NPAS, NRUN  
  WRITE(7,904) G0 
  WRITE(7,905) DX, DY 
  WRITE(7,906) DT 
  WRITE(7,*)
  900 FORMAT('  OPTION = ', I3)
  901 FORMAT('  Anisotropy AL = ', F12.6)
  902 FORMAT('# Space Stp NX = ', I8, ', NY = ', I8)
  903 FORMAT('# Time Stp : NPAS = ',I9,', NRUN = ',I9)
  904 FORMAT('  Nonlinearity G = ', F16.8)
  905 FORMAT('  Space Step DX = ', F10.6, ', DY = ', F10.6)
  906 FORMAT('  Time Step  DT = ', F10.6)
! INITIALIZE() initializes the starting normalized wave function 'CP' and
! harmonic oscillator potential V
  CALL INITIALIZE(CP)
! COEF() defines the coefficients of Crank-Nicholson Scheme.
  CALL COEF()
! NORM() calculates norm and restores normalization
  CALL NORM(CP, ZNORM)
! RAD() calculates the r.m.s radius RMS
  CALL RAD(CP, RMS)
! CHEM() calculates the chemical potential MU and energy EN.
  CALL CHEM(CP, MU, EN)
!
  DO I = 0, NX, 10  ! Writes initial wave function File 1
     DO J = 0, NY, 10
        WRITE(1, 999) X(I), Y(J), CP(I,J)
     END DO
     WRITE(1,*)
  END DO
  999 FORMAT (2F12.6, F16.8)
  WRITE (7, 1001)
  WRITE (7, 1002)
  WRITE (7, 1001)
  WRITE (7, 1003) ZNORM, MU/XOP, EN/XOP, RMS, CP(NX2, NY2)
  1001 FORMAT (20X,'---------------------------------------------------')
  1002 FORMAT (20X, 'Norm', 7X, 'Chem', 8X, 'Ener', 7X, '<r>', 7X, &
          'Psi(0)')
  1003 FORMAT ('Initial : ', 4X, F11.4, 2F12.6, 3F11.5)
  G = XOP*G0
  DO K = 1, NPAS ! 'NPAS' Iterations transient
! NU() performs time propagation with harmonic potential and nonlinear
! term (non-derivative parts)
     CALL NU(CP, DT)
! LUX() and LUY() perform the time iteration with space derivative in
! X and Y, respectively, using Crank-Nicholson scheme.
     CALL LUX(CP)
     CALL NORM(CP,ZNORM)
     CALL LUY(CP)
     CALL NORM(CP,ZNORM)
  END DO
  CALL CHEM(CP, MU, EN)
  CALL RAD(CP, RMS)
  WRITE (7, 1005) ZNORM, MU/XOP, EN/XOP, RMS, CP(NX2, NY2)
  1005 FORMAT('After NPAS iter.:', F8.4, 2F12.6, 3F11.5)
  DO I = 0, NX, 10  ! Writes intermediate wave function File 2
     DO J = 0, NY, 10
        WRITE(2, 999) X(I), Y(J), CP(I,J)
     END DO
     WRITE(2,*)
  END DO
!
  DO K = 1, NRUN 
     CALL NU(CP, DT)
     CALL LUX(CP)
     CALL NORM(CP,ZNORM)
     CALL LUY(CP)
     CALL NORM(CP,ZNORM)
  END DO
  CALL CHEM(CP, MU, EN)
  CALL RAD(CP, RMS)
  WRITE (7, 1006) ZNORM, MU/XOP, EN/XOP, RMS, CP(NX2, NY2)
  WRITE (7, 1001)
  1006 FORMAT('After NRUN iter.:',F8.4, 2F12.6, 3F11.5)
  DO I = 0, NX, 10 ! Writes final wave funtion in File 3
     DO J = 0, NY, 10
        WRITE(3,999) X(I), Y(J), CP(I,J)
     END DO
     WRITE(3,*)
  END DO
END PROGRAM GROSS_PITAEVSKII_SSCN_2D

SUBROUTINE INITIALIZE(CP)
  ! Routine that initizlizes the constant and variables.
  ! Calculates the potential term V and the initial wave function CP
  USE COMM_DATA, ONLY : NX, NY, NX2, NY2, PI
  USE GPE_DATA, ONLY : AL, OPTION, DX, DY, X, X2, Y, Y2, V, R2
  IMPLICIT NONE
  REAL (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP

  REAL (8), DIMENSION(0:NX,0:NY):: TMP2D
  REAL (8) :: SPI, SP, AL2, QR_AL
  INTEGER :: I

  SPI = SQRT(2.0D0*PI)
  SP = SQRT(PI)
  QR_AL = SQRT(SQRT(AL))
  AL2 = AL*AL
!
  FORALL (I=0:NX) X(I) = (I-NX2)*DX
  FORALL (I=0:NY) Y(I) = (I-NY2)*DY
  X2 = X*X
  Y2 = Y*Y
  FORALL (I=0:NX) R2(I,:) = X2(I) + Y2
  SELECT CASE (OPTION)
    CASE (1,2)
       FORALL (I=0:NX) 
          V(I,:) = (X2(I) + AL2*Y2)
          TMP2D(I,:) = (X2(I) + AL*Y2)/2.0D0
       END FORALL
       CP = QR_AL*EXP(-TMP2D)/SP
    CASE (3)
       FORALL (I=0:NX) 
          V(I,:) = (X2(I) + AL2*Y2)/4.0D0
          TMP2D(I,:) = (X2(I) + AL*Y2)/4.0D0
       END FORALL
       CP = QR_AL*EXP(-TMP2D)/SPI
  END SELECT
END SUBROUTINE INITIALIZE
!
SUBROUTINE COEF()
  ! Calculates the coefficients needed in subroutine LUX and LUY.
  USE COMM_DATA, ONLY : NX, NY, NXX, NYY
  USE GPE_DATA, ONLY : DX, DY, DT
  USE CN_DATA
  IMPLICIT NONE
  INTEGER :: I, J
  REAL (8) :: DX2, DY2
  REAL (8) :: DXX, DYY
!
  DX2 = DX*DX   ! Generating the Coefficients for the
  DY2 = DY*DY   ! C-N method (to solve the spatial part)           
  DXX = 1.0D0/DX2
  DYY = 1.0D0/DY2
  CA0 = 1.0D0+DT*DXX
  CA0R = 1.0D0-DT*DXX
  CB0 = 1.0D0+DT*DYY
  CB0R = 1.0D0-DT*DYY
!
  CT0X = -DT*DXX/2.0D0
  CALA(NXX) = 0.0D0
  CGAA(NXX) = -1.0D0/CA0
  DO I = NXX, 1, -1
     CALA(I-1) = CT0X*CGAA(I)
     CGAA(I-1) = -1.0D0/(CA0+CT0X*CALA(I-1))
  END DO
!
  CT0Y = -DT*DYY/2.0D0
  CALB(NYY) = 0.0D0
  CGAB(NYY) = -1.0D0/CB0
  DO J = NYY, 1, -1
     CALB(J-1) = CT0Y*CGAB(J)
     CGAB(J-1) = -1.0D0/(CB0+CT0Y*CALB(J-1))
  END DO
END SUBROUTINE COEF

SUBROUTINE NU(CP, DT) ! Exact solution
  ! Solves the partial differential equation with the potential and 
  ! the nonlinear term.
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : V, G
  IMPLICIT NONE
  REAL (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(IN) :: DT
  REAL (8), DIMENSION(0:NX,0:NY) :: P2, TMP2D
  P2 = CP*CP
  TMP2D = V + G*P2
  CP = CP*EXP(-DT*TMP2D)
END SUBROUTINE NU

SUBROUTINE LUX(CP)
  ! Solves the partial differential equation only with the X-space
  ! derivative term using the Crank-Nicholson method
  USE COMM_DATA, ONLY : NX, NY, NXX
  USE CN_DATA, ONLY : CA0, CA0R, CT0X, CALA, CGAA
  IMPLICIT NONE
  REAL (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
  INTEGER :: I
  REAL (8) :: CBE(0:NXX,0:NY), CXX(0:NY)
  CBE(NXX,:) = CP(NXX,:)
  DO I = NXX, 1, -1
     CXX = -CT0X*CP(I+1,:)+CA0R*CP(I,:)-CT0X*CP(I-1,:)
     CBE(I-1,:) = CGAA(I)*(CT0X*CBE(I,:)-CXX)
  END DO
  CP(0,:) = 0.0D0
  DO I = 0, NXX
     CP(I+1,:) = CALA(I)*CP(I,:) + CBE(I,:)
  END DO
  CP(NX,:) = 0.0D0
END SUBROUTINE LUX

SUBROUTINE LUY(CP)
  ! Solves the partial differential equation only with the Y-space
  ! derivative term using the Crank-Nicholson method
  USE COMM_DATA, ONLY : NX, NY, NYY
  USE CN_DATA, ONLY : CB0, CB0R, CT0Y, CALB, CGAB
  IMPLICIT NONE
  REAL (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
  INTEGER :: J
  REAL (8) :: CBE(0:NX,0:NYY), CYY(0:NX)
  CBE(:,NYY) = CP(:,NYY)
  DO J = NYY, 1, -1
     CYY = -CT0Y*CP(:,J+1)+CB0R*CP(:,J)-CT0Y*CP(:,J-1)
     CBE(:,J-1) = CGAB(J)*(CT0Y*CBE(:,J)-CYY)
  END DO
  CP(:,0) = 0.0D0
  DO J = 0, NYY
     CP(:,J+1) = CALB(J)*CP(:,J) + CBE(:,J)
  END DO
  CP(:,0) = CP(:,1)
  CP(:,NY) = 0.0D0
END SUBROUTINE LUY

SUBROUTINE NORM(CP, ZNORM)
  ! Calculates the normalization of the wave function and sets it to
  ! unity.
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : DX, DY
  IMPLICIT NONE
  REAL (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(OUT) :: ZNORM
   
  INTERFACE
    PURE FUNCTION SIMP(F, DX)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: SIMP
    END FUNCTION SIMP
  END INTERFACE

  INTEGER :: I
  REAL (8), DIMENSION(0:NX,0:NY) :: TMP2D
  REAL (8), DIMENSION(0:NX) :: TMP1D
  TMP2D = CP*CP
  FORALL (I = 0:NX) TMP1D(I) = SIMP(TMP2D(I,:), DY)
  ZNORM = SQRT(SIMP(TMP1D, DX))
  CP = CP/ZNORM
END SUBROUTINE NORM

SUBROUTINE RAD(CP, RMS) ! Calculates the root mean square size RMS
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : DX, DY, R2
  IMPLICIT NONE
  REAL (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(OUT) :: RMS
   
  INTERFACE
    PURE FUNCTION SIMP(F, DX)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: SIMP
    END FUNCTION SIMP
  END INTERFACE

  INTEGER :: I
  REAL (8), DIMENSION(0:NX,0:NY) :: TMP2D
  REAL (8), DIMENSION(0:NX) :: TMP1D
  TMP2D = R2*CP*CP
  FORALL (I = 0:NX) TMP1D(I) = SIMP(TMP2D(I,:), DY)
  RMS = SQRT(SIMP(TMP1D, DX))
END SUBROUTINE RAD

SUBROUTINE CHEM(CP, MU, EN)
  ! Calculates the chemical potential MU and energy EN.  CP is the wave
  ! function, V is the potential and G is the nonlinearity.
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : DX, DY, V, G
  IMPLICIT NONE
  REAL (8), DIMENSION(0:,0:), INTENT(IN) :: CP
  REAL (8), INTENT(OUT) :: MU, EN
  !----------------------------------------------------------------------
  INTERFACE
    PURE FUNCTION DIFF(P, DX) RESULT (DP)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: P
      REAL (8), INTENT(IN) :: DX
      REAL (8), DIMENSION(0:SIZE(P)-1) :: DP
    END FUNCTION DIFF
  END INTERFACE
  !----------------------------------------------------------------------
  INTERFACE
    PURE FUNCTION SIMP(F, DX)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: SIMP
    END FUNCTION SIMP
  END INTERFACE
  !----------------------------------------------------------------------
  INTEGER :: I
  REAL (8), DIMENSION(0:NX,0:NY) :: DPX, DPY, P2, GP2, DP2, TMP2D, EMP2D
  REAL (8), DIMENSION(0:NX) :: TMP1D, EMP1D

  FORALL (I=0:NX) DPY(I,:) = DIFF(CP(I,:), DY)
  FORALL (I=0:NY) DPX(:,I) = DIFF(CP(:,I), DX)

  P2 = CP*CP
  DP2 = DPX*DPX + DPY*DPY
  GP2 = G*P2
  TMP2D = (V + GP2)*P2 + DP2
  EMP2D = (V + GP2/2.0D0)*P2 + DP2
  FORALL (I = 0:NX) 
     TMP1D(I) = SIMP(TMP2D(I,:), DY)
     EMP1D(I) = SIMP(EMP2D(I,:), DY)
  END FORALL
  MU = SIMP(TMP1D, DX)
  EN = SIMP(EMP1D, DX)
END SUBROUTINE CHEM

PURE FUNCTION SIMP(F, DX)
! Does the spatial integration with Simpson's rule.
! N refer to the number of integration points, DX space step, and
! F is the function to be integrated.
  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(IN) :: F
  REAL (8), INTENT(IN) :: DX
  REAL (8) :: SIMP

  REAL (8) :: F1, F2
  INTEGER :: I, N

  N = SIZE(F) - 1

  F1 = F(1) + F(N-1) ! N EVEN
  F2 = F(2) 
  DO I = 3, N-3, 2
     F1 = F1 + F(I)
     F2 = F2 + F(I+1)
  END DO
  SIMP = DX*(F(0) + 4.0D0*F1 + 2.0D0*F2 + F(N))/3.0D0
END FUNCTION SIMP

PURE FUNCTION DIFF(P,DX) RESULT (DP)
  ! Computes the first derivative DP of P using
  ! Richardsonextrapolation formula. The derivative at the  
  ! boundaries are assumed to be zero
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
