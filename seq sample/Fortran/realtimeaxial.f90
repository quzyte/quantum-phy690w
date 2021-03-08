!  Filename : imagtimeaxial.f90
!!  Title : Fortran program for time-dependent Gross-Pitaevskii equation in 
!          a fully anisotropic trap
!  Authors : P. Muruganandam and S. K. Adhikari
! 
!  Fortran program for Gross-Pitaevskii equation in axially symmetric 
!  trap by imaginary time propagation (Fortran 90/95 Version)
!
!
MODULE COMM_DATA
! NX, NY : Number of space mesh points (X and Y)  
  INTEGER, PARAMETER :: NX = 130, NXX = NX-1 ! No. of X steps
  INTEGER, PARAMETER :: NY = 130, NYY = NY-1, NY2 = NY/2 ! No. of Y steps
! NSTP : Number of iterations to introduce the nonlinearity,
! NPAS : Number of subsequent iterations with fixed nonlinearity, and 
! NRUN : Number of final time iterations with fixed nonlinearity
  INTEGER, PARAMETER :: NSTP = 100000, NPAS = 1000, NRUN = 5000
END MODULE COMM_DATA

MODULE GPE_DATA
  USE COMM_DATA, ONLY : NX, NY
  REAL (8), PARAMETER :: PI = 3.14159265358979D0
! DX, DY : Space step and DT : Time step
  REAL (8), PARAMETER :: DX = 0.1D0, DY = 0.1D0, DT = 0.001D0
! G0 : Final nonlinearity, KAP, LAM : Anisotrophy coefficients
  REAL (8), PARAMETER :: G0 = 18.81D0 ! Final nonlinearity
  REAL (8), PARAMETER :: LAM = 4.0D0, KAP = 1.0D0 ! Lambda & Kappa
  COMPLEX (8), PARAMETER :: CI = (0.0D0, 1.0D0)
! OPTION and XOP decides which equation to be solved. 
  ! OPTION = 1 Solves -psi_xx-psi_yy+V(x,y)psi+G0|psi|^2 psi=i psi_t
  ! OPTION = 2 Solves [-psi_xx-psi_yy+V(x,y)psi]/2+G0|psi|^2 psi=i psi_t
  ! OPTION = 3 Solves -psi_xx-psi_yy+V(x,y)psi/4+G0|psi|^2 psi=i psi_t
  INTEGER, PARAMETER :: OPTION = 2
  REAL (8) :: G, XOP, TPI
! X(0:N), Y(0:N) : Space mesh, V(0:N) : Potential, and 
! CP(0:N) : Wave function (COMPLEX)
  REAL (8), DIMENSION(0:NX, 0:NY) :: V
  REAL (8), DIMENSION(0:NX) :: X, X2
  REAL (8), DIMENSION(0:NY) :: Y, Y2
  COMPLEX (8), DIMENSION(0:NX, 0:NY) :: CP
END MODULE GPE_DATA

MODULE CN_DATA
  USE COMM_DATA, ONLY : NX, NY
  COMPLEX (8), DIMENSION(0:NX) :: CALX, CGAX, CAIPX
  COMPLEX (8), DIMENSION(1:NX) :: CPX
  COMPLEX (8), DIMENSION(0:NY) :: CALY, CGAY
  COMPLEX (8) :: CAIPMX, CAIPMY
  !COMPLEX (8) :: CAI0, CAIM
END MODULE CN_DATA 

PROGRAM GROSS_PITAEVSKII_SSCN_AXIAL
  USE COMM_DATA
  USE GPE_DATA
  IMPLICIT NONE
! Subroutine INTIALIZE() used to initialize the space mesh X(I), 
! potential V(I) and the initial wave function. Subroutine COEF() used to
! generate the coefficients for the Crank-Nicholson Scheme. The routine
! NU() performs time progation for the non-derivative part and LU() performs
! time propagation of derivative part. NORM() calculates the norm and 
! normalizes the wave function, CHEM() and RAD() are used to calculate the 
! chemical potential, energy and the rms radius, respectively. The function 
! DIFF() used to calculate the space derivatives of the wave function used 
! in CHEM() and SIMP() does the integration by Simpson's rule.
!------------------------ INTERFACE BLOCKS -----------------------
  INTERFACE 
    SUBROUTINE INITIALIZE()
      IMPLICIT NONE
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
      COMPLEX (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
      REAL (8), INTENT(IN) :: DT
    END SUBROUTINE NU
  END INTERFACE
!
  INTERFACE 
    SUBROUTINE LUX(CP)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
    END SUBROUTINE LUX
  END INTERFACE
!
  INTERFACE 
    SUBROUTINE LUY(CP)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
    END SUBROUTINE LUY
  END INTERFACE
!
  INTERFACE
    SUBROUTINE NORM(CP, ZNORM)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
      REAL (8), INTENT(OUT) :: ZNORM
    END  SUBROUTINE NORM
  END INTERFACE
!
  INTERFACE
    SUBROUTINE RAD(CP, RHORMS, ZRMS)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
      REAL (8), INTENT(OUT) :: RHORMS, ZRMS
    END  SUBROUTINE RAD
  END INTERFACE
!
  INTERFACE
    SUBROUTINE CHEM(CP, MU, EN)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:,0:), INTENT(IN) :: CP
      REAL (8), INTENT(OUT) :: MU, EN
    END  SUBROUTINE CHEM
  END INTERFACE
!------------------------ END INTERFACE BLOCKS -------------------
  INTEGER :: I, J, K
  REAL (8) :: GSTP, ZNORM, RHORMS, ZRMS, MU, EN
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
  WRITE(7,901) KAP, LAM
  WRITE(7,*)  
  WRITE(7,902) NX, NY
  WRITE(7,903) NSTP, NPAS, NRUN  
  WRITE(7,904) G0 
  WRITE(7,905) DX, DY 
  WRITE(7,906) DT 
  WRITE(7,*)
  900 FORMAT('  OPTION = ', I3)
  901 FORMAT('  Anisotropy KAP = ',F12.6,', LAM = ',F12.6)
  902 FORMAT('# Space Stp NX = ', I8, ', NY = ', I8)
  903 FORMAT('# Time Stp : NSTP = ', I9, ', NPAS = ',I9,', NRUN = ',I9)
  904 FORMAT('  Nonlinearity G = ', F16.8)
  905 FORMAT('  Space Step DX = ', F10.6, ', DY = ', F10.6)
  906 FORMAT('  Time Step  DT = ', F10.6)
  CALL INITIALIZE()
  CALL COEF()
  CALL NORM(CP, ZNORM)
  CALL CHEM(CP, MU, EN)
  CALL RAD(CP, RHORMS, ZRMS)
  WRITE (7, 1001)
  WRITE (7, 1002)
  WRITE (7, 1001)
  WRITE (7, 1003) ZNORM, MU/XOP, EN/XOP, RHORMS, ZRMS, ABS(CP(0,NY2+1))
  1001 FORMAT (12X,55('-'))
  1002 FORMAT (14X,'Norm',4X,'Chem',5X,'Energy',4X,'<rho>',5X,'<z>', &
       7X,'psi(0)')
 1003 FORMAT ('Initial :  ',F6.3,5(F10.3))
  DO I = 0, NX, 4 ! Writes initial wave function
     DO J = 0, NY, 4
        WRITE (1,999) X(I), Y(J), ABS(CP(I,J))
     END DO
     WRITE (1,*)
  END DO
  999 FORMAT (2F12.3, F16.8)
! 
  GSTP = XOP*G0/NSTP
  DO K = 1, NSTP ! Introduces nonlinearity in 'NSTP' iterations
     G = G + GSTP
! NU() performs time propagation with harmonic potential and nonlinear
! term (non-derivative parts)
     CALL NU(CP, DT)
! LUX() and LUY() perform the time iteration with space derivative in
! X and Y, respectively, using Crank-Nicholson scheme.
     CALL LUX(CP)
     CALL LUY(CP)
  END DO
  CALL NORM(CP, ZNORM)
  CALL CHEM(CP, MU, EN)
  CALL RAD(CP, RHORMS, ZRMS)
  WRITE (7, 1004) ZNORM, MU/XOP, EN/XOP, RHORMS, ZRMS, ABS(CP(0,NY2+1))
 1004 FORMAT ('NSTP iter :',F6.3,  5(F10.3))
  DO I = 0, NX, 4 ! Writes intermediate wave funtion in File 2
     DO J = 0, NY, 4
        WRITE (2,999) X(I), Y(J), ABS(CP(I,J))
     END DO
     WRITE (2,*)
  END DO
!
  DO K = 1, NPAS ! 'NPAS' Iterations transient
     CALL NU(CP, DT)
     CALL LUX(CP)
     CALL LUY(CP)
  END DO
  CALL NORM(CP, ZNORM)
  CALL CHEM(CP, MU, EN)
  CALL RAD(CP, RHORMS, ZRMS)
  WRITE (7, 1005) ZNORM, MU/XOP, EN/XOP, RHORMS, ZRMS, ABS(CP(0,NY2+1))
  WRITE (7, 1001)
 1005 FORMAT ('NPAS iter :',F6.3,  5(F10.3))
  DO I = 0, NX, 4 ! Writes final wavefuntion in File 3
     DO J = 0, NY, 4
        WRITE (3,999) X(I), Y(J), ABS(CP(I,J))
     END DO
     WRITE (3,*)
  END DO
!
!  The following line defines a nonstationary problem which is studied
!  below and the time evolution written on file 8 via WRITE(8,*)
!
  G = 0.5D0*G         
  DO K = 1, NRUN 
     CALL NU(CP, DT)
     CALL LUX(CP)
     CALL LUY(CP)
     IF (MOD(K,20).EQ.0) THEN          
        CALL RAD(CP, RHORMS, ZRMS)
        WRITE(8,1999) K*DT*XOP, RHORMS, ZRMS
     END IF       
  END DO
  1999 FORMAT(F12.6, 2F16.8)
END PROGRAM GROSS_PITAEVSKII_SSCN_AXIAL

SUBROUTINE INITIALIZE()
!   Routine that initizlizes the constant and variables.
!   Calculates the potential term V and the initial wave function CP
  USE COMM_DATA, ONLY : NX, NY, NY2
  USE GPE_DATA
  IMPLICIT NONE
  REAL (8), DIMENSION(0:NX,0:NY):: TMP2D
  REAL (8) :: PI3, FAC, FAC3, KAP2, LAM2
  INTEGER :: I
  KAP2 = KAP * KAP
  LAM2 = LAM * LAM
  PI3 = PI**3
  TPI = 2.0D0*PI
  FAC = SQRT(SQRT(LAM*KAP2/(PI3)))
  FAC3 = FAC/SQRT(SQRT(8.0D0))
  FORALL (I=0:NX) X(I) = I*DX
  FORALL (I=0:NY) Y(I) = (I-NY2)*DY
  X2 = X*X
  Y2 = Y*Y
  SELECT CASE (OPTION)
    CASE (1,2)
       FORALL(I=0:NX) V(I,:) = (KAP2*X2(I) + LAM2*Y2)
       FORALL(I=0:NX) TMP2D(I,:) = (KAP*X2(I) + LAM*Y2)/2.0D0
       CP = FAC*EXP(-TMP2D)
    CASE (3)
       FORALL(I=0:NX) V(I,:) = (KAP2*X2(I) + LAM2*Y2)/4.0D0
       FORALL(I=0:NX) TMP2D(I,:) = (KAP*X2(I) + LAM*Y2)/4.0D0
       CP = FAC3*EXP(-TMP2D)
  END SELECT
END SUBROUTINE INITIALIZE
!
SUBROUTINE COEF()
!   Routine that initizlizes the constant and variables.
!   Calculates the potential term V and the initial wave function CP
  USE COMM_DATA, ONLY : NX, NY, NXX, NYY
  USE GPE_DATA, ONLY : DX, DY, DT, CI, X
  USE CN_DATA
  IMPLICIT NONE
  INTEGER :: J
  REAL (8) :: DX2, DY2
  COMPLEX (8) :: CDT, CAI0, CAIM

  CDT = DT/CI
  DX2 = DX*DX   ! Generating the Coefficients for the
  DY2 = DY*DY   ! C-N method (to solve the spatial part)           
  CAIPMX = CDT/(2.0D0*DX2)
  CAI0 = 1.D0 - CDT/DX2
  CALX(1) = 1.0D0
  CGAX(1) = -1.0D0/(CAI0+CAIPMX-CDT/(4.D0*DX*X(1)))
!
  DO J = 1, NXX
     CPX(J) = CDT/(4.0D0*DX*X(J))
     CAIPX(J) = CAIPMX-CPX(J)
     CAIM = CAIPMX+CPX(J)
     CALX(J+1)=CGAX(J)*CAIM
     CGAX(J+1)=-1.0D0/(CAI0+(CAIPMX-CDT/(4.D0*DX*X(J+1)))*CALX(J+1))
  END DO
!
  CAIPMY = CDT/(2.0D0*DY2)
  CAI0 = 1.D0 - CDT/DY2
  CALY(NYY) = 0.0D0
  CGAY(NYY) = -1.0D0/CAI0
!
  DO J = NYY, 1, -1
     CALY(J-1) = CGAY(J)*CAIPMY
     CGAY(J-1) = -1.0D0/(CAI0+CAIPMY*CALY(J-1))
  END DO
END SUBROUTINE COEF

SUBROUTINE NU(CP, DT) ! Exact solution
! Solves the partial differential equation with the potential and the
! nonlinear term.
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : V, G, CI
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(IN) :: DT
  REAL (8), DIMENSION(0:NX,0:NY) :: P, P2, GP2, TMP2D
  P = ABS(CP)
  P2 = P*P
  GP2 = G*P2
  TMP2D = -DT*(V + GP2)
  CP = CP*EXP(CI*TMP2D)
END SUBROUTINE NU

SUBROUTINE LUX(CP)
! Solves the partial differential equation only with the X-space
! derivative term using the Crank-Nicholson method
  USE COMM_DATA, ONLY : NX, NY, NXX
  USE CN_DATA, ONLY : CALX, CGAX, CAIPX, CAIPMX, CPX
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
  COMPLEX (8), DIMENSION (NX,0:NY) :: CBE
  COMPLEX (8), DIMENSION (0:NY) :: CXX
  INTEGER :: I
  CBE(1,:) = 0.0D0
  DO I = 1, NXX
     CXX = CP(I,:)-CAIPMX*(CP(I+1,:)-2.0D0*CP(I,:)+CP(I-1,:)) &
              -CPX(I)*(CP(I+1,:)-CP(I-1,:))
     CBE(I+1,:) = CGAX(I)*(CAIPX(I)*CBE(I,:)-CXX)
  END DO
  CP(NX,:) = 0.0D0
  DO I = NX, 1, -1
     CP(I-1,:) = CALX(I)*CP(I,:)+CBE(I,:)
  END DO
END SUBROUTINE LUX

SUBROUTINE LUY(CP)
!  Solves the partial differential equation only with the Y-space
!  derivative term using the Crank-Nicholson method
  USE COMM_DATA, ONLY : NX, NY, NYY
  USE CN_DATA, ONLY : CALY, CGAY, CAIPMY
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
  COMPLEX (8), DIMENSION (0:NX,0:NY) :: CBE
  COMPLEX (8), DIMENSION (0:NX) :: CYY
  INTEGER :: J
  CBE(:,NYY) = CP(:,NY)
  DO J = NYY, 1, -1
     CYY = CP(:,J)-CAIPMY*(CP(:,J+1)-2.0D0*CP(:,J)+CP(:,J-1))
     CBE(:,J-1) = CGAY(J)*(CAIPMY*CBE(:,J)-CYY)
  END DO
  DO J = 0, NYY
     CP(:,J+1) = CALY(J)*CP(:,J)+CBE(:,J)
  END DO
END SUBROUTINE LUY

SUBROUTINE NORM(CP, ZNORM)
!  Calculates the normalization of the wave function and sets it to
!  unity.
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : DX, DY, X, TPI
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(OUT) :: ZNORM
!-----------------------------------------------------------
  INTERFACE
    PURE FUNCTION SIMP(F, DX) RESULT(VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
!-----------------------------------------------------------
  INTEGER :: I
  REAL (8), DIMENSION(0:NX,0:NY) :: P, TMP2D
  REAL (8), DIMENSION(0:NX) :: TMP1D
  P = ABS(CP)
  TMP2D = P*P
  FORALL (I = 0:NX) TMP1D(I) = X(I)*SIMP(TMP2D(I,:), DY)
  ZNORM = SQRT(TPI*SIMP(TMP1D, DX))
  CP = CP/ZNORM
END SUBROUTINE NORM

SUBROUTINE RAD(CP, RHORMS, ZRMS)
 ! Calculates the root mean square size RMS
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : DX, DY, X, X2, Y2, TPI
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(OUT) :: RHORMS, ZRMS
!-----------------------------------------------------------
  INTERFACE
    PURE FUNCTION SIMP(F, DX) RESULT(VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
!-----------------------------------------------------------
  INTEGER :: I
  REAL (8), DIMENSION(0:NX,0:NY) :: P, TMP2D, TMPY1, TMPY2
  REAL (8), DIMENSION(0:NX) :: TMP1D
  P = ABS(CP)
  TMP2D = P*P
  FORALL (I = 0:NX) TMPY1(I,:) = X2(I)*TMP2D(I,:)
  FORALL (I = 0:NY) TMPY2(:,I) = Y2(I)*TMP2D(:,I)
  FORALL (I = 0:NX) TMP1D(I) = X(I)*SIMP(TMPY1(I,:), DY)
  RHORMS = SQRT(TPI*SIMP(TMP1D, DX))
  FORALL (I = 0:NX) TMP1D(I) = X(I)*SIMP(TMPY2(I,:), DY)
  ZRMS = SQRT(TPI*SIMP(TMP1D, DX))
END SUBROUTINE RAD

SUBROUTINE CHEM(CP, MU, EN)
  ! Calculates the chemical potential MU and energy EN.  CP is the wave
  ! function, V is the potential and G is the nonlinearity.
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : DX, DY, X, V, G, TPI
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:,0:), INTENT(IN) :: CP
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
    PURE FUNCTION SIMP(F, DX) RESULT(VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
!-----------------------------------------------------------
  INTEGER :: I
  REAL (8), DIMENSION(0:NX,0:NY) :: DPX, DPY, P, P2, GP2, DP2, TMP2D
  REAL (8), DIMENSION(0:NX) :: TMP1D
  P = ABS(CP)
  FORALL (I=0:NY) DPX(:,I) = DIFF(P(:,I), DX)
  FORALL (I=0:NX) DPY(I,:) = DIFF(P(I,:), DY)
  P2 = P*P
  GP2 = G*P2
  DP2 = DPX*DPX + DPY*DPY
  TMP2D = (V + GP2)*P2 + DP2
  FORALL (I = 0:NX) TMP1D(I) = X(I)*SIMP(TMP2D(I,:), DY)
  MU = TPI*SIMP(TMP1D, DX)
  TMP2D = (V + GP2/2.0D0)*P2 + DP2
  FORALL (I = 0:NX) TMP1D(I) = X(I)*SIMP(TMP2D(I,:), DY)
  EN = TPI*SIMP(TMP1D, DX)
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
