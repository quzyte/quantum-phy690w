!  Filename : imagtime1d.F
!  Title : Fortran program for time-dependent Gross-Pitaevskii equation in 
!          a fully anisotropic trap
!  Authors : P. Muruganandam and S. K. Adhikari
! 
! Gross-Pitaevskii equation in a one-dimensional parabolic trap by 
! imaginary time propagation (FORTRAN 77 version)
!
      PROGRAM MAIN 
      IMPLICIT NONE
      INTEGER N, N2, NX, NPAS, NRUN
! N : Number of space mesh points
      PARAMETER (N = 8000, N2 = N/2, NX = N-1)
! NPAS : Number of subsequent iterations with fixed nonlinearity, and 
! NRUN : Number of final time iterations with fixed nonlinearity
      PARAMETER (NPAS = 200000, NRUN = 20000)
! X(0:N) : Space mesh, V(0:N) : Potential, and CP(0:N) : Wave function 
      DOUBLE PRECISION X(0:N), X2(0:N), V(0:N), CP(0:N)
      DOUBLE PRECISION G, DX, DT, MU, ZNORM, RMS, EN
      DOUBLE PRECISION CAIPM, CAL(0:NX), CGA(0:NX)
! OPTION and XOP decides which equation to be solved. 
      DOUBLE PRECISION G0, XOP, GSTP
      INTEGER I, K, OPTION
! DX : Space step and DT : Time step
       DATA DX/0.0025D0/, DT/0.00002D0/
       DATA G0/62.742D0/ ! Final nonlinearity
! Subroutine INTIALIZE() used to initialize the space mesh X(I), 
! potential V(I) and the initial wave function. Subroutine COEF() used to
! generate the coefficients for the Crank-Nicholson Scheme. The routine
! NU() performs time progation for the non-derivative part and LU() performs
! time propagation of derivative part. NORM() calculates the norm and 
! normalizes the wave function, CHEM() and RAD() are used to calculate the 
! chemical potential, energy and the rms radius, respectively. The subroutine 
! NONLIN() calculates the nonlinear term, and functions DIFF(), used to
! calculate the space derivatives of the wave function used in CHEM(). The
! function SIMP() does the integration by Simpson's rule.
!      OPTION=1 ! Solves -psi_xx+V(x)psi+G0|psi|^2 psi=i psi_t
      OPTION=2 !Solves [-psi_xx+V(x)psi]/2+G0|psi|^2 psi=i psi_t
!      OPTION=3 ! Solves -psi_xx+V(x)psi/4+G0|psi|^2 psi=i psi_t
      IF (OPTION.EQ.1) XOP = 1.D0
      IF (OPTION.EQ.2) XOP = 2.D0
      IF (OPTION.EQ.3) XOP = 1.D0
      G = 0.D0
      WRITE(7,900) OPTION
      WRITE(7,*)  
      WRITE(7,901) N
      WRITE(7,902) NPAS, NRUN  
      WRITE(7,903) G0 
      WRITE(7,904) DX, DT 
      WRITE(7,*)
  900 FORMAT('  OPTION = ', I3)
  901 FORMAT('# Space Stp N = ', I8)
  902 FORMAT('# Time Stp :  NPAS = ',I9,', NRUN = ',I9)
  903 FORMAT('  Nonlinearity G = ', F16.8)
  904 FORMAT('  Space Step DX = ', F12.6, ', Time Step DT = ', F12.6)
  905 FORMAT(F12.6, F16.8)
! INITIALIZE() initializes the starting normalized wave function 'CP' and
! harmonic oscillator potential V
      CALL INITIALIZE(N, N2,OPTION, X,X2, V, CP, DX)
! COEF() defines the coefficients of Crank-Nicholson Scheme.
      CALL COEF(NX, DX, DT, CAIPM, CAL, CGA)
      WRITE (1,905)(X(I),CP(I),I=0,N)!Writes initial wave function in File  1
! NORM() calculates norm and restores normalization
      CALL NORM(N, CP, DX, ZNORM)
! CHEM() calculates the chemical potential MU and energy EN.
      CALL CHEM(N, V, CP, DX, G, MU, EN)
! RAD() calculates the r.m.s radius RMS
      CALL RAD(N,X2, CP, DX, RMS)
      WRITE (7, 1001)
      WRITE (7, 1002)
      WRITE (7, 1001)
      WRITE (7, 1003) ZNORM, MU/XOP, EN/XOP, RMS, CP(N2)
 1001 FORMAT (19X,52('-'))
 1002 FORMAT (20X, 'Norm', 7X, 'Chem', 8X, 'Ener', 7X, '<x>', 7X, 
     &    'Psi(0)')
 1003 FORMAT ('Initial : ', 4X, F11.4, 2F12.6, 3F11.5)
!
      G = XOP*G0
!
      DO K = 1, NPAS  ! NPAS iterations transient
! NU() performs time propagation with harmonic potential and nonlinear
! term (non-derivative parts)
         CALL NU(N, V, CP, G, DT)
! LU() performs the time iteration with space derivative alone using
! Crank-Nicholson scheme.
         CALL LU(N, NX, CAIPM, CAL, CGA, CP)
         CALL NORM(N, CP, DX, ZNORM)
      END DO
      CALL CHEM(N, V, CP, DX, G, MU, EN)
      CALL RAD(N,X2, CP, DX, RMS)
      WRITE (7, 1005) ZNORM, MU/XOP, EN/XOP, RMS, CP(N2)
 1005 FORMAT('After NPAS iter.:', F8.4, 2F12.6, 3F11.5)
      !Writes intermediate wavefuntion in File 2
      WRITE(2,905) (X(I),CP(I),I=0,N) 
!
      DO K = 1, NRUN ! NRUN iterations to check convergence
         CALL NU(N, V, CP, G, DT)
         CALL LU(N, NX, CAIPM, CAL, CGA, CP)
         CALL NORM(N, CP, DX, ZNORM)
      END DO
      CALL CHEM(N, V, CP, DX, G, MU, EN)
      CALL RAD(N,X2, CP, DX, RMS)
      WRITE (7, 1006) ZNORM, MU/XOP, EN/XOP, RMS, CP(N2)
      WRITE (7, 1001)
 1006 FORMAT('After NRUN iter.:',F8.4, 2F12.6, 3F11.5)
      WRITE (3,905)(X(I), CP(I),I=0,N)! Writes final wavefuntion in File 3
      STOP
      END
!
      SUBROUTINE INITIALIZE(N, N2, OPTION, X, X2, V, CP, DX)
!  Routine that initizlizes the constant and variables.
!  Calculates the potential term V and the initial wave function CP
      IMPLICIT NONE
      INTEGER N, N2, OPTION
      DOUBLE PRECISION X(0:N), X2(0:N), V(0:N), CP(0:N)
      DOUBLE PRECISION DX      
      DOUBLE PRECISION PI, PI4, TMP,PI2
      INTEGER I      
      PI = 4.0D0*ATAN(1.0D0)
      PI4 = SQRT(SQRT(PI))
      PI2 = SQRT(SQRT(2.D0*PI))
      DO I = 0, N ! Setting up initial condition
         X(I) = DFLOAT(I-N2)*DX
         X2(I) = X(I)*X(I)
         IF(OPTION.EQ.3)THEN
           V(I) = X2(I)/4.D0
           TMP = V(I)
           CP(I) = EXP(-TMP)/PI2
         ELSE
           V(I) = X2(I)
           TMP = X2(I)/2.0D0
           CP(I) = EXP(-TMP)/PI4
         END IF
      END DO
      RETURN
      END
!
      SUBROUTINE COEF(NX, DX, DT, CAIPM, CAL, CGA)
!   Calculates the coefficients needed in subroutine LU.
      IMPLICIT NONE
      INTEGER NX
      DOUBLE PRECISION DX, DT, CAIPM, CAL(0:NX), CGA(0:NX)      
      INTEGER J
      DOUBLE PRECISION DX2, CAI0    
      DX2 = DX*DX
      CAIPM = -DT/(2.0D0*DX2)
      CAI0 = 1.0D0 + DT/DX2
      CAL(NX) = 0.0D0
      CGA(NX) = -1.0D0/CAI0
      DO J = NX, 1, -1
         CAL(J-1) = CGA(J)*CAIPM
         CGA(J-1) = -1.0D0/(CAI0+CAIPM*CAL(J-1))
      END DO
      RETURN
      END
!
      SUBROUTINE NU(N, V, CP, G, DT)
!   Solves the partial differential equation with the potential and the
!   nonlinear term.
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION V(0:N), CP(0:N)
      DOUBLE PRECISION G, DT, P3      
      INTEGER I
      DOUBLE PRECISION TMP, TMP2
      DO I = 1, N
         TMP = CP(I)*CP(I)
            CALL NONLIN(TMP,P3,G)
         TMP2 = DT*(V(I) + P3)
         CP(I) = CP(I)*EXP(-TMP2)
      END DO
      RETURN
      END
!
      SUBROUTINE LU(N, NX, CAIPM, CAL, CGA, CP)
!  Solves the partial differential equation only with the space
!  derivative term using the Crank-Nicholson method
      IMPLICIT NONE
      INTEGER N, NX
      DOUBLE PRECISION CAIPM, CAL(0:NX), CGA(0:NX), CP(0:N)
      INTEGER I
      DOUBLE PRECISION CBE(0:NX), CXX
      CBE(NX) = CP(N)
      DO I = NX, 1, -1
         CXX = CP(I)-(CP(I+1)-2.0D0*CP(I)+CP(I-1))*CAIPM
         CBE(I-1) = CGA(I)*(CAIPM*CBE(I)-CXX)
      END DO
      CP(0) = 1.0D-16
      DO I = 0, NX
         CP(I+1) = CAL(I)*CP(I)+CBE(I)
      END DO      
      RETURN
      END
!
      SUBROUTINE NORM(N, CP, DX, ZNORM)
!  Calculates the normalization of the wave function and sets it to
!  unity.
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION CP(0:N), DX, ZNORM
      INTEGER I
      DOUBLE PRECISION CP2(0:N)
      DOUBLE PRECISION SIMP
      EXTERNAL SIMP 
      DO I = 0, N
         CP2(I) = CP(I)*CP(I)
      END DO
      ZNORM = SQRT(SIMP(N, CP2, DX))
      DO I = 0, N
         CP(I) = CP(I)/ZNORM
      END DO
      RETURN
      END
!
      SUBROUTINE RAD(N,X2, CP, DX, RMS)
!  Calculates the root mean square size RMS
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION CP(0:N), DX, RMS
      INTEGER I
      DOUBLE PRECISION CP2(0:N),X2(0:N)
      DOUBLE PRECISION SIMP
      EXTERNAL SIMP 
      DO I = 0, N
         CP2(I) = CP(I)*CP(I)*X2(I)
      END DO
      RMS = SQRT(SIMP(N, CP2, DX))
      RETURN
      END
!
      SUBROUTINE CHEM(N, V, CP, DX, G, MU, EN)
!  Calculates the chemical potential MU and energy EN.  CP is the wave
!  function defined at N+1 points at steps of DX, V is the potential and
!  G is the nonlinearity.
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION V(0:N), CP(0:N), DX, G, MU, EN      
      INTEGER I
      DOUBLE PRECISION P2, DP2, DCP(0:N), TMP1D(0:N), TMP2D(0:N)
      DOUBLE PRECISION SIMP,P3
      EXTERNAL SIMP 
      CALL DIFF(N, CP, DCP, DX)
      DO I = 0, N
         P2 = CP(I)*CP(I)
         DP2 = DCP(I)*DCP(I)
         CALL NONLIN(P2,P3,G)
         TMP1D(I) = (V(I) + P3)*P2 + DP2
         TMP2D(I) = (V(I) + P3/2.D0)*P2 + DP2
      END DO
      MU = SIMP(N, TMP1D, DX)
      EN = SIMP(N, TMP2D, DX)
      RETURN
      END 
!
      DOUBLE PRECISION FUNCTION SIMP(N, F, DX)
! Does the spatial integration with Simpson's rule.
! N refer to the number of integration points, DX space step, and
! F is the function to be integrated.
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION F(0:N), DX
      DOUBLE PRECISION F1, F2
      INTEGER I
      F1 = F(1) + F(N-1)  ! N even
      F2 = F(2)
      DO I = 3, N-3, 2
         F1 = F1 + F(I)
         F2 = F2 + F(I+1)
      END DO
      SIMP = DX*(F(0) + 4.0D0*F1 + 2.0D0*F2 + F(N))/3.0D0
      RETURN
      END
!
      SUBROUTINE NONLIN(P2,P3,G)
! Calculates the nonlinear term
      IMPLICIT NONE
      DOUBLE PRECISION P2,P3,G
      P3=G*P2
      RETURN
      END
!
      SUBROUTINE DIFF(N, P, DP, DX)
! Computes the first derivative DP of P using
! Richardson extrapolation formula. The derivative at the  
! boundaries are assumed to be zero
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION P(0:N), DP(0:N), DX
      INTEGER I
      DP(0) = 0.0D0
      DP(1) = (P(2) - P(0))/(2.0D0*DX)
      DO I = 2, N-2
        DP(I) = (P(I-2)-8.0D0*P(I-1)+8.0D0*P(I+1)-P(I+2))/(12.0D0*DX)
      END DO
      DP(N-1) = (P(N) - P(N-2))/(2.0D0*DX)
      DP(N) = 0.0D0
      RETURN      
      END     
