!  Filename : realtimesph.F
!  Title : Fortran program for time-dependent Gross-Pitaevskii equation in 
!          a fully anisotropic trap
!  Authors : P. Muruganandam and S. K. Adhikari
! 
! Gross-Pitaevskii equation in a three-dimensional spherical trap by 
! real time propagation (FORTRAN 77 version)
!
      PROGRAM MAIN 
      IMPLICIT NONE
      INTEGER N, NX, NSTP, NPAS, NRUN
! N : Number of space mesh points
      PARAMETER (N = 2000, NX = N-1)
! NSTP : Number of iterations to introduce the nonlinearity
! NPAS : Number of subsequent iterations with fixed nonlinearity, and 
! NRUN : Number of final time iterations with fixed nonlinearity
      PARAMETER (NSTP = 1000000, NPAS = 1000, NRUN = 40000)      
! X(0:N) : Space mesh, V(0:N) : Potential
! CP(0:N) : Wave function (COMPLEX)
      DOUBLE PRECISION X(0:N), X2(0:N), V(0:N)
      DOUBLE COMPLEX CI, CP(0:N)
      DOUBLE PRECISION G, DX, DT, MU, EN, ZNORM
      DOUBLE COMPLEX CAIPM, CAL(0:NX), CGA(0:NX)
! OPTION and XOP decides which equation to be solved. 
      DOUBLE PRECISION G0, GSTP, PI, RMS, XOP
      INTEGER I, K, OPTION
! DX : Space step and DT : Time step
      DATA DX/0.01D0/, DT/0.0001D0/
      DATA G0/125.484D0/   ! Final nonlinearity
      DATA CI/(0.0D0,1.0D0)/
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
      !OPTION=1 ! Solves -psi_xx+V(x)psi+G0|psi|^2 psi=i psi_t
      OPTION=2 ! Solves [-psi_xx+V(x)psi]/2+G0|psi|^2 psi=i psi_t
      !OPTION=3 ! Solves -psi_xx+V(x)psi/4+G0|psi|^2 psi=i psi_t
      IF(OPTION.EQ.1) XOP = 1.D0
      IF(OPTION.EQ.2) XOP = 2.D0
      IF(OPTION.EQ.3) XOP = 1.D0
      WRITE(7,900) OPTION
      WRITE(7,*)  
      WRITE(7,901) N
      WRITE(7,902) NSTP, NPAS, NRUN  
      WRITE(7,903) G0 
      WRITE(7,904) DX, DT 
      WRITE(7,*)
  900 FORMAT('  OPTION = ', I3)
  901 FORMAT('# Space Stp N = ', I8)
  902 FORMAT('# Time Stp : NSTP = ', I9, ', NPAS = ',I9,', NRUN = ',I9)
  903 FORMAT('  Nonlinearity G = ', F16.8)
  904 FORMAT('  Space Step DX = ', F12.6, ', Time Step DT = ', F12.6)
! INITIALIZE() initializes the starting normalized wave function 'CP' and
! harmonic oscillator potential V
      CALL INITIALIZE(N, OPTION, X, X2, V, DX, PI, CP)
! COEF() defines the coefficients of Crank-Nicholson Scheme.
      CALL COEF(NX, DX, DT, CI, CAIPM, CAL, CGA)
      DO I = 0, N ! Writes initial wave function File 1
         IF(I.NE.0) WRITE(1,905) X(I), ABS(CP(I))/X(I)
         IF(I.EQ.0) WRITE(1,905) X(I), ABS(CP(1))/X(1)
      END DO
  905 FORMAT(F12.6, F16.8)
! NORM() calculates norm and restores normalization
      CALL NORM(N, CP, DX, ZNORM, PI)
! CHEM() calculates the chemical potential MU and energy EN.
      CALL CHEM(N, X2, V, DX, G, MU, EN, PI, CP)
! RAD() calculates the r.m.s radius RMS
      CALL RAD(N, X2, DX, RMS, PI, CP)
      WRITE (7, 1001)
      WRITE (7, 1002)
      WRITE (7, 1001)
      WRITE (7, 1003) ZNORM, MU/XOP, EN/XOP, RMS, ABS(CP(1))/X(1)
 1001 FORMAT (19X,52('-'))
 1002 FORMAT (20X, 'Norm', 7X, 'Chem', 8X, 'Ener', 7X, '<r>', 7X, 
     &    'Psi(0)')
 1003 FORMAT ('Initial : ', 4X, F11.4, 2F12.3, 3F11.3)
      GSTP = XOP*G0/DFLOAT(NSTP)
      DO K = 1, NSTP ! Introduces nonlinearity in 'NSTP' iterations
         G = G + GSTP
! NU() performs time propagation with harmonic potential and nonlinear
! term (non-derivative parts)
         CALL NU(N, X2, V, G, DT, CP, CI)
! LU() performs the time iteration with space derivative alone using
! Crank-Nicholson scheme.
         CALL LU(N, NX, CAIPM, CAL, CGA, CP)
      END DO
      CALL NORM(N, CP, DX, ZNORM,PI)
      CALL CHEM(N, X2, V, DX, G, MU, EN, PI, CP)
      CALL RAD(N, X2, DX, RMS, PI, CP)
      WRITE (7, 1004) ZNORM, MU/XOP, EN/XOP, RMS, ABS(CP(1))/X(1)
 1004 FORMAT('After NSTP iter.:', F8.4, 2F12.3, 3F11.3)
      DO I = 0, N ! Writes intermediate wave funtion in File 2
         IF (I.NE.0) WRITE(2,905) X(I), ABS(CP(I))/X(I)
         IF (I.EQ.0) WRITE(2,905) X(I), ABS(CP(1))/X(1)
      END DO
!
      DO K = 1, NPAS ! NPAS iterations transient
         CALL NU(N, X2, V, G, DT, CP, CI)
         CALL LU(N, NX, CAIPM, CAL, CGA, CP)
      END DO
      CALL NORM(N, CP, DX, ZNORM,PI)
      CALL CHEM(N, X2, V, DX, G, MU, EN, PI, CP)
      CALL RAD(N, X2, DX, RMS, PI, CP)
      WRITE (7, 1005) ZNORM, MU/XOP, EN/XOP, RMS, ABS(CP(1))/X(1)
 1005 FORMAT('After NPAS iter.:', F8.4, 2F12.3, 3F11.3)
      WRITE (7, 1001)
      DO I = 0, N ! Writes final wave funtion in File 3
         IF (I.NE.0) WRITE(3,905) X(I), ABS(CP(I))/X(I)
         IF (I.EQ.0) WRITE(3,905) X(I), ABS(CP(1))/X(1)
      END DO
!  The following line defines a nonstationary problem which is studied
!  below and the time evolution written on file 8 via WRITE(8,*)
      G = 0.5D0*G                     
      DO K = 1, NRUN ! NRUN iterations to study nonlinear dynamics
         CALL NU(N, X2, V, G, DT, CP, CI)
         CALL LU(N, NX, CAIPM, CAL, CGA, CP)
         IF (MOD(K,100).EQ.0) THEN    
            CALL RAD(N, X2, DX, RMS, PI, CP)
            WRITE(8,905) DFLOAT(K)*DT*XOP, RMS 
         END IF                               
      END DO
      STOP
      END
!
      SUBROUTINE INITIALIZE(N, OPTION, X, X2, V, DX, PI, CP)
!  Routine that initializes the constant and variables.
!  Calculates the potential term V and the initial wave function CP
      IMPLICIT NONE
      INTEGER N, OPTION
      DOUBLE PRECISION X(0:N), X2(0:N), V(0:N), DX, PI
      DOUBLE COMPLEX CP(0:N)
      DOUBLE PRECISION PI2, PI4, TMP
      INTEGER I
      PI = 4.D0*ATAN(1.D0)
      PI2 = SQRT(2.0D0*PI*SQRT(2.0D0*PI)) ! (2*PI)^(3/4)
      PI4 = SQRT(PI*SQRT(PI)) ! (PI)^(3/4)
      DO I = 0, N ! Setting up initial condition
         X(I) = DFLOAT(I)*DX
         X2(I) = X(I)*X(I)
         IF (OPTION.EQ.3) THEN
           V(I) = X2(I)/4.0D0
           TMP = V(I)
           CP(I) = X(I)*EXP(-TMP)/PI2
         ELSE
           V(I) = X2(I)
           TMP = V(I)/2.0D0
           CP(I) = X(I)*EXP(-TMP)/PI4
         END IF
      END DO
      RETURN
      END
!
      SUBROUTINE COEF(NX, DX, DT, CI, CAIPM, CAL, CGA)
!   Calculates the coefficients needed in subroutine LU.
      IMPLICIT NONE
      INTEGER NX
      DOUBLE PRECISION DX, DT
      DOUBLE COMPLEX CI, CAIPM, CAL(0:NX), CGA(0:NX)      
      INTEGER J
      DOUBLE PRECISION DX2
      DOUBLE COMPLEX CDT, CAI0
      CDT = DT/CI       ! Generating the Coefficients for the
      DX2 = DX*DX       ! C-N method (to solve the spatial part)
      CAIPM = CDT/2.0D0/DX2
      CAI0 = 1.D0-CDT/DX2
      CAL(NX) = 0.0D0
      CGA(NX) = -1.0D0/CAI0
      DO J = NX, 1, -1
         CAL(J-1) = CGA(J)*CAIPM
         CGA(J-1) = -1.0D0/(CAI0+CAIPM*CAL(J-1))
      END DO
      RETURN
      END
!
      SUBROUTINE NU(N, X2, V, G, DT, CP, CI)
!   Solves the partial differential equation with the potential and the
!   nonlinear term.
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION X2(0:N), V(0:N), G, DT
      DOUBLE COMPLEX CP(0:N), CI
      INTEGER I
      DOUBLE PRECISION TMP, P2, GP2
      DO I = 1, N
         TMP = ABS(CP(I))
         P2 = TMP*TMP/X2(I)
         CALL NONLIN(P2, GP2, G)
         TMP = DT*(V(I)+GP2)
         CP(I) = CP(I)*EXP(-CI*TMP)
      END DO
      RETURN
      END
!
      SUBROUTINE NONLIN(P2, P3, G)
! Calculates the nonlinear term.
      IMPLICIT NONE
      DOUBLE PRECISION P2, P3, G
      P3=G*P2
      RETURN
      END
!
      SUBROUTINE LU(N, NX, CAIPM, CAL, CGA, CP)
!  Solves the partial differential equation only with the space
!  derivative term using the Crank-Nicholson method
      IMPLICIT NONE
      INTEGER N, NX
      DOUBLE COMPLEX CAIPM, CAL(0:NX), CGA(0:NX), CP(0:N)
      INTEGER I
      DOUBLE COMPLEX CBE(0:NX), CXX
      CBE(NX) = CP(N)
      DO I = NX, 1, -1
         CXX = CP(I)-(CP(I+1)-2.0D0*CP(I)+CP(I-1))*CAIPM
         CBE(I-1) = CGA(I)*(CAIPM*CBE(I)-CXX)
      END DO
      CP(0) = 0.0D0
      DO I = 0, NX
         CP(I+1) = CAL(I)*CP(I)+CBE(I)
      END DO
      RETURN
      END
!
      SUBROUTINE NORM(N, CP, DX, ZNORM, PI)
!  Calculates the normalization of the wave function and sets it to
!  unity.
      IMPLICIT NONE
      INTEGER N
      DOUBLE COMPLEX CP(0:N)
      DOUBLE PRECISION DX, ZNORM, PI
      INTEGER I
      DOUBLE PRECISION TMP, CP2(0:N)
      DOUBLE PRECISION SIMP
      EXTERNAL SIMP 
      DO I = 0, N
         TMP = ABS(CP(I))
         CP2(I) = TMP*TMP
      END DO
      ZNORM = SQRT(4.0D0*PI*SIMP(N, CP2, DX))
      DO I = 0, N
         CP(I) = CP(I)/ZNORM
      END DO
      RETURN
      END
!
      SUBROUTINE RAD(N, X2, DX, RMS, PI, CP)
!  Calculates the root mean square radius RMS
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION X2(0:N), DX, RMS, PI
      DOUBLE COMPLEX CP(0:N)
      INTEGER I
      DOUBLE PRECISION TMP, CP2(0:N)
      DOUBLE PRECISION SIMP
      EXTERNAL SIMP 
      DO I = 0, N
         TMP = ABS(CP(I))
         CP2(I) = TMP*TMP*X2(I)
      END DO
      RMS = SQRT(4.D0*PI*SIMP(N, CP2, DX))
      RETURN
      END
!
      SUBROUTINE CHEM(N, X2, V, DX, G, MU, EN, PI, CP)
!  Calculates the chemical potential MU and energy EN.  CP is the wave
!  function defined at N+1 points at steps of DX, V is the potential and
!  G is the nonlinearity.
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION X2(0:N), V(0:N), DX, G, MU, EN, PI
      DOUBLE COMPLEX CP(0:N)
      INTEGER I
      DOUBLE PRECISION P2, P22, GP22, DP2, TMP1D(0:N), EMP1D(0:N)
      DOUBLE PRECISION  DCP(0:N), XCP(0:N)
      DOUBLE PRECISION SIMP
      EXTERNAL SIMP 
      DO I = 0,N
         XCP(I) = ABS(CP(I))
      END DO
      CALL DIFF(N, DX, XCP, DCP)
      DO I = 1, N
         P2 = XCP(I)**2
         DP2 = DCP(I)**2
         P22 = P2/X2(I)
         CALL NONLIN(P22, GP22, G)
         TMP1D(I) = (V(I) + GP22)*P2 + DP2
         EMP1D(I) = (V(I) + GP22/2.0D0)*P2 + DP2
      END DO
      TMP1D(0) = TMP1D(1)
      EMP1D(0) = EMP1D(1)
      MU = 4.0D0*PI*SIMP(N, TMP1D, DX)
      EN = 4.0D0*PI*SIMP(N, EMP1D, DX)
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
      SUBROUTINE DIFF(N, DX, P, DP)
! Computes the first derivative DP of P using
! Richardson extrapolation formula. The derivative at the
! boundaries are assumed to be zero
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION DX
      DOUBLE PRECISION P(0:N), DP(0:N)
      INTEGER I
      DP(0)=0.0D0
      DO I = 1, 2
      DP(I) = (P(I+1) - P(I-1))/(2.0D0*DX)
      END DO
      DO I = 3, N-3
      DP(I) = (P(I-2)-8.0D0*P(I-1)+8.0D0*P(I+1)-P(I+2))/(12.0D0*DX)
      END DO
      DO I = N-2, N-1
         DP(I) = (P(I+1) - P(I-1))/(2.0D0*DX)
      END DO
      DP(N) = 0.0D0
      RETURN
      END

