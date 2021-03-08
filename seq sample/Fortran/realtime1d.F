!  Filename : realtime1d.F
!  Title : Fortran program for time-dependent Gross-Pitaevskii equation in 
!          a fully anisotropic trap
!  Authors : P. Muruganandam and S. K. Adhikari
! 
! Gross-Pitaevskii equation in a one-dimensional parabolic trap by 
! real time propagation (FORTRAN 77 version)
!
      PROGRAM MAIN 
      IMPLICIT NONE
      INTEGER N, N2, NX, NPAS, NRUN,NSTP
! N : Number of space mesh points
      PARAMETER (N = 5000, N2 = N/2, NX = N-1)
! NSTP : Number of iterations to introduce the nonlinearity
! NPAS : Number of subsequent iterations with fixed nonlinearity, and 
! NRUN : Number of final time iterations with fixed nonlinearity
      PARAMETER (NSTP = 1000000, NPAS = 1000, NRUN = 40000)
! X(0:N) : Space mesh, V(0:N) : Potential, and CP(0:N) : Wave function  
      DOUBLE PRECISION X(0:N), X2(0:N), V(0:N)
      DOUBLE COMPLEX CI, CP(0:N)
      DOUBLE PRECISION G, DX, DT, MU, EN, ZNORM, RMS
      DOUBLE COMPLEX CAIPM, CAL(0:NX), CGA(0:NX)
! OPTION and XOP decides which equation to be solved. 
      DOUBLE PRECISION G0, GSTP, XOP
      INTEGER I, K, OPTION
! DX : Space step and DT : Time step
      DATA DX/0.01D0/, DT/0.0001D0/
      DATA G0/62.742D0/   ! Final nonlinearity
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
      WRITE(7,902) NSTP, NPAS, NRUN  
      WRITE(7,903) G0 
      WRITE(7,904) DX, DT 
      WRITE(7,*)
  900 FORMAT('  OPTION = ', I3)
  901 FORMAT('# Space Stp N = ', I8)
  902 FORMAT('# Time Stp : NSTP', I9,' , NPAS = ', I9,', NRUN = ', I9)
  903 FORMAT('  Nonlinearity G = ', F16.8)
  904 FORMAT('  Space Step DX = ', F12.6, ', Time Step DT = ', F12.6)
  905 FORMAT(F12.6, F16.8)
! INITIALIZE() initializes the starting normalized wave function 'CP' and
! harmonic oscillator potential V
      CALL INITIALIZE(N, N2, OPTION, DX, X, X2, V, CI, CP)
! COEF() defines the coefficients of Crank-Nicholson Scheme.
      CALL COEF(NX, DX, DT, CI, CAIPM, CAL, CGA)
      ! Writes initial wave function in File  1
      WRITE (1,905) (X(I), ABS(CP(I)), I = 0, N)
! NORM() calculates norm and restores normalization
      CALL NORM(N, CP, DX, ZNORM)
! CHEM() calculates the chemical potential MU and energy EN.
      CALL CHEM(N, V, CP, DX, G, MU, EN)
! RAD() calculates the r.m.s radius RMS
      CALL RAD(N, X2, CP, DX, RMS)
      WRITE (7, 1001)
      WRITE (7, 1002)
      WRITE (7, 1001)
      WRITE (7, 1003) ZNORM, MU/XOP, EN/XOP, RMS, ABS(CP(N2))
 1001 FORMAT (19X,52('-'))
 1002 FORMAT (20X, 'Norm', 7X, 'Chem', 8X, 'Ener', 7X, '<x>', 7X, 
     &    'Psi(0)')
 1003 FORMAT ('Initial : ', 4X, F11.4, 2F12.3, 3F11.3)
      GSTP = XOP*G0/DFLOAT(NSTP)
      DO K = 1, NSTP ! Introduces nonlinearity in 'NSTP' iterations
         G = G + GSTP
! NU() performs time propagation with harmonic potential and nonlinear
! term (non-derivative parts)
         CALL NU(N, V, CP, CI, G, DT) 
! LU() performs the time iteration with space derivative alone using
! Crank-Nicholson scheme.
         CALL LU(N, NX, CAIPM, CAL, CGA, CP) ! 
      END DO
      CALL NORM(N, CP, DX, ZNORM)
      CALL CHEM(N, V, CP, DX, G, MU, EN)
      CALL RAD(N, X2, CP, DX, RMS)
      WRITE (7, 1004) ZNORM, MU/XOP, EN/XOP, RMS, ABS(CP(N2))
 1004 FORMAT('After NSTP iter.:', F8.4, 2F12.3, 3F11.3)
      ! Write intermediate wavefuntion in File 2
      WRITE(2,905) (X(I), ABS(CP(I)), I = 0, N)  
!
      DO K = 1, NPAS ! NPAS iterations transient
         CALL NU(N, V, CP, CI, G, DT)
         CALL LU(N, NX, CAIPM, CAL, CGA, CP)
      END DO
      CALL NORM(N, CP, DX, ZNORM)
      CALL CHEM(N, V, CP, DX, G, MU, EN)
      CALL RAD(N, X2, CP, DX, RMS)
      WRITE (7, 1005) ZNORM, MU/XOP, EN/XOP, RMS, ABS(CP(N2))
 1005 FORMAT('After NPAS iter.:', F8.4, 2F12.3, 3F11.3)
      WRITE (7, 1001)
      ! Writes final wavefuntion in File 3
      WRITE (3,905) (X(I), ABS(CP(I)), I=0,N)
!
!  The following line defines a nonstationary problem which is studied
!  below and the time evolution written on file 8 via WRITE(8,*)
!
      G = 0.5D0*G                                      
      DO K = 1, NRUN ! NRUN iterations to study nonlinear dynamics
         CALL NU(N, V, CP, CI, G, DT)
         CALL LU(N, NX, CAIPM, CAL, CGA, CP)
         IF (MOD(K,100).EQ.0) THEN          
            CALL RAD(N, X2, CP, DX, RMS)                  
            WRITE(8,905) DFLOAT(K)*DT*XOP, RMS     
         END IF                                   
      END DO
      STOP
      END
!
      SUBROUTINE INITIALIZE(N, N2, OPTION, DX, X, X2, V, CI, CP)
!  Calculates the potential term V and the initial wave function CP
      IMPLICIT NONE
      INTEGER N, N2, OPTION
      DOUBLE PRECISION DX, X(0:N), X2(0:N), V(0:N)
      DOUBLE COMPLEX CI, CP(0:N)      
      DOUBLE PRECISION PI, PI2, PI4, TMP
      INTEGER I
      CI = CMPLX(0.0D0, 1.0D0)
      PI = 4.0D0*ATAN(1.0D0)
      PI2 = SQRT(SQRT(2.0D0*PI))
      PI4 = SQRT(SQRT(PI))
      DO I = 0, N ! Setting up initial condition
         X(I) = DFLOAT(I-N2)*DX
         X2(I) = X(I)*X(I)
         IF (OPTION.EQ.3) THEN
            V(I) = X2(I)/4.0D0
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
      SUBROUTINE NU(N, V, CP, CI, G, DT)
!   Solves the partial differential equation with the potential and the 
!   nonlinear term.  
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION V(0:N)
      DOUBLE COMPLEX CP(0:N), CI
      DOUBLE PRECISION G, DT,P3      
      INTEGER I
      DOUBLE PRECISION TMP, TMP2
      DO I = 1, N
         TMP = ABS(CP(I))
         TMP2 = TMP*TMP
          CALL NONLIN(TMP2,P3,G)
         TMP = DT*(V(I)+P3)
         CP(I) = CP(I)*EXP(-CI*TMP)
      END DO
      RETURN
      END
!
      SUBROUTINE NONLIN(P2,P3,G)
!  Calculates the nonlinear term.
      IMPLICIT NONE
      DOUBLE PRECISION P2,P3,G
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
      CP(0) = CMPLX(0.0D0, 0.0D0)
      DO I = 0, NX
         CP(I+1) = CAL(I)*CP(I)+CBE(I)
      END DO
      RETURN
      END
!
      SUBROUTINE NORM(N, CP, DX, ZNORM)
!   Calculates the normalization of the wave function and sets it to 
!   unity. 
      IMPLICIT NONE
      INTEGER N
      DOUBLE COMPLEX CP(0:N)
      DOUBLE PRECISION DX, ZNORM
      INTEGER I
      DOUBLE PRECISION CP2(0:N)
      DOUBLE PRECISION SIMP
      EXTERNAL SIMP 
      DO I = 0, N
         CP2(I) = ABS(CP(I))**2
      END DO
      ZNORM = SQRT(SIMP(N, CP2, DX))
      DO I = 0, N
         CP(I) = CP(I)/ZNORM
      END DO
      RETURN
      END
!
      SUBROUTINE RAD(N, X2, CP, DX, RMS)
!  Calculates the root mean square size RMS.
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION X2(0:N)
      DOUBLE COMPLEX CP(0:N)
      DOUBLE PRECISION DX, RMS
      INTEGER I
      DOUBLE PRECISION CP2(0:N)
      DOUBLE PRECISION SIMP
      EXTERNAL SIMP 
      DO I = 0, N
         CP2(I) = ABS(CP(I))**2*X2(I)
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
      DOUBLE PRECISION V(0:N)
      DOUBLE COMPLEX CP(0:N)
      DOUBLE PRECISION DX, G, MU, EN
      INTEGER I
      DOUBLE PRECISION P2, DP2, GP2, TMP1D(0:N), EMP1D(0:N)
      DOUBLE PRECISION DCP(0:N), XCP(0:N)
      DOUBLE PRECISION SIMP
      EXTERNAL SIMP 
      DO I=0,N
      XCP(I)= ABS(CP(I))
      END DO
      CALL DIFF(N, DX, XCP, DCP)
      DO I = 0, N
         P2 = XCP(I)**2
         DP2 = DCP(I)**2
         CALL NONLIN(P2,GP2,G)
         TMP1D(I) = (V(I) + GP2)*P2 + DP2
         EMP1D(I) = (V(I) + GP2/2.0D0)*P2 + DP2
      END DO
      MU = SIMP(N, TMP1D, DX)
      EN = SIMP(N, EMP1D, DX)
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
      DOUBLE PRECISION P(0:N), DP(0:N)
      DOUBLE PRECISION DX
      INTEGER I
      DP(0) = 0.0D0
      DP(1) = (P(2) - P(0))/(2.0D0*DX)
      DO I = 2, N-2
        DP(I) = (P(I-2)-8.0D0*P(I-1)+8.0D0*P(I+1)
     c   -P(I+2))/(12.0D0*DX)
      END DO
      DP(N-1) = (P(N) - P(N-2))/(2.0D0*DX)
      DP(N) = 0.0D0
      RETURN      
      END     
