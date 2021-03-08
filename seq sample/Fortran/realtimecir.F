!  Filename : realtimecir.F
!  Title : Fortran program for time-dependent Gross-Pitaevskii equation in 
!          a fully anisotropic trap
!  Authors : P. Muruganandam and S. K. Adhikari
! 
! Gross-Pitaevskii equation in a two-dimensional circular symmetric trap 
! by real time propagation (FORTRAN 77 version)
!
! L = 0 Case 
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
      DOUBLE PRECISION X(0:N), X2(0:N), X3(0:N), V(0:N)
      DOUBLE COMPLEX CI, CP(0:N)
      DOUBLE PRECISION G, DX, DT, ZNORM, RMS, MU, EN
      DOUBLE COMPLEX CAIPMX, CAIM, CTT(N), CAIMX(0:N), CAL(N), CGA(N)
! OPTION and XOP decides which equation to be solved. 
      DOUBLE PRECISION G0, GSTP, PI, XOP
      INTEGER I, K, OPTION
! DX : Space step and DT : Time step
      DATA DX/0.01D0/, DT/0.0001D0/
      DATA G0/12.5484D0/ ! Final nonlinearity
      DATA CI/(0.0D0,1.0D0)/
      DATA PI/3.14159265358979D0/
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
      OPTION=2 !Solves [-psi_xx+V(x)psi]/2+G0|psi|^2 psi=i psi_t
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
      CALL INITIALIZE(N, OPTION, DX, PI, X, X2, X3, V, CP)
! COEF() defines the coefficients of Crank-Nicholson Scheme.
      CALL COEF(N,NX,DX,DT,X,CI,CAIPMX,CAIM,CTT,CAIMX,CAL,CGA)
      !CALL COEF(NX, DX, DT, CI, CAIPM, CAL, CGA)
      DO I = 0, N ! Writes initial wave function File 1
         WRITE(1,999) X(I), ABS(CP(I))
      END DO
  999 FORMAT(F12.6, F16.8)
! NORM() calculates norm and restores normalization
      CALL NORM(N, DX, ZNORM, PI, X, CP)
! CHEM() calculates the chemical potential MU and energy EN.
      CALL CHEM(N, X, V, DX, G, PI, MU, EN, CP)
! RAD() calculates the r.m.s radius RMS
      CALL RAD(N, X3, DX, RMS, PI, CP)
      WRITE (7, 1001)
      WRITE (7, 1002)
      WRITE (7, 1001)
      WRITE (7, 1003) ZNORM, MU/XOP, EN/XOP, RMS, ABS(CP(0))
 1001 FORMAT (19X,52('-'))
 1002 FORMAT (20X, 'Norm', 7X, 'Chem', 8X, 'Ener', 7X, '<r>', 7X, 
     &    'Psi(0)')
 1003 FORMAT ('Initial : ', 4X, F11.4, 2F12.3, 3F11.3)
      GSTP = XOP*G0/DFLOAT(NSTP)
      DO K = 1, NSTP ! Introduces nonlinearity in 'NSTP' iterations
         G = G + GSTP
! NU() performs time propagation with harmonic potential and nonlinear
! term (non-derivative parts)
         CALL NU(N, V, G, DT, CI, CP)
! LU() performs the time iteration with space derivative alone using
! Crank-Nicholson scheme.
         CALL LU(N, NX, CAIPMX, CTT, CAIMX, CAL, CGA, CP)
      END DO
      CALL NORM(N, DX, ZNORM, PI, X, CP)
      CALL CHEM(N, X, V, DX, G, PI, MU, EN, CP)
      CALL RAD(N, X3, DX, RMS, PI, CP)
      WRITE (7, 1004) ZNORM, MU/XOP, EN/XOP, RMS, ABS(CP(0))
 1004 FORMAT('After NSTP iter.:', F8.4, 2F12.3, 3F11.3)
      DO I = 0, N ! Writes intermediate wave funtion in File 2
         WRITE(2,999) X(I), ABS(CP(I))
      END DO
!
      DO K = 1, NPAS ! NPAS iterations transient
         CALL NU(N, V, G, DT, CI, CP)
         CALL LU(N, NX, CAIPMX, CTT, CAIMX, CAL, CGA, CP)
      END DO
      CALL NORM(N, DX, ZNORM, PI, X, CP)
      CALL CHEM(N, X, V, DX, G, PI, MU, EN, CP)
      CALL RAD(N, X3, DX, RMS, PI, CP)
      WRITE (7, 1005) ZNORM, MU/XOP, EN/XOP, RMS, ABS(CP(0))
 1005 FORMAT('After NPAS iter.:', F8.4, 2F12.3, 3F11.3)
      WRITE (7, 1001)
      DO I = 0, N ! Writes final wave funtion in File 3
         WRITE(3,999) X(I), ABS(CP(I))
      END DO
!  The following line defines a nonstationary problem which is studied
!  below and the time evolution written on file 8 via WRITE(8,*)
      G = 0.5D0*G                     
      DO K = 1, NRUN ! NRUN iterations to study nonlinear dynamics
         CALL NU(N, V, G, DT, CI, CP)
         CALL LU(N, NX, CAIPMX, CTT, CAIMX, CAL, CGA, CP)
         IF (MOD(K,100).EQ.0) THEN    
            CALL RAD(N, X3, DX, RMS, PI, CP)
            WRITE(8,999) DFLOAT(K)*DT*XOP, RMS 
         END IF                               
      END DO
      STOP
      END
!
      SUBROUTINE INITIALIZE(N, OPTION, DX, PI, X, X2, X3, V, CP)
!  Routine that initializes the constant and variables.
!  Calculates the potential term V and the initial wave function CP
      IMPLICIT NONE
      INTEGER N, OPTION
      DOUBLE PRECISION DX, PI, X(0:N), X2(0:N), X3(0:N), V(0:N)
      DOUBLE COMPLEX CP(0:N)
      DOUBLE PRECISION PI2, PI4, TMP
      INTEGER I
      PI4 = SQRT(PI)
      PI2 = SQRT(2.D0*PI*SQRT(2.D0*PI)) ! (2 Pi)^(3/4)
      DO I = 0, N ! Setting up initial condition
         X(I) = DFLOAT(I)*DX
         X2(I) = X(I)*X(I)
         X3(I) = X2(I)*X(I)
         IF (OPTION.EQ.3) THEN
           V(I) = X2(I)/4.0D0
           TMP = V(I)
           CP(I) = EXP(-TMP)/PI2
         ELSE
           V(I) = X2(I)
           TMP = V(I)/2.0D0
           CP(I) = EXP(-TMP)/PI4
         END IF
      END DO
      RETURN
      END
!
      SUBROUTINE COEF(N, NX, DX, DT, X, CI, CAIPMX, CAIM, CTT, CAIMX,  
     & CAL, CGA)
!   Calculates the coefficients needed in subroutine LU.
      IMPLICIT NONE
      INTEGER N, NX
      DOUBLE PRECISION DX, DT, X(0:N)
      DOUBLE COMPLEX CI, CAIPMX, CAIM
      DOUBLE COMPLEX CTT(N), CAIMX(0:N), CAL(N), CGA(N)      
      INTEGER J
      DOUBLE PRECISION DX2
      DOUBLE COMPLEX CDT, CAI0
      CDT = DT/CI       ! Generating the Coefficients for the
      DX2 = DX*DX       ! C-N method (to solve the spatial part)
      CAIPMX = CDT/(2.0D0*DX2)
      CAI0 = 1.0D0 - CDT/DX2
      CAL(1) = 1.0D0
      CGA(1) = -1.D0/(CAI0+CAIPMX-CDT/(4.D0*DX*X(1)))
      DO J = 1, NX
        CTT(J) = CDT/(4.0D0*DX*X(J))
        CAIMX(J) = CAIPMX - CTT(J)
        CAIM = CAIPMX + CTT(J)
        CAL(J+1) = CGA(J)*CAIM
        CGA(J+1) = -1.0D0/(CAI0+(CAIPMX-CDT/(4.0D0*DX*X(J+1)))*CAL(J+1))
      END DO
      END
!
      SUBROUTINE NU(N, V, G, DT, CI, CP)
!   Solves the partial differential equation with the potential and the
!   nonlinear term.
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION V(0:N), G, DT
      DOUBLE COMPLEX CI, CP(0:N)
      INTEGER I
      DOUBLE PRECISION TMP, P2, GP2
      DO I = 0, N
         TMP = ABS(CP(I))
         P2 = TMP*TMP
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
      SUBROUTINE LU(N, NX, CAIPMX, CTT, CAIMX, CAL, CGA, CP)
!  Solves the partial differential equation only with the space
!  derivative term using the Crank-Nicholson method
      IMPLICIT NONE
      INTEGER N, NX
      DOUBLE COMPLEX CAIPMX
      DOUBLE COMPLEX CAIMX(0:N), CTT(N), CAL(N), CGA(N), CP(0:N)
      INTEGER I
      DOUBLE COMPLEX CXX, CBE(N)
      CBE(1) = 0.D0
      DO I = 1, NX
         CXX = CP(I)-(CP(I+1)-2.0D0*CP(I)+CP(I-1))*CAIPMX
     $            -CTT(I)*(CP(I+1)-CP(I-1))            
         CBE(I+1) = CGA(I)*(CAIMX(I)*CBE(I)-CXX)
      END DO
      CP(N) = 0.0D0
      DO I = N, 1, -1
         CP(I-1) = CAL(I)*CP(I) + CBE(I)
      END DO      
      RETURN
      END
!
      SUBROUTINE NORM(N, DX, ZNORM, PI, X, CP)
!  Calculates the normalization of the wave function and sets it to
!  unity.
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION DX, ZNORM, PI, X(0:N)
      DOUBLE COMPLEX CP(0:N)
      INTEGER I
      DOUBLE PRECISION TMP, CP2(0:N)
      DOUBLE PRECISION SIMP
      EXTERNAL SIMP 
      DO I = 0, N
         TMP = ABS(CP(I))
         CP2(I) = TMP*TMP*X(I)
      END DO
      ZNORM = SQRT(2.D0*PI*SIMP(N, CP2, DX))
      DO I = 0, N
         CP(I) = CP(I)/ZNORM
      END DO
      RETURN
      END
!
      SUBROUTINE RAD(N, X3, DX, RMS, PI, CP)
!  Calculates the root mean square radius RMS
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION X3(0:N), DX, RMS, PI
      DOUBLE COMPLEX CP(0:N)
      INTEGER I
      DOUBLE PRECISION TMP, CP2(0:N)
      DOUBLE PRECISION SIMP
      EXTERNAL SIMP 
      DO I = 0, N
         TMP = ABS(CP(I))
         CP2(I) = TMP*TMP*X3(I)
      END DO
      RMS = SQRT(2.D0*PI*SIMP(N, CP2, DX))
      RETURN
      END
!
      SUBROUTINE CHEM(N, X, V, DX, G, PI, MU, EN, CP)
!  Calculates the chemical potential MU and energy EN.  CP is the wave
!  function defined at N+1 points at steps of DX, V is the potential and
!  G is the nonlinearity.
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION X(0:N), V(0:N), DX, G, PI, MU, EN
      DOUBLE COMPLEX CP(0:N)
      INTEGER I
      DOUBLE PRECISION P2, GP2, DP2
      DOUBLE PRECISION  TMP1D(0:N), EMP1D(0:N), DCP(0:N), XCP(0:N)
      DOUBLE PRECISION SIMP
      EXTERNAL SIMP 
      DO I = 0,N
         XCP(I) = ABS(CP(I))
      END DO
      CALL DIFF(N, DX, XCP, DCP)
      DO I = 0, N
         P2 = XCP(I)*XCP(I)
         DP2 = DCP(I)*DCP(I)
         CALL NONLIN(P2, GP2, G)
         TMP1D(I) = (V(I) + GP2)*P2 + DP2
         EMP1D(I) = (V(I) + GP2/2.0D0)*P2 + DP2
      END DO
      DO I = 0, N
        TMP1D(I)=X(I)*TMP1D(I)
        EMP1D(I)=X(I)*EMP1D(I)
      END DO
      MU = 2.0D0*PI*SIMP(N, TMP1D, DX)
      EN = 2.0D0*PI*SIMP(N, EMP1D, DX)
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

