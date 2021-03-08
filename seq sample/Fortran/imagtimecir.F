!  Filename : imagtimecir.F
!  Title : Fortran program for time-dependent Gross-Pitaevskii equation in 
!          a fully anisotropic trap
!  Authors : P. Muruganandam and S. K. Adhikari
!
! Gross-Pitaevskii equation in a two-dimensional circular symmetric trap 
! by imaginary time propagation (FORTRAN 77 version)
!
! L = 0 Case 
!
      PROGRAM MAIN  
      IMPLICIT NONE
      INTEGER N, NX, NPAS, NRUN
! N : Number of space mesh points
      PARAMETER (N = 2000, NX = N-1)
! NPAS : Number of subsequent iterations with fixed nonlinearity, and 
! NRUN : Number of final time iterations with fixed nonlinearity
      PARAMETER (NPAS = 200000, NRUN = 20000)      
! X(0:N) : Space mesh, V(0:N) : Potential, and CP(0:N) : Wave function 
      DOUBLE PRECISION X(0:N), X2(0:N), X3(0:N), V(0:N), CP(0:N)
      DOUBLE PRECISION G, DX, DT, MU, EN, ZNORM
      DOUBLE PRECISION CAIPMX, CAIM
      DOUBLE PRECISION CTT(N), CAIMX(0:N), CAL(N), CGA(N)
      DOUBLE PRECISION G0, XOP, PI, RMS
! OPTION and XOP decides which equation to be solved. 
      INTEGER I, K, OPTION
! DX : Space step and DT : Time step
      DATA DX/0.0025D0/, DT/0.00002D0/
      DATA G0/-2.5097D0/ ! Final nonlinearity
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
      IF(OPTION.EQ.1) XOP=1.D0
      IF(OPTION.EQ.2) XOP=2.D0
      IF(OPTION.EQ.3) XOP=1.D0
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
  902 FORMAT('# Time Stp : ',', NPAS = ',I9,', NRUN = ',I9)
  903 FORMAT('  Nonlinearity G = ', F16.8)
  904 FORMAT('  Space Step DX = ', F12.6, ', Time Step DT = ', F12.6)
      PI=4.D0*ATAN(1.D0)
! INITIALIZE() initializes the starting normalized wave function 'CP' and
! harmonic oscillator potential V
      CALL INITIALIZE(N, OPTION, DX, X, X2, X3, V, CP)
! COEF() defines the coefficients of Crank-Nicholson Scheme.
      CALL COEF(N, NX, DX, DT, CAIPMX, CAIM, X, CTT, CAIMX, CAL, CGA)
      DO I = 0, N ! Writes initial wave function File 1
         WRITE(1,905) X(I), CP(I)
      END DO
  905 FORMAT(F12.6, F16.8)
! NORM() calculates norm and restores normalization
      CALL NORM(N, X, CP, DX, ZNORM,PI)
! CHEM() calculates the chemical potential MU and energy EN.
      CALL CHEM(N, X, V, CP, DX, G, MU, EN)
! RAD() calculates the r.m.s radius RMS
      CALL RAD(N, X3, CP, DX, RMS, PI)
      WRITE (7, 1001)
      WRITE (7, 1002)
      WRITE (7, 1001)
      WRITE (7, 1003) ZNORM, MU/XOP, EN/XOP, RMS, CP(0)
 1001 FORMAT (19X,52('-'))
 1002 FORMAT (20X, 'Norm', 7X, 'Chem', 8X, 'Ener', 7X, '<r>', 7X, 
     &    'Psi(0)')
 1003 FORMAT ('Initial : ', 4X, F11.4, 2F12.6, 3F11.5)
!
      G = XOP*G0
!
      DO K = 1, NPAS ! NPAS iterations transient
! NU() performs time propagation with harmonic potential and nonlinear
! term (non-derivative parts)
         CALL NU(N, V, CP, G, DT)
! LU() performs the time iteration with space derivative alone using
! Crank-Nicholson scheme.
         CALL LU(N, NX, CAIPMX, CTT, CAIMX, CAL, CGA, CP)
         CALL NORM(N, X, CP, DX, ZNORM,PI)
      END DO
      CALL CHEM(N, X, V, CP, DX, G, MU, EN)
      CALL RAD(N, X3, CP, DX, RMS,PI)
      WRITE (7, 1005) ZNORM, MU/XOP, EN/XOP, RMS, CP(0)
 1005 FORMAT('After NPAS iter.:', F8.4, 2F12.6, 3F11.5)
      DO I = 0, N ! Writes intermediate wave funtion in File 2
         WRITE(2,905) X(I), CP(I)
      END DO
!
      DO K = 1, NRUN ! NRUN iterations to check convergence
         CALL NU(N, V, CP, G, DT)
         CALL LU(N, NX, CAIPMX, CTT, CAIMX, CAL, CGA, CP)
         CALL NORM(N, X, CP, DX, ZNORM,PI)
      END DO
      CALL CHEM(N, X, V, CP, DX, G, MU, EN)
      CALL RAD(N, X3, CP, DX, RMS,PI)
      WRITE (7, 1006) ZNORM, MU/XOP, EN/XOP, RMS, CP(0)
      WRITE (7, 1001)
 1006 FORMAT('After NRUN iter.:',F8.4, 2F12.6, 3F11.5)
      DO I = 0, N ! Writes final wave funtion in File 3
         WRITE(3,905) X(I), CP(I)
      END DO
      STOP
      END
!
      SUBROUTINE INITIALIZE(N, OPTION, DX, X, X2, X3, V, CP)
!  Routine that initizlizes the constant and variables.
!  Calculates the potential term V and the initial wave function CP
      IMPLICIT NONE
      INTEGER N, OPTION
      DOUBLE PRECISION DX
      DOUBLE PRECISION X(0:N), X2(0:N), X3(0:N), V(0:N), CP(0:N)
      INTEGER I
      DOUBLE PRECISION PI, PI2, PI4, TMP      
      PI = 4.0D0*ATAN(1.0D0)
      PI4 = SQRT(PI)
      PI2 = SQRT(2.D0*PI*SQRT(2.D0*PI)) ! (2 Pi)^(3/4)
      DO I = 0, N ! Setting up initial condition
         X(I) = DFLOAT(I)*DX
         X2(I) = X(I)*X(I)
         X3(I) = X2(I)*X(I)
         IF(OPTION.EQ.3)THEN
           V(I) = X2(I)/4.D0
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
      SUBROUTINE COEF(N,NX,DX,DT,CAIPMX,CAIM,X,CTT,CAIMX,CAL,CGA)
!   Calculates the coefficients needed in subroutine LU.
      IMPLICIT NONE
      INTEGER N, NX
      DOUBLE PRECISION DX, DT, CAIPMX, CAIM
      DOUBLE PRECISION X(0:N), CTT(N), CAIMX(0:N), CAL(N), CGA(N)
      INTEGER J
      DOUBLE PRECISION DX2, CAI0
      DX2 = DX*DX 
      CAIPMX = -DT/(2.0D0*DX2)
      CAI0 = 1.0D0 + DT/DX2
      CAL(1) = 1.0D0
      CGA(1) = -1.D0/(CAI0+CAIPMX+DT/(4.D0*DX*X(1)))
      DO J = 1, NX
        CTT(J) = -DT/(4.0D0*DX*X(J))
        CAIMX(J) = CAIPMX - CTT(J)
        CAIM = CAIPMX + CTT(J)
        CAL(J+1) = CGA(J)*CAIM
        CGA(J+1) = -1.0D0/(CAI0+(CAIPMX+DT/(4.0D0*DX*X(J+1)))*CAL(J+1))
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
      DOUBLE PRECISION G, DT
      INTEGER I
      DOUBLE PRECISION TMP, P2, P3
      DO I = 0, N
         P2 = CP(I)*CP(I)
         CALL NONLIN(P2,P3,G)
         TMP = DT*(V(I) + P3)
         CP(I) = CP(I)*EXP(-TMP)
      END DO
      RETURN
      END
!
      SUBROUTINE LU(N, NX, CAIPMX, CTT, CAIMX, CAL, CGA, CP)
!  Solves the partial differential equation only with the space
!  derivative term using the Crank-Nicholson method
      IMPLICIT NONE
      INTEGER N, NX
      DOUBLE PRECISION CAIPMX
      DOUBLE PRECISION CAIMX(0:N), CTT(N), CAL(N), CGA(N), CP(0:N)
      INTEGER I
      DOUBLE PRECISION CXX, CBE(N)
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
      SUBROUTINE NORM(N, X, CP, DX, ZNORM, PI)
!  Calculates the normalization of the wave function and sets it to
!  unity.
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION X(0:N), CP(0:N), DX, ZNORM
      INTEGER I
      DOUBLE PRECISION CP2(0:N)
      DOUBLE PRECISION SIMP,PI
      EXTERNAL SIMP 
      DO I = 0, N
         CP2(I) = CP(I)*CP(I)*X(I)
      END DO
      ZNORM = SQRT(2.D0*PI*SIMP(N, CP2, DX))
      DO I = 0, N
         CP(I) = CP(I)/ZNORM
      END DO
      RETURN
      END
!
      SUBROUTINE RAD(N, X3, CP, DX, RMS, PI)
!  Calculates the root mean square radius RMS.
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION CP(0:N), X3(0:N), DX, RMS, PI
      INTEGER I
      DOUBLE PRECISION CP2(0:N), TMP1D(0:N)
      DOUBLE PRECISION SIMP
      EXTERNAL SIMP 
      DO I = 0, N
         CP2(I) = CP(I)*CP(I)
         TMP1D(I) = CP2(I)*X3(I)
      END DO
      RMS = SQRT(2.0D0*PI*SIMP(N, TMP1D, DX))
      RETURN
      END
!
      SUBROUTINE CHEM(N, X, V, CP, DX, G, MU, EN)
!  Calculates the chemical potential MU and energy EN.  CP is the wave
!  function defined at N+1 points at steps of DX, V is the potential and
!  G is the nonlinearity.
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION X(0:N), V(0:N), CP(0:N)
      DOUBLE PRECISION DX, G, MU, EN
      INTEGER I
      DOUBLE PRECISION P2, DP2, DCP(0:N), TMP1D(0:N), TMP2D(0:N)
      DOUBLE PRECISION SIMP, P3, PI, P4
      EXTERNAL SIMP
      PI=4.D0*ATAN(1.D0)
      CALL DIFF(N, CP, DCP, DX)
      !TMP1D(0) = DCP(1)*DCP(1)
      DO I = 0, N
         P2 = CP(I)*CP(I)
         P4 = CP(I)*CP(I)
         DP2 = DCP(I)*DCP(I)
         CALL NONLIN(P4,P3,G)
         TMP1D(I) = (V(I)+P3)*P2 + DP2
         TMP2D(I) = (V(I)+P3/2.0D0)*P2 + DP2
      END DO
      DO I = 0, N
        TMP1D(I)=X(I)*TMP1D(I)
        TMP2D(I)=X(I)*TMP2D(I)
      END DO
      MU = 2.D0*PI*SIMP(N, TMP1D, DX)
      EN = 2.D0*PI*SIMP(N, TMP2D, DX)
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
! Richardsonextrapolation formula. The derivative at the
! boundaries are assumed to be zero
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION P(0:N), DP(0:N), DX
      INTEGER I
      DP(0)=0.D0
      DO I = 1, 2
      DP(i) = (P(i+1) - P(i-1))/(2.0d0*dx)
      END DO
      DO I = 3, N-3
         DP(I) = (P(i-2)-8.0d0*P(i-1)+8.0d0*P(i+1)
     &             -P(i+2))/(12.0d0*dx)
      END DO
      DO I = N-2, N-1
         DP(i) = (P(i+1) - P(i-1))/(2.0d0*dx)
      END DO
      DP(N) = 0.0D0
      RETURN
      END

