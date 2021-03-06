!  Filename : realtimeaxial.F
!  Title : Fortran program for time-dependent Gross-Pitaevskii equation in 
!          a fully anisotropic trap
!  Authors : P. Muruganandam and S. K. Adhikari
! 
!  Fortran program for Gross-Pitaevskii equation in axially symmetric 
!  trap by real time propagation (FORTRAN 77 Version)
!
!
      IMPLICIT NONE
      INTEGER NX, NY, NXX, NYY, NY2, NSTP, NPAS, NRUN
! NX, NY : Number of space mesh points (X and Y)  
      PARAMETER (NX = 130, NXX = NX-1)
      PARAMETER (NY = 130, NYY = NY-1, NY2 = NY/2)
! NSTP : Number of iterations to introduce the nonlinearity,
! NPAS : Number of subsequent iterations with fixed nonlinearity, and 
! NRUN : Number of final time iterations with fixed nonlinearity
      PARAMETER (NSTP = 100000, NPAS = 1000, NRUN = 20000)
! X(0:N), Y(0:N) : Space mesh, V(0:N) : Potential, and 
! CP(0:N) : Wave function (COMPLEX)
      DOUBLE PRECISION X(0:NX), Y(0:NY), X2(0:NX), Y2(0:NY)
      DOUBLE PRECISION V(0:NX,0:NY)
      DOUBLE COMPLEX CP(0:NX,0:NY)
      DOUBLE PRECISION DX, DY, DT
!
      DOUBLE COMPLEX CALX(0:NX), CGAX(0:NX), CAIPX(0:NX), CPX(1:NX)
      DOUBLE COMPLEX CALY(0:NY), CGAY(0:NY)
      DOUBLE COMPLEX CI, CAIPMX, CAIPMY
!
      DOUBLE PRECISION G, GSTP, MU, EN, ZNORM, RHORMS, ZRMS
! G0 : Final nonlinearity, KAP, LAM : Anisotrophy coefficients
! OPTION and XOP decides which equation to be solved. 
      DOUBLE PRECISION G0, KAP, LAM, XOP, PI
      INTEGER I, J, K, OPTION
!       
! DX, DY : Space step and DT : Time step
      DATA DX/0.1D0/, DY/0.1D0/, DT/0.001D0/
      DATA G0/18.81d0/, KAP/1.0D0/, LAM/4.0D0/
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
      !OPTION=1 ! Solves -psi_xx-psi_yy+V(x,y)psi+G0|psi|^2 psi=i psi_t
      OPTION=2 !Solves [-psi_xx-psi_yy+V(x,y)psi]/2+G0|psi|^2 psi=i psi_t
      !OPTION=3 ! Solves -psi_xx-psi_yy+V(x,y)psi/4+G0|psi|^2 psi=i psi_t
      IF(OPTION.EQ.1) XOP = 1.D0  
      IF(OPTION.EQ.2) XOP = 2.D0  
      IF(OPTION.EQ.3) XOP = 1.D0  
      G = 0.0D0  ! NONLINEARITY COEFFICIENT
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
! INITIALIZE() initializes the starting normalized wave function 'CP' and
! harmonic oscillator potential V
!
      CALL INITIALIZE(NX, NY, OPTION, DX, DY, LAM, KAP, PI, X, Y, X2, 
     & Y2, V, CP)
! COEF() defines the coefficients of Crank-Nicholson Scheme.
      CALL COEF(NX, NY, DX, DY, DT, X, CI, CAIPMX, CAIPMY, CAIPX, 
     $  CPX, CALX, CGAX, CALY, CGAY)
!      
! NORM() calculates norm and restores normalization
      CALL NORM(NX, NY, DX, DY, ZNORM, PI, X, CP)
! CHEM() calculates the chemical potential MU and energy EN.
      CALL CHEM(NX, NY, DX, DY, G, MU, EN, PI, X, V, CP)
! RAD() calculates the r.m.s radius RMS
      CALL RAD(NX, NY, DX, DY, RHORMS, ZRMS, PI, X, X2, Y2, CP)
      WRITE (7, 1001)
      WRITE (7, 1002)
      WRITE (7, 1001)
      WRITE (7, 1003) ZNORM, MU/XOP, EN/XOP, RHORMS, ZRMS, 
     &  ABS(CP(0,NY2+1))
 1001 FORMAT (12X,55('-'))
 1002 FORMAT(14X,'Norm',4X,'Chem',5X,'Energy',4X,'<rho>',5X,'<z>',
     &  7X,'psi(0)')
 1003 FORMAT ('Initial :  ',F6.3,5(F10.3))
      DO I = 0, NX, 4 ! Writes initial wave function File 1
         DO J = 0, NY, 4
            WRITE (1,999) X(I), Y(J), ABS(CP(I,J))
         END DO
         WRITE (1,*)
      END DO 
  999 FORMAT (2F12.6, F16.8)
!
      GSTP = XOP*G0/NSTP
      DO K = 1, NSTP ! Introduces nonlinearity in 'NSTP' iterations
         G = G + GSTP
! NU() performs time propagation with harmonic potential and nonlinear
! term (non-derivative parts)
         CALL NU(NX, NY, DT, G, V, CI, CP) 
! LUX() and LUY() perform the time iteration with space derivative in
! X and Y, respectively, using Crank-Nicholson scheme.
         CALL LUX(NX, NY, CP, CAIPX, CPX, CALX, CGAX, CAIPMX)
         CALL LUY(NX, NY, CP, CALY, CGAY, CAIPMY)
      END DO
!     
      CALL NORM(NX, NY, DX, DY, ZNORM, PI, X, CP)
      CALL CHEM(NX, NY, DX, DY, G, MU, EN, PI, X, V, CP)
      CALL RAD(NX, NY, DX, DY, RHORMS, ZRMS, PI, X, X2, Y2, CP)
      WRITE (7, 1004) ZNORM, MU/XOP, EN/XOP, RHORMS, ZRMS, 
     &  ABS(CP(0,NY2+1))
 1004 FORMAT ('NSTP iter :',F6.3,  5(F10.3))
!
      DO I = 0, NX, 4 ! Writes intermediate wave funtion in File 2
         DO J = 0, NY, 4
            WRITE(2,999) X(I), Y(J), ABS(CP(I,J))
         END DO
         WRITE(2,*)
      END DO 
      DO K = 1, NPAS! NPAS iterations transient
         CALL NU(NX, NY, DT, G, V, CI, CP) 
         CALL LUX(NX, NY, CP, CAIPX, CPX, CALX, CGAX, CAIPMX)
         CALL LUY(NX, NY, CP, CALY, CGAY, CAIPMY)
      END DO
!     
      CALL NORM(NX, NY, DX, DY, ZNORM, PI, X, CP)
      CALL CHEM(NX, NY, DX, DY, G, MU, EN, PI, X, V, CP)
      CALL RAD(NX, NY, DX, DY, RHORMS, ZRMS, PI, X, X2, Y2, CP)
      WRITE (7, 1005) ZNORM, MU/XOP, EN/XOP, RHORMS, ZRMS, 
     &  ABS(CP(0,NY2+1))
      WRITE (7, 1001)
 1005 FORMAT ('NPAS iter :',F6.3,  5(F10.3))
!
      DO I = 0, NX, 4 ! Writes final wavefuntion in File 3
         DO J = 0, NY, 4
            WRITE(3,999) X(I), Y(J), ABS(CP(I,J))
         END DO
         WRITE(3,*)
      END DO 
!
!  The following line defines a nonstationary problem which is studied
!  below and the time evolution written on file 8 via WRITE(8,*)
!
      G = 0.5D0*G           
      DO K = 1, NRUN ! NRUN iterations to study nonlinear dynamics
         CALL NU(NX, NY, DT, G, V, CI, CP)
         CALL LUX(NX, NY, CP, CAIPX, CPX, CALX, CGAX, CAIPMX)
         CALL LUY(NX, NY, CP, CALY, CGAY, CAIPMY)
         IF (MOD(K,100).EQ.0) THEN          
            CALL RAD(NX, NY, DX, DY, RHORMS, ZRMS, PI, X, X2, Y2, CP)
            WRITE(8,1999) DFLOAT(K)*DT*XOP, RHORMS, ZRMS
         END IF                                   
      END DO
 1999 FORMAT(F12.6, 2F16.8)
      STOP
      END
!
      SUBROUTINE INITIALIZE(NX, NY, OPTION, DX, DY, LAM, KAP, PI, X, Y, 
     & X2, Y2, V, CP)
!   Routine that initizlizes the constant and variables.
!   Calculates the potential term V and the initial wave function CP
      IMPLICIT NONE
      INTEGER NX, NY, OPTION
      DOUBLE PRECISION DX, DY, LAM, KAP, PI
      DOUBLE PRECISION X(0:NX), Y(0:NY), X2(0:NX), Y2(0:NY)
      DOUBLE PRECISION V(0:NX,0:NY)
      DOUBLE COMPLEX CP(0:NX,0:NY)
      DOUBLE PRECISION XLY(0:NX,0:NY)

      INTEGER I, J, NY2
      DOUBLE PRECISION FAC, LAM2, KAP2, PI3, FAC3
    
      NY2 = NY/2
      LAM2 = LAM*LAM
      KAP2 = KAP*KAP

      PI3 = PI**3
!
      FAC = SQRT(SQRT(LAM*KAP2/(PI3)))
      FAC3 = FAC/SQRT(SQRT(8.D0))
!
      DO I = 0, NX
         X(I) = I*DX
         X2(I) = X(I)*X(I)
      END DO
      DO J = 0, NY
         Y(J) = (J-NY2)*DY
         Y2(J) = Y(J)*Y(J)
      END DO
!
      IF (OPTION .EQ. 3) THEN ! Setting up initial condition
        DO I = 0, NX
          DO J = 0, NY
             V(I,J) = (KAP2*X2(I) + LAM2*Y2(J))/4.0D0
             XLY(I,J) = (KAP*X2(I) + LAM*Y2(J))/4.0D0
             CP(I,J) = FAC3*EXP(-XLY(I,J))
          END DO
        END DO
      ELSE
        DO I = 0, NX
          DO J = 0, NY
              V(I,J) = (KAP2*X2(I) + LAM2*Y2(J))
              XLY(I,J) = (KAP*X2(I) + LAM*Y2(J))/2.0D0
              CP(I,J) = FAC*EXP(-XLY(I,J))
          END DO
        END DO
      END IF
!
      END
!
      SUBROUTINE COEF(NX, NY, DX, DY, DT, X, CI, CAIPMX, CAIPMY,   
     $ CAIPX, CPX, CALX, CGAX, CALY, CGAY)
!   Calculates the coefficients needed in subroutine LUX and LUY.
      IMPLICIT NONE
      INTEGER NX, NY
      DOUBLE PRECISION DX, DY, DT, X(0:NX)
      DOUBLE COMPLEX CI, CAIPMX, CAIPMY
      DOUBLE COMPLEX CALX(0:NX), CGAX(0:NX), CAIPX(0:NX), CPX(1:NX)
      DOUBLE COMPLEX CALY(0:NY), CGAY(0:NY)
C
      INTEGER NYY, NXX, J
      DOUBLE PRECISION DX2, DY2
      DOUBLE COMPLEX CAI0, CAIM, CDT
C     
      NXX = NX - 1       
      DX2 = DX*DX
      CDT = DT/CI
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
      NYY = NY - 1
      DY2 = DY*DY
      CAIPMY = CDT/(2.0D0*DY2)
      CAI0 = 1.D0 - CDT/DY2
      CALY(NYY) = 0.0D0
      CGAY(NYY) = -1.0D0/CAI0
!
      DO J = NYY, 1, -1
         CALY(J-1) = CGAY(J)*CAIPMY
         CGAY(J-1) = -1.0D0/(CAI0+CAIPMY*CALY(J-1))
      END DO
!
      END 
!
      SUBROUTINE NU(NX, NY, DT, G, V, CI, CP)  ! Exact solution
!   Solves the partial differential equation with the potential and the
!   nonlinear term.
      IMPLICIT NONE
      INTEGER NX, NY
      DOUBLE PRECISION DT, G,P3
      DOUBLE PRECISION V(0:NX, 0:NY)
      DOUBLE COMPLEX CI, CP(0:NX, 0:NY)
C
      INTEGER I, J
      DOUBLE PRECISION P2, TMP
C
      DO I = 0, NX
         DO J = 0, NY
            TMP = ABS(CP(I,J))
            P2 = TMP*TMP
            CALL NONLIN(P2,P3,G)
            TMP = -DT*(V(I,J)+P3)
            CP(I,J) = CP(I,J)*EXP(CI*TMP)
         END DO
      END DO
C      
      RETURN 
      END
!
      SUBROUTINE LUX(NX, NY, CP, CAIPX, CPX, CALX, CGAX, CAIPMX)
!  Solves the partial differential equation only with the X-space
!  derivative term using the Crank-Nicholson method
      IMPLICIT NONE
      INTEGER NX, NY
      DOUBLE COMPLEX CP(0:NX,0:NY)
      DOUBLE COMPLEX CALX(0:NX), CGAX(0:NX) 
      DOUBLE COMPLEX CAIPX(0:NX), CPX(NX)
      DOUBLE COMPLEX CAIPMX

      DOUBLE COMPLEX CBE(NX), CXX
      INTEGER I, J
!
      DO J = 0, NY            
         CBE(1) = 0.D0
         DO I = 1,(NX-1)
            CXX = CP(I,J)-CAIPMX*(CP(I+1,J)-2.0D0*CP(I,J)+CP(I-1,J))
     $            -CPX(I)*(CP(I+1,J)-CP(I-1,J))
            CBE(I+1) = CGAX(I)*(CAIPX(I)*CBE(I)-CXX)
         END DO
C     
         CP(NX,J)=0.D0
         DO I = NX,1,-1
            CP(I-1,J) = CALX(I)*CP(I,J)+CBE(I)
         END DO
      END DO
C      
      RETURN 
      END
!
      SUBROUTINE LUY(NX, NY, CP, CALY, CGAY, CAIPMY)
!  Solves the partial differential equation only with the Y-space
!  derivative term using the Crank-Nicholson method
      IMPLICIT NONE
      INTEGER NX, NY
      DOUBLE COMPLEX CP(0:NX,0:NY)
      DOUBLE COMPLEX CALY(0:NY), CGAY(0:NY)
      DOUBLE COMPLEX CAIPMY
C
      INTEGER I, J
      DOUBLE COMPLEX CBE(0:NY), CYY
C
      DO I = 0, NX
         CBE(NY-1) = CP(I, NY)
         DO J = NY-1, 1, -1
            CYY = CP(I,J)-CAIPMY*(CP(I,J+1)-2.0D0*CP(I,J)+CP(I,J-1))
            CBE(J-1) = CGAY(J)*(CAIPMY*CBE(J)-CYY)
         END DO
C     
         DO J = 0, NY-1
            CP(I,J+1) = CALY(J)*CP(I,J)+CBE(J)
         END DO
      END DO
C      
      RETURN 
      END
!
      SUBROUTINE NORM(NX, NY, DX, DY, ZNORM, PI, X, CP)
!  Calculates the normalization of the wave function and sets it to
!  unity.
      IMPLICIT NONE
      INTEGER NX, NY
      DOUBLE PRECISION DX, DY, ZNORM, PI
      DOUBLE PRECISION X(0:NX)
      DOUBLE COMPLEX CP(0:NX,0:NY)
C
      DOUBLE PRECISION TMPX(0:NX), TMPY(0:NY)
      INTEGER I, J
C
      DOUBLE PRECISION SIMP
      EXTERNAL SIMP 
!
      DO I = 0, NX
         DO J = 0, NY
            TMPY(J) = ABS(CP(I,J))**2
         END DO
         TMPX(I) = SIMP(NY, TMPY, DY)*X(I)
      END DO
!
      ZNORM = SQRT(2.0D0*PI*SIMP(NX, TMPX, DX))
!
      DO I = 0, NX
         DO J = 0, NY
            CP(I,J) = CP(I,J)/ZNORM
         END DO
      END DO
!
      RETURN
      END
!
      SUBROUTINE RAD(NX, NY, DX, DY, RHORMS, ZRMS, PI, X, X2, Y2, CP)
!  Calculates the root mean square size RMS
      IMPLICIT NONE
      INTEGER NX, NY
      DOUBLE PRECISION DX, DY, RHORMS, ZRMS, PI
      DOUBLE PRECISION X(0:NX), X2(0:NX), Y2(0:NY)
      DOUBLE COMPLEX CP(0:NX,0:NY)
!
      DOUBLE PRECISION TMPX1(0:NX), TMPY1(0:NY), TMPY(0:NY)
      DOUBLE PRECISION TMPX2(0:NX), TMPY2(0:NY)
      INTEGER I, J
!
      DOUBLE PRECISION SIMP
      EXTERNAL SIMP 
!
      DO I = 0, NX
         DO J = 0, NY
            TMPY(J) = ABS(CP(I,J))**2
            TMPY1(J) = TMPY(J)*X2(I)
            TMPY2(J) = TMPY(J)*Y2(J)
          END DO
          TMPX1(I) = SIMP(NY, TMPY1, DY)*X(I)
          TMPX2(I) = SIMP(NY, TMPY2, DY)*X(I)
      END DO
!
      RHORMS = SQRT(2.0D0*PI*SIMP(NX, TMPX1, DX))
      ZRMS = SQRT(2.0D0*PI*SIMP(NX, TMPX2, DX))
!
      RETURN
      END
!
      SUBROUTINE CHEM(NX, NY, DX, DY, G, MU, EN, PI, X, V, CP)
!  Calculates the chemical potential MU and energy EN.  CP is the wave
!  function defined at N+1 points at steps of DX, DY, V is the potential 
!  and G is the nonlinearity.
      IMPLICIT NONE
      INTEGER NX, NY
      DOUBLE PRECISION DX, DY, G, MU, PI,EN
      DOUBLE PRECISION X(0:NX), V(0:NX,0:NY)
      DOUBLE COMPLEX CP(0:NX,0:NY)
      
      INTEGER I, J
      DOUBLE PRECISION DCPX(0:NX), DCPY(0:NY), TMPX(0:NX), TMPY(0:NY)
      DOUBLE PRECISION EMPX(0:NX), EMPY(0:NY)
      DOUBLE PRECISION PX1(0:NX,0:NY), PY1(0:NX,0:NY), P2(0:NX,0:NY)
      DOUBLE PRECISION DP2, P3
C
      DOUBLE PRECISION SIMP
      EXTERNAL SIMP
!
      DO J = 0, NY
         DO I = 0, NX
            TMPX(I) = ABS(CP(I,J))
            P2(I,J) = TMPX(I)*TMPX(I)
         END DO
         CALL DIFF(NX, TMPX, DCPX, DX)
         DO I = 0, NX
            PX1(I,J) = DCPX(I)
         END DO
      END DO
!
      DO I = 0, NX
         DO J = 0, NY
            TMPY(J)=ABS(CP(I,J))
         END DO
         CALL DIFF(NY, TMPY, DCPY, DY)
         DO J = 0, NY
            PY1(I,J) = DCPY(J)
         END DO
      END DO
!
      DO I = 0, NX
         DO J = 0, NY
            CALL NONLIN(P2(I,J), P3, G)
            DP2 = PX1(I,J)*PX1(I,J)+PY1(I,J)*PY1(I,J)
            TMPY(J) = (V(I,J) + P3)*P2(I,J) + DP2
            EMPY(J) = (V(I,J) + P3/2.0D0)*P2(I,J) + DP2
         END DO
         TMPX(I) = SIMP(NY, TMPY, DY)*X(I)
         EMPX(I) = SIMP(NY, EMPY, DY)*X(I)
      END DO
      MU = 2.0D0*PI*SIMP(NX, TMPX, DX)
      EN = 2.0D0*PI*SIMP(NX, EMPX, DX)
!
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
! Computes the first derivative DP of P using  Richardson extrapolation
! formula. The derivative at the boundaries are assumed to be zero
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
