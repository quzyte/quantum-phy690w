!  Filename : realtime3d.F
!
!  Title : Fortran program for time-dependent Gross-Pitaevskii equation in 
!          a fully anisotropic trap
!  Authors : P. Muruganandam and S. K. Adhikari
! 
!  Fortran program for Gross-Pitaevskii equation in three-dimensional 
!  anisotropic trap by real time propagation (FORTRAN 77 Version)
!
      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER NX, NXX, NX2, NY, NYY, NY2, NZ, NZZ, NZ2
      INTEGER NSTP, NPAS, NRUN
! NX, NY, NZ : Number of space mesh points (X, Y and Z)  
      PARAMETER (NX=200, NXX = NX-1, NX2 = NX/2)
      PARAMETER (NY=160, NYY = NY-1, NY2 = NY/2)
      PARAMETER (NZ=120, NZZ = NZ-1, NZ2 = NZ/2)
! NSTP : Number of iterations to introduce the nonlinearity
! NPAS : Number of subsequent iterations with fixed nonlinearity, and 
! NRUN : Number of final time iterations with fixed nonlinearity
      PARAMETER (NSTP = 60000, NPAS = 1000, NRUN = 4000)
! X(0:NX), Y(0:NY), Z(0:NZ) : Space mesh, V(0:NX,0:NY,0:NZ) : Potential 
! CP(0:NX,0:NY,0:NZ) : Wave function (COMPLEX)
      DOUBLE PRECISION X(0:NX), Y(0:NY), Z(0:NZ)
      DOUBLE PRECISION X2(0:NX), Y2(0:NY), Z2(0:NZ)
      DOUBLE PRECISION V(0:NX, 0:NY, 0:NZ), R2(0:NX, 0:NY, 0:NZ)
      DOUBLE PRECISION G, DX, DY, DZ, DT, MU, EN, ZNORM, RMS
      DOUBLE COMPLEX CI, CP(0:NX,0:NY,0:NZ)
      DOUBLE COMPLEX CT0X, CALA(0:NX), CGAA(0:NX)
      DOUBLE COMPLEX CT0Y, CALB(0:NY), CGAB(0:NY)
      DOUBLE COMPLEX CT0Z, CALC(0:NZ), CGAC(0:NZ)
      DOUBLE COMPLEX CA0, CA0R, CB0, CB0R, CC0, CC0R
! G0 : Final nonlinearity, AL, BL : Anisotropy coefficients
! OPTION and XOP decides which equation to be solved.
      DOUBLE PRECISION AL, BL
      DOUBLE PRECISION G0,  GSTP, XOP
      INTEGER K, OPTION, I, J
! DX, DY, DZ : Space step and DT : Time step 
      DATA DX /0.1D0/, DY /0.1D0/,  DZ /0.1D0/, DT/0.002D0/
      DATA CI /(0.0D0,1.0D0)/   ! Complex i
      DATA G0 /22.454D0/        ! Final nonlinearity
!
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
!
!      OPTION=1!Solves -psi_xx-psi_yy+V(x,y)psi+G0|psi|^2 psi=i psi_t
      OPTION=2!Solves [-psi_xx-psi_yy+V(x,y)psi]/2+G0|psi|^2 psi=i psi_t
!      OPTION=3! Solves -psi_xx-psi_yy+V(x,y)psi/4+G0|psi|^2 psi=i psi_t
      IF (OPTION.EQ.1) XOP = 1.D0
      IF (OPTION.EQ.2) XOP = 2.D0
      IF (OPTION.EQ.3) XOP = 1.D0
      AL = SQRT(2.D0)  !  Anisotropy 'LAMBDA' of the trap
      BL = 2.D0  !  Anisotropy 'KAPPA' of the trap
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
      CALL INITIALIZE(NX, NY, NZ, OPTION, X, Y, Z, X2, Y2, Z2, 
     &  R2, V, DX, DY, DZ, AL, BL, CP)
! COEF() defines the coefficients of Crank-Nicholson Scheme.
      CALL COEF(NXX, NYY, NZZ, DX, DY, DZ, DT, CI,  
     &  CT0X, CALA, CGAA, CT0Y, CALB, CGAB, CT0Z, CALC, CGAC,
     &  CA0, CA0R, CB0, CB0R, CC0, CC0R)
! NORM() calculates norm and restores normalization
      CALL NORM(NX, NY, NZ, DX, DY, DZ, ZNORM, CP)
! RAD() calculates the r.m.s radius RMS
      CALL RAD(NX, NY, NZ, DX, DY, DZ, RMS, R2, CP, X2, Y2, Z2)
! CHEM() calculates the chemical potential MU and energy EN.
      CALL CHEM(NX, NY, NZ, DX, DY, DZ, G, MU, EN, V, CP)
  999 FORMAT (2F12.6, F16.8)
      WRITE (7, 1001)
      WRITE (7, 1002)
      WRITE (7, 1001)
      WRITE (7, 1003) ZNORM,MU/XOP,EN/XOP,RMS,ABS(CP(NX2, NY2, NZ2))
 1001 FORMAT (18X, 54('-'))
 1002 FORMAT (20X, 'Norm', 6X, 'Chem', 8X, 'Ener', 8X, '<r>', 5X, 
     &    'Psi(0,0,0)')
 1003 FORMAT ('Initial :  ', 3X, F11.3, 2F12.3, 3F11.3)
 1000 FORMAT(F12.4,F16.6)
!
            DO I = 0, NX
                WRITE(1,1000)X(I),ABS(CP(I,NY/2,NZ/2))
             END DO
             DO J = 0, NY
                WRITE(2,1000)Y(J),ABS(CP(NX/2,J,NZ/2))
             END DO
             DO K = 0, NZ
                WRITE(3,1000)Z(K),ABS(CP(NX/2,NY/2,K))
             END DO



      GSTP = XOP*G0/DFLOAT(NSTP)     
      DO K = 1, NSTP ! NSTP iterations sets the Nonlinear term
         G = G + GSTP
         CALL NU(NX, NY, NZ, V, DT, G, CP, CI)  
         CALL LUX(NX, NY, NZ, NXX, CA0R, CP, CT0X, CALA, CGAA)
         CALL LUY(NX, NY, NZ, NYY, CB0R, CP, CT0Y, CALB, CGAB)
         CALL LUZ(NX, NY, NZ, NZZ, CC0R, CP, CT0Z, CALC, CGAC)
      END DO
      CALL NORM(NX, NY, NZ, DX, DY, DZ, ZNORM, CP)
      CALL CHEM(NX, NY, NZ, DX, DY, DZ, G, MU, EN, V, CP)
      CALL RAD(NX, NY, NZ, DX, DY, DZ, RMS, R2, CP, X2, Y2, Z2)
      WRITE (7, 1004) ZNORM,MU/XOP,EN/XOP,RMS,ABS(CP(NX2, NY2, NZ2))
 1004 FORMAT('After NSTP iter.:', F8.3, 2F12.3, 3F11.3)
!          
      DO K = 1, NPAS ! NPAS iterations transient
         CALL NU(NX, NY, NZ, V, DT, G, CP, CI)  
         CALL LUX(NX, NY, NZ, NXX, CA0R, CP, CT0X, CALA, CGAA)
         CALL LUY(NX, NY, NZ, NYY, CB0R, CP, CT0Y, CALB, CGAB)
         CALL LUZ(NX, NY, NZ, NZZ, CC0R, CP, CT0Z, CALC, CGAC)
      END DO
      CALL NORM(NX, NY, NZ, DX, DY, DZ, ZNORM, CP)
      CALL CHEM(NX, NY, NZ, DX, DY, DZ, G, MU, EN, V, CP)
      CALL RAD(NX, NY, NZ, DX, DY, DZ, RMS, R2, CP, X2, Y2, Z2)
      WRITE (7, 1005) ZNORM,MU/XOP,EN/XOP,RMS,ABS(CP(NX2, NY2, NZ2))
      WRITE (7, 1001)
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
         CALL NU(NX, NY, NZ, V, DT, G, CP, CI)
         CALL LUX(NX, NY, NZ, NXX, CA0R, CP, CT0X, CALA, CGAA)
         CALL LUY(NX, NY, NZ, NYY, CB0R, CP, CT0Y, CALB, CGAB)
         CALL LUZ(NX, NY, NZ, NZZ, CC0R, CP, CT0Z, CALC, CGAC)
         IF (MOD(K,4).EQ.0) THEN
            CALL RAD(NX, NY, NZ, DX, DY, DZ, RMS, R2, CP, X2, Y2, Z2)
            WRITE(8,999) DFLOAT(K)*DT*XOP, RMS 
         END IF
      END DO
      STOP
      END
!
      SUBROUTINE INITIALIZE(NX, NY, NZ, OPTION, X, Y, Z, X2, Y2, Z2, 
     &  R2, V, DX, DY, DZ, AL, BL, CP)
!   Routine that initizlizes the constant and variables.
!   Calculates the potential term V and the initial wave function CP
      IMPLICIT NONE
      INTEGER NX, NY, NZ, OPTION
      DOUBLE PRECISION X(0:NX), Y(0:NY), Z(0:NZ)
      DOUBLE PRECISION X2(0:NX), Y2(0:NY), Z2(0:NZ)
      DOUBLE PRECISION R2(0:NX, 0:NY, 0:NZ), V(0:NX, 0:NY, 0:NZ)
      DOUBLE PRECISION DX, DY, DZ, AL, BL
      DOUBLE COMPLEX CP(0:NX,0:NY,0:NZ)
      DOUBLE PRECISION PI, SPI, TMP, SP, AL2, BL2
      INTEGER I, J, K, NX2, NY2, NZ2
      PI = 4.0D0*ATAN(1.0D0)
      SP = SQRT(PI*SQRT(PI/(AL*BL)))
      SPI = SP*SQRT(2.0D0*SQRT(2.0D0))
      NX2 = NX/2
      NY2 = NY/2
      NZ2 = NZ/2
      AL2 = AL*AL
      BL2 = BL*BL
      DO I = 0, NX ! Setting up initial condition
         X(I) = DFLOAT(I-NX2)*DX
         X2(I) = X(I)*X(I)
         DO J = 0, NY
            Y(J) = DFLOAT(J-NY2)*DY
            Y2(J) = Y(J)*Y(J)
            DO K = 0, NZ
               Z(K) = DFLOAT(K-NZ2)*DZ
               Z2(K) = Z(K)*Z(K)
               R2(I,J,K) = X2(I) + Y2(J) + Z2(K)
               IF (OPTION .EQ. 3) THEN
                  V(I,J,K) = (X2(I)+AL2*Y2(J)+BL2*Z2(K))/4.0D0
                  TMP = EXP(-(X2(I)+AL*Y2(J)+BL*Z2(K))/4.0D0)
                  CP(I,J,K) = CMPLX(TMP/SPI,0.0D0)
               ELSE
                  V(I,J,K) = (X2(I)+AL2*Y2(J)+BL2*Z2(K))
                  TMP = EXP(-(X2(I)+AL*Y2(J)+BL*Z2(K))/2.0D0)
                  CP(I,J,K) = CMPLX(TMP/SP,0.0D0)
               END IF
            END DO
         END DO
      END DO
      RETURN
      END
!
      SUBROUTINE COEF(NXX, NYY, NZZ, DX, DY, DZ, DT, CI,  
     &  CT0X, CALA, CGAA, CT0Y, CALB, CGAB, CT0Z, CALC, CGAC,
     &  CA0, CA0R, CB0, CB0R, CC0, CC0R)
!   Calculates the coefficients needed in subroutines LUX, LUY, LUZ      
      IMPLICIT NONE
      INTEGER NXX, NYY, NZZ
      DOUBLE PRECISION DX, DY, DZ, DT
      DOUBLE COMPLEX CI, CT0X, CALA(0:NXX), CGAA(0:NXX)
      DOUBLE COMPLEX CT0Y, CALB(0:NYY), CGAB(0:NYY)
      DOUBLE COMPLEX CT0Z, CALC(0:NZZ), CGAC(0:NZZ)
      DOUBLE COMPLEX CA0, CA0R, CB0, CB0R, CC0, CC0R
      DOUBLE PRECISION DX2, DY2, DZ2, DXX, DYY, DZZ
      INTEGER I, J, K
      DOUBLE COMPLEX CDT
      DX2 = DX*DX       ! Generating the Coefficients for the
      DY2 = DY*DY       ! C-N method (to solve the spatial part)
      DZ2 = DZ*DZ
      DXX = 1.0D0/DX2
      DYY = 1.0D0/DY2
      DZZ = 1.0D0/DZ2
      CDT = CI*DT
      CA0 = 1.0D0+CDT*DXX
      CA0R = 1.0D0-CDT*DXX
      CB0 = 1.0D0+CDT*DYY
      CB0R = 1.0D0-CDT*DYY
      CC0 = 1.0D0+CDT*DZZ
      CC0R = 1.0D0-CDT*DZZ
      CT0X = -CDT*DXX/2.0D0
      CALA(NXX) = 0.0D0
      CGAA(NXX) = -1.0D0/CA0
      DO I = NXX, 1, -1
         CALA(I-1) = CT0X*CGAA(I)
         CGAA(I-1) = -1.0D0/(CA0+CT0X*CALA(I-1))
      END DO
      CT0Y = -CDT*DYY/2.0D0
      CALB(NYY) = 0.0D0
      CGAB(NYY) = -1.0D0/CB0
      DO J = NYY, 1, -1
         CALB(J-1) = CT0Y*CGAB(J)
         CGAB(J-1) = -1.0D0/(CB0+CT0Y*CALB(J-1))
      END DO
      CT0Z = -CDT*DZZ/2.0D0
      CALC(NZZ) = 0.0D0
      CGAC(NZZ) = -1.0D0/CC0
      DO K = NZZ, 1, -1
         CALC(K-1) = CT0Z*CGAC(K)
         CGAC(K-1) = -1.0D0/(CC0+CT0Z*CALC(K-1))
      END DO
      RETURN
      END
!
      SUBROUTINE NU(NX, NY, NZ, V, DT, G, CP, CI)  ! Exact solution
!   Solves the partial differential equation with the potential and the
!   nonlinear term.
      IMPLICIT NONE
      INTEGER NX, NY, NZ
      DOUBLE PRECISION V(0:NX,0:NY,0:NZ), DT, G
      DOUBLE COMPLEX CP(0:NX,0:NY,0:NZ), CI
      DOUBLE PRECISION P2, TMP,P3
      INTEGER I, J, K
      DO I = 0, NX
         DO K = 0, NZ
            DO J = 0, NY
               P2 = ABS(CP(I,J,K))
               P2 = P2*P2
               CALL NONLIN(P2, P3, G)
               TMP = DT*(V(I,J,K) + P3)
               CP(I,J,K) = CP(I,J,K)*EXP(-CI*TMP)
            END DO
         END DO
      END DO
      RETURN
      END
!
      SUBROUTINE LUX(NX, NY, NZ, NXX, CA0R, CP, CT0X, CALA, CGAA)
!  Solves the partial differential equation only with the X space
!  derivative term using the Crank-Nicholson method
      IMPLICIT NONE
      INTEGER NX, NY, NZ, NXX
      DOUBLE COMPLEX CA0R, CP(0:NX,0:NY,0:NZ)
      DOUBLE COMPLEX CT0X, CALA(0:NXX), CGAA(0:NXX)
      DOUBLE COMPLEX CBE(0:NXX), CXX
      INTEGER I, J, K
      DO J = 0, NY            
         DO K = 0, NZ            
            CBE(NXX) = CP(NX,J,K)
            DO I = NXX, 1, -1
               CXX = -CT0X*CP(I+1,J,K)+CA0R*CP(I,J,K)-CT0X*CP(I-1,J,K)
               CBE(I-1) = CGAA(I)*(CT0X*CBE(I)-CXX)
            END DO
            CP(0,J,K) = 0.0D0     
            DO I = 0, NXX
               CP(I+1,J,K) = CALA(I)*CP(I,J,K)+CBE(I)
            END DO
            CP(NX,J,K) = 0.0D0    ! CP(NXX,J,K)
         END DO
      END DO      
      RETURN 
      END
!
      SUBROUTINE LUY(NX, NY, NZ, NYY, CB0R, CP, CT0Y, CALB, CGAB)
!  Solves the partial differential equation only with the Y space
!  derivative term using the Crank-Nicholson method
      IMPLICIT NONE
      INTEGER NX, NY, NZ, NYY
      DOUBLE COMPLEX CB0R, CP(0:NX,0:NY,0:NZ)
      DOUBLE COMPLEX CT0Y, CALB(0:NYY), CGAB(0:NYY)
      DOUBLE COMPLEX CBE(0:NYY), CYY
      INTEGER I, J, K
      DO I = 0, NX            
         DO K = 0, NZ            
            CBE(NYY) = CP(I,NY,K)
            DO J = NYY, 1, -1
               CYY = -CT0Y*CP(I,J+1,K)+CB0R*CP(I,J,K)-CT0Y*CP(I,J-1,K)
               CBE(J-1) = CGAB(J)*(CT0Y*CBE(J)-CYY)
            END DO
            CP(I,0,K) = 0.0D0      
            DO J = 0, NYY
               CP(I,J+1,K) = CALB(J)*CP(I,J,K)+CBE(J)
            END DO
            CP(I,NY,K) = 0.D0 
         END DO
      END DO
      RETURN 
      END
!
      SUBROUTINE LUZ(NX, NY, NZ, NZZ, CC0R, CP, CT0Z, CALC, CGAC)
!  Solves the partial differential equation only with the Z space
!  derivative term using the Crank-Nicholson method
      IMPLICIT NONE
      INTEGER NX, NY, NZ, NZZ
      DOUBLE COMPLEX CC0R, CP(0:NX,0:NY,0:NZ)
      DOUBLE COMPLEX CT0Z, CALC(0:NZZ), CGAC(0:NZZ)
      DOUBLE COMPLEX CBE(0:NZZ), CZZ
      INTEGER I, J, K
       DO I = 0, NX
          DO J = 0, NY            
             CBE(NZZ) = CP(I,J,NZ)
             DO K = NZZ, 1, -1
                CZZ = -CT0Z*CP(I,J,K+1)+CC0R*CP(I,J,K)-CT0Z*CP(I,J,K-1)
                CBE(K-1) = CGAC(K)*(CT0Z*CBE(K)-CZZ)
             END DO
             CP(I,J,0) = 0.0D0      
             DO K = 0, NZZ
                CP(I,J,K+1) = CALC(K)*CP(I,J,K)+CBE(K)
             END DO
             CP(I,J,NZ) = 0.D0  ! CP(I,NYY)
          END DO
       END DO      
      RETURN 
      END
!
      SUBROUTINE NORM(NX, NY, NZ, DX, DY, DZ, ZNORM, CP)
!   Calculates the normalization of the wave function and sets it to
!   unity.
      IMPLICIT NONE
      INTEGER NX, NY, NZ
      DOUBLE PRECISION DX, DY, DZ, ZNORM
      DOUBLE COMPLEX CP(0:NX,0:NY,0:NZ)
      DOUBLE PRECISION P2, TMPX(0:NX), TMPY(0:NY), TMPZ(0:NZ)
      INTEGER I, J, K
      DOUBLE PRECISION SIMP
      EXTERNAL SIMP 
      DO I = 0, NX
         DO J = 0, NY
            DO K = 0, NZ
               P2 = ABS(CP(I,J,K))
               TMPZ(K) = P2*P2
            END DO
            TMPY(J) = SIMP(NZ, TMPZ, DZ)
         END DO
         TMPX(I) = SIMP(NY, TMPY, DY)
      END DO
      ZNORM = SQRT(SIMP(NX, TMPX, DX))
      DO I = 0, NX
         DO J = 0, NY
            DO K = 0, NZ
               CP(I,J,K) = CP(I,J,K)/ZNORM
            END DO
         END DO
      END DO
      RETURN
      END
!
      SUBROUTINE CHEM(NX, NY, NZ, DX, DY, DZ, G, MU, EN, V, CP)
!  Calculates the chemical potential MU and energy EN.  CP is the wave
!  function, V is the potential and G is the nonlinearity.
      IMPLICIT NONE
      INTEGER NX, NY, NZ
      DOUBLE PRECISION DX, DY, DZ, G, MU, EN
      DOUBLE PRECISION V(0:NX,0:NY,0:NZ)
      DOUBLE COMPLEX  CP(0:NX,0:NY,0:NZ)
      DOUBLE PRECISION XCP(0:NX),XDPX(0:NX)
      DOUBLE PRECISION YCP(0:NY),YDPY(0:NY)
      DOUBLE PRECISION ZCP(0:NZ),ZDPZ(0:NZ)
      DOUBLE PRECISION P2, DP2, TMPX(0:NX), TMPY(0:NY), TMPZ(0:NZ)
      DOUBLE PRECISION EMPX(0:NX), EMPY(0:NY), EMPZ(0:NZ)
      DOUBLE PRECISION DPX(0:NX,0:NY,0:NZ), DPY(0:NX,0:NY,0:NZ)
      DOUBLE PRECISION DPZ(0:NX,0:NY,0:NZ)
      DOUBLE PRECISION P3
      INTEGER I, J, K
      DOUBLE PRECISION SIMP
      EXTERNAL SIMP
      DO J = 0, NY
         DO K = 0, NZ
            DO I = 0, NX
               XCP(I)=ABS(CP(I,J,K))
            END DO
            CALL DIFF(NX, XCP, XDPX, DX)
            DO I = 0, NX
               DPX(I,J,K)=XDPX(I)
            END DO
         END DO
      END DO
!
      DO I = 0, NX
         DO K = 0, NZ
            DO J = 0, NY
               YCP(J)=ABS(CP(I,J,K))
            END DO
            CALL DIFF(NY,YCP, YDPY, DY)
            DO J = 0, NY
               DPY(I,J,K)=YDPY(J)
            END DO
         END DO
      END DO
!
      DO I = 0, NX
         DO J = 0, NY
            DO K = 0, NZ
               ZCP(K)=ABS(CP(I,J,K))
            END DO
            CALL DIFF(NZ, ZCP, ZDPZ, DZ)
            DO K = 0, NZ
               DPZ(I,J,K)=ZDPZ(K)
            END DO
         END DO
      END DO
!
      DO I = 0, NX
         DO J = 0, NY
            DO K = 0, NZ
               P2 = ABS(CP(I,J,K))*ABS(CP(I,J,K))
               DP2 = DPX(I,J,K)**2 + DPY(I,J,K)**2 + DPZ(I,J,K)**2
               CALL NONLIN(P2,P3,G)
               TMPZ(K) = (V(I,J,K) + P3)*P2 + DP2
               EMPZ(K) = (V(I,J,K) + P3/2.0D0)*P2 + DP2
            END DO
            TMPY(J) = SIMP(NZ, TMPZ, DZ)
            EMPY(J) = SIMP(NZ, EMPZ, DZ)
         END DO
         TMPX(I) = SIMP(NY, TMPY, DY)
         EMPX(I) = SIMP(NY, EMPY, DY)
      END DO
      MU = SIMP(NX, TMPX, DX)
      EN = SIMP(NX, EMPX, DX)
      RETURN
      END
!
      SUBROUTINE RAD(NX, NY, NZ, DX, DY, DZ, RMS, R2, CP, X2, Y2, Z2)
!  Calculates the root mean square size RMS
      IMPLICIT NONE
      INTEGER NX, NY, NZ
      DOUBLE PRECISION DX, DY, DZ, RMS
      DOUBLE PRECISION X2(0:NX),Y2(0:NY),Z2(0:NZ)
      DOUBLE PRECISION R2(0:NX,0:NY,0:NZ)
      DOUBLE COMPLEX CP(0:NX,0:NY,0:NZ)
      DOUBLE PRECISION P2, TMPX(0:NX), TMPY(0:NY), TMPZ(0:NZ)
      INTEGER I, J, K
      DOUBLE PRECISION SIMP
      EXTERNAL SIMP
      DO I = 0, NX
         DO J = 0, NY
            DO K = 0, NZ
               P2 = ABS(CP(I,J,K))**2
               TMPZ(K) = (X2(I)+Y2(J)+Z2(K))*P2 
            END DO
            TMPY(J) = SIMP(NZ, TMPZ, DZ)
         END DO
         TMPX(I) = SIMP(NY, TMPY, DY)
      END DO
      RMS = SQRT(SIMP(NX, TMPX, DX))
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
      SUBROUTINE DIFF(N, P, DP, DX)
      ! Computes the first derivative DP of P using
      ! Richardsonextrapolation formula. The derivative at the
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
!
      SUBROUTINE NONLIN(P2,P3,G)
! Calculates the nonlinear term
      IMPLICIT NONE
      DOUBLE PRECISION P2,P3,G
      P3=G*P2
      RETURN
      END
!
     
