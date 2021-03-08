!  Filename : imagtime2d.F
!  Title : Fortran program for time-dependent Gross-Pitaevskii equation in 
!          a fully anisotropic trap
!  Authors : P. Muruganandam and S. K. Adhikari
! 
!  Fortran program for Gross-Pitaevskii equation in two-dimensional 
!  anisotropic trap by imaginary time propagation (FORTRAN 77 Version)
!
      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER NX, NXX, NX2, NY, NYY, NY2
      INTEGER NPAS, NRUN
! NX, NY : Number of space mesh points (X and Y)  
      PARAMETER (NX = 800, NXX = NX-1, NX2 = NX/2)
      PARAMETER (NY = 800, NYY = NY-1, NY2 = NY/2)
! NPAS : Number of subsequent iterations with fixed nonlinearity, and 
! NRUN : Number of final time iterations with fixed nonlinearity
       PARAMETER (NPAS = 30000, NRUN = 5000)
! X(0:N), Y(0:N) : Space mesh, V(0:N) : Potential, and 
! CP(0:N) : Wave function 
      DOUBLE PRECISION X(0:NX), Y(0:NY)
      DOUBLE PRECISION X2(0:NX), Y2(0:NY)
      DOUBLE PRECISION V(0:NX, 0:NY)
      DOUBLE PRECISION G, DX, DY, DT
      DOUBLE PRECISION CP(0:NX,0:NY)
      DOUBLE PRECISION CT0X, CALA(0:NX), CGAA(0:NX)
      DOUBLE PRECISION CT0Y, CALB(0:NY), CGAB(0:NY)
      DOUBLE PRECISION CA0, CA0R, CB0, CB0R
! G0 : Final nonlinearity, AL : Anisotrophy coefficient
! OPTION and XOP decides which equation to be solved. 
      DOUBLE PRECISION G0, ZNORM, RMS, MU, EN, AL, XOP
      INTEGER I, J, K, OPTION      
! DX, DY : Space step and DT : Time step
      DATA DX /0.02D0/, DY /0.02D0/,  DT/0.0001D0/
      DATA G0 /12.5484D0/ ! Final nonlinearity
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
      !OPTION = 1 ! Solves -psi_xx-psi_yy+V(x,y)psi+G0|psi|^2 psi=i psi_t
      OPTION = 2 !Solves [-psi_xx-psi_yy+V(x,y)psi]/2+G0|psi|^2 psi=i psi_t
      !OPTION = 3 ! Solves -psi_xx-psi_yy+V(x,y)psi/4+G0|psi|^2 psi=i psi_t
      IF(OPTION.EQ.1) XOP = 1.D0  
      IF(OPTION.EQ.2) XOP = 2.D0  
      IF(OPTION.EQ.3) XOP = 1.D0  
      G = 0.0D0
      AL = 2.D0
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
      CALL INITIALIZE(NX, NY, OPTION, X, Y, X2, Y2, V, DX, DY, AL, CP)
! COEF() defines the coefficients of Crank-Nicholson Scheme.
      CALL COEF(NXX, NYY, DX, DY, DT, CT0X, CALA, CGAA, CT0Y, 
     & CALB, CGAB, CA0, CA0R, CB0, CB0R)
! NORM() calculates norm and restores normalization
      CALL NORM(NX, NY, DX, DY, ZNORM, CP)
! RAD() calculates the r.m.s radius RMS
      CALL RAD(NX, NY, X2, Y2, DX, DY, RMS, CP)
! CHEM() calculates the chemical potential MU and energy EN.
      CALL CHEM(NX, NY, DX, DY, G, MU, EN, V, CP)
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
 1002 FORMAT (20X, 'Norm', 7X, 'Chem', 8X, 'Ener', 7X, '<r>', 7X, 
     &    'Psi(0,0)')
 1003 FORMAT ('Initial : ', 4X, F11.4, 2F12.6, 3F11.5)
      G = XOP*G0
!
      DO K = 1, NPAS ! NPAS iterations transient
! NU() performs time propagation with harmonic potential and nonlinear
! term (non-derivative parts)
         CALL NU(NX, NY, V, DT, G, CP)  
! LUX() and LUY() perform the time iteration with space derivative in
! X and Y, respectively, using Crank-Nicholson scheme.
         CALL LUX(NX, NY, NXX, CA0R, CP, CT0X, CALA, CGAA)
         CALL NORM(NX, NY, DX, DY, ZNORM, CP)
         CALL LUY(NX, NY, NYY, CB0R, CP, CT0Y, CALB, CGAB)
         CALL NORM(NX, NY, DX, DY, ZNORM, CP)
      END DO
!
      CALL RAD(NX, NY, X2, Y2, DX, DY, RMS, CP)
      CALL CHEM(NX, NY, DX, DY, G, MU, EN, V, CP)          
      WRITE (7, 1005) ZNORM, MU/XOP, EN/XOP, RMS, CP(NX2, NY2)
 1005 FORMAT('After NPAS iter.:', F8.4, 2F12.6, 3F11.5)
      DO I = 0, NX, 10 ! Writes intermediate wave funtion in File 2
          DO J = 0, NY, 10
             WRITE(2,999) X(I), Y(J), CP(I,J)
          END DO
          WRITE(2,*)
      END DO
      DO K = 1, NRUN ! NRUN iterations to check convergence
         CALL NU(NX, NY, V, DT, G, CP)  
         CALL LUX(NX, NY, NXX, CA0R, CP, CT0X, CALA, CGAA)
         CALL NORM(NX, NY, DX, DY, ZNORM, CP)
         CALL LUY(NX, NY, NYY, CB0R, CP, CT0Y, CALB, CGAB)
         CALL NORM(NX, NY, DX, DY, ZNORM, CP)
      END DO
!
      CALL RAD(NX, NY,X2,Y2, DX, DY, RMS, CP)
      CALL CHEM(NX, NY, DX, DY, G, MU, EN, V, CP)
      WRITE (7, 1006) ZNORM, MU/XOP, EN/XOP, RMS, CP(NX2, NY2)
      WRITE (7, 1001)
 1006 FORMAT('After NRUN iter.:',F8.4, 2F12.6, 3F11.5)
      DO I = 0, NX, 10 ! Writes final wave funtion in File 3
          DO J = 0, NY, 10
             WRITE(3,999) X(I), Y(J), CP(I,J)
          END DO
          WRITE(3,*)
      END DO
      STOP
      END
!
      SUBROUTINE INITIALIZE(NX,NY,OPTION,X,Y,X2,Y2,V,DX,DY,AL,CP)
!   Routine that initizlizes the constant and variables.
!   Calculates the potential term V and the initial wave function CP
      IMPLICIT NONE
      INTEGER NX, NY,OPTION
      DOUBLE PRECISION X(0:NX), Y(0:NY)
      DOUBLE PRECISION X2(0:NX), Y2(0:NY)
      DOUBLE PRECISION V(0:NX, 0:NY)
      DOUBLE PRECISION DX, DY,AL
      DOUBLE PRECISION CP(0:NX,0:NY)
      DOUBLE PRECISION PI, SPI, TMP,SP
      INTEGER I, J, NX2, NY2
      PI = 4.0D0*ATAN(1.0D0)
      SPI = SQRT(2.D0*PI)
      SP=SQRT(PI)
      NX2 = NX/2
      NY2 = NY/2
      DO I = 0, NX      ! Setting up initial condition
         X(I) = DFLOAT(I-NX2)*DX
         X2(I) = X(I)*X(I)
         DO J = 0, NY
            Y(J) = DFLOAT(J-NY2)*DY
            Y2(J) = Y(J)*Y(J)
           IF(OPTION.EQ.3)THEN
            V(I,J) = (X2(I)+AL*AL*Y2(J))/4.0D0  
            TMP = EXP(-(X2(I)+AL*Y2(J))/4.D0)
            CP(I,J) = TMP*AL**0.25D0/SPI       
           ELSE
            V(I,J) = (X2(I)+AL*AL*Y2(J))  
            TMP = EXP(-(X2(I)+AL*Y2(J))/2.D0)
            CP(I,J) = TMP*AL**0.25D0/SP        
           END IF 
         END DO
      END DO
      RETURN
      END
! 
      SUBROUTINE COEF(NXX, NYY, DX, DY, DT, CT0X, CALA, CGAA, 
     &  CT0Y, CALB, CGAB, CA0, CA0R, CB0, CB0R)
!   Calculates the coefficients needed in subroutine LUX and LUY.
      IMPLICIT NONE
      INTEGER NXX, NYY
      DOUBLE PRECISION DX, DY, DT
      DOUBLE PRECISION CT0X, CALA(0:NXX), CGAA(0:NXX)
      DOUBLE PRECISION CT0Y, CALB(0:NYY), CGAB(0:NYY)
      DOUBLE PRECISION CA0, CA0R, CB0, CB0R
      DOUBLE PRECISION DX2, DY2, DXX, DYY
      INTEGER I, J           
      DX2 = DX*DX       ! Generating the Coefficients for the
      DY2 = DY*DY       ! C-N method (to solve the spatial part)      
      DXX = 1.0D0/DX2
      DYY = 1.0D0/DY2
      CA0 = 1.0D0+DT*DXX
      CA0R = 1.0D0-DT*DXX
      CB0 = 1.0D0+DT*DYY
      CB0R = 1.0D0-DT*DYY 
      CT0X = -DT*DXX/2.0D0
      CALA(NXX) = 0.0D0
      CGAA(NXX) = -1.0D0/CA0
      DO I = NXX, 1, -1
         CALA(I-1) = CT0X*CGAA(I)
         CGAA(I-1) = -1.0D0/(CA0+CT0X*CALA(I-1))
      END DO
      CT0Y = -DT*DYY/2.0D0
      CALB(NYY) = 0.0D0
      CGAB(NYY) = -1.0D0/CB0
      DO J = NYY, 1, -1
         CALB(J-1) = CT0Y*CGAB(J)
         CGAB(J-1) = -1.0D0/(CB0+CT0Y*CALB(J-1))
      END DO
      RETURN
      END
!
      SUBROUTINE NU(NX, NY, V, DT, G, CP)  ! Exact solution
!   Solves the partial differential equation with the potential and the
!   nonlinear term.
      IMPLICIT NONE
      INTEGER NX, NY
      DOUBLE PRECISION V(0:NX,0:NY), DT, G
      DOUBLE PRECISION CP(0:NX,0:NY)
      DOUBLE PRECISION P2, P3,TMP
      INTEGER I, J            
      DO I = 0, NX
         DO J = 0, NY
            P2 = CP(I,J)*CP(I,J)
            CALL NONLIN(P2,P3,G)
            TMP = DT*(V(I,J) + P3)
            CP(I,J) = CP(I,J)*EXP(-TMP)
         END DO
      END DO
      RETURN
      END
!
      SUBROUTINE LUX(NX, NY, NXX, CA0R, CP, CT0X, CALA, CGAA)
!  Solves the partial differential equation only with the X-space
!  derivative term using the Crank-Nicholson method
      IMPLICIT NONE
      INTEGER NX, NY, NXX
      DOUBLE PRECISION CA0R, CP(0:NX,0:NY)
      DOUBLE PRECISION CT0X, CALA(0:NXX), CGAA(0:NXX)
      DOUBLE PRECISION CBE(0:NXX), CXX
      INTEGER I, J
      DO J = 0, NY            
         CBE(NXX) = CP(NX,J)
         DO I = NXX, 1, -1
            CXX = -CT0X*CP(I+1,J)+CA0R*CP(I,J)-CT0X*CP(I-1,J)
            CBE(I-1) = CGAA(I)*(CT0X*CBE(I)-CXX)
         END DO
         CP(0,J) = 0.0D0     
         DO I = 0, NXX
            CP(I+1,J) = CALA(I)*CP(I,J)+CBE(I)
         END DO
         CP(NX,J) = 0.0D0    ! CP(NXX,J)
      END DO      
      RETURN 
      END
!
      SUBROUTINE LUY(NX, NY, NYY, CB0R, CP, CT0Y, CALB, CGAB)
!  Solves the partial differential equation only with the Y-space
!  derivative term using the Crank-Nicholson method
      IMPLICIT NONE
      INTEGER NX, NY, NYY
      DOUBLE PRECISION CB0R, CP(0:NX,0:NY)
      DOUBLE PRECISION CT0Y, CALB(0:NYY), CGAB(0:NYY)
      DOUBLE PRECISION CBE(0:NYY), CYY
      INTEGER I, J
      DO I = 0, NX            
         CBE(NYY) = CP(I,NY)
         DO J = NYY, 1, -1
            CYY = -CT0Y*CP(I,J+1)+CB0R*CP(I,J)-CT0Y*CP(I,J-1)
            CBE(J-1) = CGAB(J)*(CT0Y*CBE(J)-CYY)
         END DO
         CP(I,0) = 0.0D0      
         DO J = 0, NYY
            CP(I,J+1) = CALB(J)*CP(I,J)+CBE(J)
         END DO
         CP(I,0) = 0.D0 
         CP(I,NY) = 0.D0 
      END DO
      RETURN 
      END
!
      SUBROUTINE NORM(NX, NY, DX, DY, ZNORM, CP)
!  Calculates the normalization of the wave function and sets it to
!  unity.
      IMPLICIT NONE
      INTEGER NX, NY
      DOUBLE PRECISION DX, DY, ZNORM
      DOUBLE PRECISION CP(0:NX,0:NY)
      DOUBLE PRECISION  TMPX(0:NX), TMPY(0:NY)
      INTEGER I, J
      DOUBLE PRECISION SIMP
      EXTERNAL SIMP 
      DO I = 0, NX
         DO J = 0, NY
            TMPY(J) = CP(I,J)*CP(I,J)
         END DO
         TMPX(I) = SIMP(NY, TMPY, DY)
      END DO
      ZNORM = SQRT(SIMP(NX, TMPX, DX))
      DO I = 0, NX
         DO J = 0, NY
               CP(I,J) = CP(I,J)/ZNORM
         END DO
      END DO
      RETURN
      END
!
      SUBROUTINE RAD(NX, NY, X2, Y2, DX, DY, RMS, CP)
!  Calculates the root mean square size RMS
      IMPLICIT NONE
      INTEGER NX, NY
      DOUBLE PRECISION DX, DY, RMS
      DOUBLE PRECISION CP(0:NX,0:NY),X2(0:NX),Y2(0:NY)
      DOUBLE PRECISION P2, TMPX(0:NX), TMPY(0:NY)
      INTEGER I, J
      DOUBLE PRECISION SIMP
      EXTERNAL SIMP 
      DO I = 0, NX
         DO J = 0, NY
            P2 = CP(I,J)* CP(I,J)
            TMPY(J) = P2*(X2(I)+Y2(J))
         END DO
         TMPX(I) = SIMP(NY, TMPY, DY)
      END DO
      RMS = SQRT(SIMP(NX, TMPX, DX))
      RETURN
      END
!
      SUBROUTINE CHEM(NX, NY, DX, DY, G, MU, EN, V, CP)
!  Calculates the chemical potential MU and energy EN.  CP is the wave
!  function, V is the potential and G is the nonlinearity.
      IMPLICIT NONE
      INTEGER NX, NY
      DOUBLE PRECISION DX, DY, G, MU, EN, V(0:NX,0:NY)
      DOUBLE PRECISION CP(0:NX,0:NY)
      DOUBLE PRECISION XCP(0:NX),YCP(0:NY)
      DOUBLE PRECISION P2, P3, DP2, TMPX(0:NX),TMPY(0:NY)
      DOUBLE PRECISION EMPX(0:NX),EMPY(0:NY)
      DOUBLE PRECISION DPX(0:NX,0:NY), DPY(0:NX,0:NY)
      DOUBLE PRECISION XDPX(0:NX), YDPY(0:NY)
      INTEGER I, J
      DOUBLE PRECISION SIMP
      EXTERNAL SIMP 
      DO J = 0, NY
      DO I = 0, NX
          XCP(I)=CP(I,J)
      END DO
         CALL DIFF(NX,XCP, XDPX, DX)
      DO I = 0, NX
         DPX(I,J)=XDPX(I)
      END DO
      END DO
      DO I = 0, NX
      DO J = 0, NY
          YCP(J)=CP(I,J)
         END DO
         CALL DIFF(NY,YCP, YDPY, DY)
      DO J = 0, NY
         DPY(I,J)=YDPY(J)
      END DO
      END DO
      DO I = 0, NX
         DO J = 0, NY
               P2 = CP(I,J)*CP(I,J)
               DP2 = DPX(I,J)*DPX(I,J) + DPY(I,J)*DPY(I,J)
               CALL NONLIN(P2,P3,G)
               TMPY(J) = (V(I,J) + P3)*P2 + DP2
               EMPY(J) = (V(I,J) + P3/2.0D0)*P2 + DP2
         END DO
         TMPX(I) = SIMP(NY, TMPY, DY)
         EMPX(I) = SIMP(NY, EMPY, DY)
      END DO
      MU = SIMP(NX, TMPX, DX)
      EN = SIMP(NX, EMPX, DX)
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
      DP(0) = 0.0D0
      DP(1) = (P(2) - P(0))/(2.0D0*DX)
      DO I = 2, N-2
        DP(I) = (P(I-2)-8.0D0*P(I-1)+8.0D0*P(I+1)-P(I+2))/(12.0D0*DX)
      END DO
      DP(N-1) = (P(N) - P(N-2))/(2.0D0*DX)
      DP(N) = 0.0D0
      RETURN      
      END     

