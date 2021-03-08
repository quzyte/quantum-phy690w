!#  File name : real2d.f90	!# Last modified : 14 Mar 2016
!#  Fortran program for Gross-Pitaevskii equation in two-dimensional 
!#  anisotropic trap by real time propagation (Fortran 90/95 Version)
!
!# BEC-GP-OMP codes are developed and (c)opyright-ed by:
!#
!# Luis E. Young-S., Sadhan K. Adhikari
!# (UNESP - Sao Paulo State University, Brazil)
!#
!# Paulsamy Muruganandam
!# (Bharathidasan University, Tamil Nadu, India)
!#
!# Dusan Vudragovic, Antun Balaz
!# (Scientific Computing Laboratory, Institute of Physics Belgrade, Serbia)
!#
!# Public use and modification of this code are allowed provided that the 
!# following three papers are cited:
!#
!# [1] L. E. Young-S. et al., Comput. Phys. Commun. 204 (2016) 209. 
!# [2] P. Muruganandam, S. K. Adhikari, Comput. Phys. Commun. 180 (2009) 1888.
!# [3] D. Vudragovic et al., Comput. Phys. Commun. 183 (2012) 2021.
!#
!# The authors would be grateful for all information and/or comments
!# regarding the use of the code.
!
!# To compile :
!# (1) Intel Fortran Compiler
! ifort -O3 -parallel -par-report -V -mcmodel medium -shared-intel 
!
!# (2) GNU Fortran (gfortran)
! gfortran -O3 -fopenmp
!
MODULE COMM_DATA
! NX, NY : Number of space mesh points (X and Y)  
  INTEGER, PARAMETER :: NX = 200, NXX = NX-1, NX2 = NX/2 ! No. of X steps
  INTEGER, PARAMETER :: NY = 200, NYY = NY-1, NY2 = NY/2 ! No. of Y steps
! NSTP : Number of iterations to introduce the nonlinearity. 
! NSTP=0 reads the wave function, /=0 calculates the wave function.
! NPAS : Number of subsequent iterations with fixed nonlinearity.
! NRUN : Number of final time iterations with fixed nonlinearity.
  INTEGER, PARAMETER :: NSTP = 100000, NPAS = 100, NRUN = 5000
  REAL (8), PARAMETER :: PI = 3.14159265358979D0, SQR_2PI = 2.506628274631D0
  REAL (8), PARAMETER :: TWO_PI = 2.D0*PI, FOU_PI3 = 4.D0*PI/3.D0
  INTEGER, PARAMETER :: NN = NX*NY
END MODULE COMM_DATA

MODULE GPE_DATA
  USE COMM_DATA, ONLY : NX, NY, PI, SQR_2PI
  REAL (8), PARAMETER :: AHO = 1.D-6     		  ! Unit of length (l = 1 MICRON)            
  REAL (8), PARAMETER :: Bohr_a0 =  5.2917720859D-11/AHO  ! Bohr radius (scaled with AHO)
  COMPLEX (8), PARAMETER :: CI = (0.D0,1.D0) 	 	  ! Complex i
!  
  REAL (8), PARAMETER :: DX = 0.1D0, DY = 0.1D0, DT = 0.001D0	! DX, DY : SPACE STEPS and DT : TIME STEP 
  INTEGER, PARAMETER  :: NATOMS = 500			  ! Number of Atoms
  REAL (8), PARAMETER :: AS = 94.60134*Bohr_a0	 	  ! Scattering length (in units of Bohr_a0)
  REAL (8), PARAMETER :: GAMMA = 1.D0, NU = 1.D0	  ! GAMMA and NU : Parameteres of Trap
  REAL (8), PARAMETER :: D_Z = 1.D0			  ! D_Z = Axial Gaussian Width = l/SQRT(Lambda)
!
  REAL (8), PARAMETER :: G_2D = 4.D0*PI*AS*NATOMS/(SQR_2PI*D_Z)	! G_2D : Nonlinearity in the two-dimensional GP equation
  REAL (8), PARAMETER :: G_3D = G_2D*SQR_2PI*D_Z		! G_3D : Three-dimensional nonlinearity 
  REAL (8), PARAMETER :: GPAR = 0.5D0 				! Change for dynamics
!
! OPTION  decides which equation to be solved.
! OPTION=1 Solves -psi_xx-psi_yy+V(x,y)psi+G_2D|psi|^2 psi =i psi_t
! OPTION=2 Solves [-psi_xx-psi_yy+V(x,y)psi]/2+G_2D|psi|^2 psi =i psi_t
  INTEGER, PARAMETER :: OPTION = 2 
! X(0:NX), Y(0:NY): Space mesh, V(0:NX,0:NY) : Potential, CP(0:NX,0:NY): Wave function 
  REAL (8), DIMENSION(0:NX, 0:NY) :: V, R2
  COMPLEX (8), DIMENSION(0:NX, 0:NY) :: CP
  REAL (8), DIMENSION(0:NX) :: X, X2
  REAL (8), DIMENSION(0:NY) :: Y, Y2
  REAL (8), DIMENSION(NX, NY) :: VDD
  REAL (8) :: G, GSTP, XOP, DZ2, GAM2, ANU2 
END MODULE GPE_DATA

MODULE CN_DATA
  USE COMM_DATA, ONLY : NX, NY
  COMPLEX (8), DIMENSION(0:NX) :: CALA, CGAA
  COMPLEX (8), DIMENSION(0:NY) :: CALB, CGAB
  COMPLEX (8) :: CT0X, CT0Y
  COMPLEX (8) :: CA0, CB0, CA0R, CB0R
END MODULE CN_DATA 

PROGRAM GROSS_PITAEVSKII_SSCN_2D
   USE COMM_DATA, ONLY : NX, NY, NX2, NY2, NSTP, NPAS, NRUN
   USE GPE_DATA
   IMPLICIT NONE
! Subroutine INITIALIZE() used to initialize the space mesh X(I) and the initial wave function. 
! Subroutine CALCULATE_TRAP() used to initialize the harmonic oscillator potential V. 
! Subroutine COEF() used to generate the coefficients for the Crank-Nicholson Scheme. 
! The routine CALCNU() performs time progation for the non-derivative part and LU() performs
! time propagation of derivative part. NORM() calculates the norm and 
! normalizes the wave function, CHEM() and RAD() are used to calculate the 
! chemical potential, energy and the rms radius, respectively. The function 
! DIFF() used to calculate the space derivatives of the wave function used 
! in CHEM() and SIMP() does the integration by Simpson's rule.
!------------------------ INTERFACE BLOCKS -----------------------
   INTERFACE 
     SUBROUTINE INITIALIZE(CP)
       IMPLICIT NONE
       COMPLEX (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
     END SUBROUTINE INITIALIZE
   END INTERFACE
!------------------------
   INTERFACE 
     SUBROUTINE CALCULATE_TRAP()
     END SUBROUTINE CALCULATE_TRAP
   END INTERFACE
!------------------------
   INTERFACE 
     SUBROUTINE COEF()
       IMPLICIT NONE
     END SUBROUTINE COEF
   END INTERFACE
!------------------------
   INTERFACE 
     SUBROUTINE CALCNU(CP, DT)
       IMPLICIT NONE
       COMPLEX (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
       REAL (8), INTENT(IN) :: DT
     END SUBROUTINE CALCNU
   END INTERFACE
!------------------------
   INTERFACE 
     SUBROUTINE LUX(CP)
       IMPLICIT NONE
       COMPLEX (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
     END SUBROUTINE LUX
   END INTERFACE
!------------------------
   INTERFACE 
     SUBROUTINE LUY(CP)
       IMPLICIT NONE
       COMPLEX (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
     END SUBROUTINE LUY
   END INTERFACE
!------------------------
   INTERFACE
     SUBROUTINE NORM(CP, ZNORM)
       IMPLICIT NONE
       COMPLEX (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
       REAL (8), INTENT(OUT) :: ZNORM
     END  SUBROUTINE NORM
   END INTERFACE
!------------------------
   INTERFACE
     SUBROUTINE RAD(CP2, RMS)
       IMPLICIT NONE
       REAL (8), DIMENSION(0:,0:), INTENT(IN) :: CP2
       REAL (8), DIMENSION(:), INTENT(OUT) :: RMS
     END  SUBROUTINE RAD
   END INTERFACE
!------------------------
   INTERFACE
     SUBROUTINE CHEM(CP, MU, EN)
       IMPLICIT NONE
       COMPLEX (8), DIMENSION(0:,0:), INTENT(IN) :: CP
       REAL (8), INTENT(OUT) :: MU, EN
     END  SUBROUTINE CHEM
   END INTERFACE
!------------------------
  INTERFACE
    SUBROUTINE WRITE_DENSITY(FUNIT, U2)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FUNIT
      REAL (8), DIMENSION(0:, 0:), INTENT(IN) :: U2
    END SUBROUTINE WRITE_DENSITY
  END INTERFACE
!------------------------
  INTERFACE
    SUBROUTINE DENSITY_1DX(FUNIT, U2, DE1DX)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FUNIT
      REAL (8), DIMENSION(0:,0:), INTENT(IN) :: U2
      REAL (8), DIMENSION(0:), INTENT(OUT) :: DE1DX
    END  SUBROUTINE DENSITY_1DX
  END INTERFACE
!------------------------
  INTERFACE
    SUBROUTINE DENSITY_1DY(FUNIT, U2, DE1DY)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FUNIT
      REAL (8), DIMENSION(0:,0:), INTENT(IN) :: U2
      REAL (8), DIMENSION(0:), INTENT(OUT) :: DE1DY
    END  SUBROUTINE DENSITY_1DY
  END INTERFACE
!------------------------ END INTERFACE BLOCKS -------------------
  INTEGER :: K
  REAL (8) :: ZNORM, MU, EN, T, T1, T2
  REAL (8), DIMENSION(0:NX, 0:NY) :: CP2
!  <x> = RMS(1), <y> = RMS(2), <rho> = RMS(3)
  REAL (8), DIMENSION(3) :: RMS
  REAL (8), DIMENSION(0:NX) :: DEN1X
  REAL (8), DIMENSION(0:NY) :: DEN1Y
  INTEGER :: CLCK_COUNTS_BEG, CLCK_COUNTS_END, CLCK_RATE
!
  CALL SYSTEM_CLOCK ( CLCK_COUNTS_BEG, CLCK_RATE )
  CALL CPU_TIME(T1)
! 
  SELECT CASE (OPTION)
    CASE (1)
      XOP = 1.D0
    CASE (2)
      XOP = 2.D0
    CASE (:0,3:)
      PRINT *, 'ERROR: Wrong OPTION', OPTION
      STOP
  END SELECT
!
! INITIALIZE() initializes the starting normalized wave function 'CP'.
  CALL INITIALIZE(CP)
!
  OPEN(17, FILE = 'real2d-out.txt')
  OPEN(14, FILE = 'real2d-rms.txt')  
!
  WRITE(17,900) OPTION 
  WRITE(14,900) OPTION 
  WRITE(17,*) 
  WRITE(14,*) 
  WRITE(17,901) NATOMS, AHO
  WRITE(17,902) AS/Bohr_a0 
  WRITE(17,910) G_3D 
  WRITE(17,903) G_2D 
  WRITE(17,904) GAMMA, NU
  WRITE(17,909) D_Z
  WRITE(17,*)
  WRITE(17,905) NX, NY
  WRITE(17,906) DX, DY 
  WRITE(17,907) NSTP, NPAS, NRUN  
  WRITE(17,908) DT 
  WRITE(17,*)
  WRITE(17, 1014) GPAR 
  WRITE(17,*)
  900 FORMAT(' Real time propagation 2d,   OPTION = ',I3 )
  901 FORMAT('  Number of Atoms N =',I10,', Unit of length AHO =',F12.8,' m')
  902 FORMAT('  Scattering length a = ',F9.2,'*a0')
  910 FORMAT('  Nonlinearity G_3D =',F16.6 )
  903 FORMAT('  Nonlinearity G_2D =',F16.6 )
  904 FORMAT('  Parameters of trap: GAMMA =',F7.2, ', NU =',F7.2)
  909 FORMAT('  Axial Gaussian Width =',F7.2)
  905 FORMAT(' # Space Stp: NX = ', I8, ', NY = ', I8)
  906 FORMAT('  Space Step: DX = ', F10.6, ', DY = ', F10.6)
  907 FORMAT(' # Time Stp : NSTP = ',I9,', NPAS = ',I9,', NRUN = ',I9)
  908 FORMAT('   Time Step:   DT = ', F10.6)
! CALCULATE_TRAP() initializes the harmonic oscillator potential V.
  CALL CALCULATE_TRAP()
! COEF() defines the coefficients of Crank-Nicholson Scheme.
  CALL COEF()
!
  WRITE (17, 1001)
  WRITE (17, 1002)
  WRITE (17, 1001)
  WRITE (14, 1001)
  WRITE (14, 1011)
  WRITE (14, 1001)
  1001 FORMAT (19X,'-----------------------------------------------------')
  1002 FORMAT (20X, 'Norm', 6X, 'Chem', 8X, 'Ener/N', 5X, '<rho>', 5X, '|Psi(0,0)|^2')
  1011 FORMAT ('Values of rms size:', 8X, '<x>', 13X, '<y>', 13X, '<rho>')
  1003 FORMAT ('Initial    : ', 4X,  F8.4, F12.3, F12.4, 2F11.3)
  1013 FORMAT (11X, 'Initial:', F14.5, 2F16.5)
!  
   CP2 = ABS(CP)**2
!   OPEN(21, FILE = 'real2d-den-init.txt')
!   CALL WRITE_DENSITY(21, CP2)
!   CLOSE(21)
! 
!   OPEN(11, FILE = 'real2d-den-init1d_x.txt')
!   OPEN(12, FILE = 'real2d-den-init1d_y.txt')
!   CALL DENSITY_1DX(11, CP2, DEN1X)
!   CALL DENSITY_1DY(12, CP2, DEN1Y)
!   CLOSE(11)
!   CLOSE(12)
!
  IF(NSTP /= 0)THEN
     GSTP = XOP*G_2D/DFLOAT(NSTP)
     G = 0.D0
     CALL NORM(CP,ZNORM)	! NORM() calculates norm and restores normalization.
     CALL CHEM(CP, MU, EN)	! CHEM() calculates the chemical potential MU and energy EN.
     CP2 = ABS(CP)**2
     CALL RAD(CP2, RMS)		! RAD() calculates the r.m.s radius RMS.
     WRITE (17, 1003) ZNORM, MU/XOP, EN/XOP, RMS(3), CP2(NX2, NY2)
     WRITE (14, 1013) RMS(1:3)
! 
     DO K = 1, NSTP 		! NSTP iterations to introduce the nonlinearity
       G = G + GSTP
       CALL CALCNU(CP, DT)
       CALL LUX(CP)
       CALL LUY(CP)
     END DO
     CALL NORM(CP,ZNORM)
     CALL CHEM(CP, MU, EN)
     CP2=ABS(CP)**2
     CALL RAD(CP2, RMS)
     WRITE (17, 1005) ZNORM, MU/XOP, EN/XOP, RMS(3), CP2(NX2, NY2)
     WRITE (14, 1015) RMS(1:3)
     1005 FORMAT('After NSTP iter.:',  F8.4, F12.3, F12.4, 2F11.3)
     1015 FORMAT (2X, 'After NSTP iter.:', F14.5, 2F16.5) 
!  
!   OPEN(21, FILE = 'real2d-den-nstp.txt')
!   CALL WRITE_DENSITY(21, CP2)
!   CLOSE(21)
! 
!   OPEN(11, FILE = 'real2d-den-nstp1d_x.txt')
!   OPEN(12, FILE = 'real2d-den-nstp1d_y.txt')
!   CALL DENSITY_1DX(11, CP2, DEN1X)
!   CALL DENSITY_1DY(12, CP2, DEN1Y)
!   CLOSE(11)
!   CLOSE(12)
!
  ELSE
    G = XOP*G_2D
    CALL NORM(CP, ZNORM)
    CALL CHEM(CP, MU, EN) 
    CP2 = ABS(CP)**2
    CALL RAD(CP2, RMS)  
    WRITE (17, 1003) ZNORM, MU/XOP, EN/XOP, RMS(3), CP2(NX2, NY2)
    WRITE (14, 1013) RMS(1:3) 
  END IF
!
  T=0.D0
  IF((NPAS /= 0).OR.(NRUN /= 0)) OPEN(8, FILE = 'real2d-dyna.txt')
! 
  DO K = 1, NPAS ! NPAS iterations transient
     T = T + DT
     CALL CALCNU(CP, DT)
     CALL LUX(CP)
     CALL LUY(CP)
      IF (MOD(K,20).EQ.0) THEN
       CP2=ABS(CP)**2
       CALL RAD(CP2, RMS)
       WRITE(8,1008) T*XOP, RMS(1), RMS(2), RMS(3)
      END IF
  END DO
  IF(NPAS /= 0)THEN 
    CALL NORM(CP,ZNORM)
    CALL CHEM(CP, MU, EN)
    CP2 = ABS(CP)**2
    CALL RAD(CP2, RMS)
    WRITE (17, 1006) ZNORM, MU/XOP, EN/XOP, RMS(3), CP2(NX2, NY2)
    WRITE (14, 1016) RMS(1:3)
    1006 FORMAT('After NPAS iter.:', F8.4, F12.3, F12.4, 2F11.3)
    1016 FORMAT (2X, 'After NPAS iter.:', F14.5, 2F16.5)
  END IF
!  
  OPEN(21, FILE = 'real2d-den.txt')
    CALL WRITE_DENSITY(21, CP2)
  CLOSE(21)
  OPEN(11, FILE = 'real2d-den1d_x.txt')
  OPEN(12, FILE = 'real2d-den1d_y.txt')
    CALL DENSITY_1DX(11, CP2, DEN1X)
    CALL DENSITY_1DY(12, CP2, DEN1Y)
  CLOSE(11)
  CLOSE(12)
!
!  The following line defines a problem which is studied in the time evolution.
  G = GPAR*G
!
  DO K = 1, NRUN	! NRUN iterations to study a nonstationary problem
    T = T + DT 
    CALL CALCNU(CP, DT)
    CALL LUX(CP)
    CALL LUY(CP)
     IF (MOD(K,20).EQ.0) THEN
      CP2 = ABS(CP)**2
      CALL RAD(CP2, RMS)
      WRITE(8,1008) T*XOP, RMS(1), RMS(2), RMS(3)
     END IF
  END DO
  CLOSE(8)
  IF(NRUN /= 0)THEN 
    CALL NORM(CP,ZNORM)
    CALL CHEM(CP, MU, EN)
    CP2=ABS(CP)**2
    CALL RAD(CP2, RMS)
    WRITE (17, 1007) ZNORM, MU/XOP, EN/XOP, RMS(3), CP2(NX2, NY2)
    WRITE (14, 1017) RMS(1:3)
    1007 FORMAT('After NRUN iter.:', F8.4, F12.3, F12.4, 2F11.3)
    1017 FORMAT (2X, 'After NRUN iter.:', F14.5, 2F16.5)
    1014 FORMAT(' * Change for dynamics: GPAR = ',F11.3,' *')
!  
!   OPEN(21, FILE = 'real2d-den-nrun.txt')
!   CALL WRITE_DENSITY(21, CP2)
!   CLOSE(21)
! 
!   OPEN(11, FILE = 'real2d-den-nrun1d_x.txt')
!   OPEN(12, FILE = 'real2d-den-nrun1d_y.txt')
!   CALL DENSITY_1DX(11, CP2, DEN1X)
!   CALL DENSITY_1DY(12, CP2, DEN1Y)
!   CLOSE(11)
!   CLOSE(12)
! 
  END IF
!
  CALL SYSTEM_CLOCK (CLCK_COUNTS_END, CLCK_RATE)
  CALL CPU_TIME(T2)
  WRITE (17, 1001)
  WRITE (14, 1001)
  WRITE (17,*)
  WRITE (17,'(A,I7,A)') ' Clock Time: ', (CLCK_COUNTS_END - CLCK_COUNTS_BEG)/INT (CLCK_RATE,8), ' seconds'
  WRITE (17,'(A,I7,A)') '   CPU Time: ', INT(T2-T1), ' seconds' 
  1008 FORMAT(4F10.4)
! 
END PROGRAM GROSS_PITAEVSKII_SSCN_2D

SUBROUTINE INITIALIZE(CP)
! Routine that initializes the constant and variables and calculates the initial wave function CP.
  USE COMM_DATA, ONLY : NX, NY, NX2, NY2, PI, NSTP
  USE GPE_DATA, ONLY : GAMMA, DX, DY, X, X2, Y, Y2, R2, NU, GAM2, ANU2
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
  REAL (8), DIMENSION(0:NX,0:NY):: TMP2D
  REAL (8) ::  QR_GAM, TX, TY, TMP
  INTEGER :: I, J
!  
  GAM2 = GAMMA*GAMMA
  ANU2 = NU*NU
  QR_GAM = SQRT(SQRT(NU*GAMMA)/PI)
!
  FORALL (I=0:NX) X(I) = (I-NX2)*DX
  FORALL (I=0:NY) Y(I) = (I-NY2)*DY
  X2 = X*X
  Y2 = Y*Y
  FORALL (I=0:NX) R2(I,:) = X2(I) + Y2
!         
  IF (NSTP == 0) THEN
    WRITE(*,'(a)') "Run the program using the input file to read. e.g.: ./real2d < imag2d-den.txt"
    DO I = 0,NX
      DO J = 0,NY 
        READ (*, *)TX, TY, TMP
        CP(I,J) = SQRT(TMP)
      END DO
    END DO
  ELSE 
    FORALL (I=0:NX) TMP2D(I,:) = (GAMMA*X2(I) + NU*Y2)/2.D0
    CP = QR_GAM*EXP(-TMP2D)
  END IF
END SUBROUTINE INITIALIZE

SUBROUTINE CALCULATE_TRAP()
! Routine that  initializes the harmonic oscillator potential V.
  USE COMM_DATA, ONLY : NX
  USE GPE_DATA, ONLY : XOP, GAM2, X2, Y2, V, ANU2
  IMPLICIT NONE
  INTEGER :: I
!
  DO I = 0, NX
    V(I,:) = (GAM2*X2(I) + ANU2*Y2)
  END DO
END SUBROUTINE CALCULATE_TRAP

SUBROUTINE COEF()
!   Calculates the coefficients needed in subroutine LUX and LUY
  USE COMM_DATA, ONLY : NXX, NYY
  USE GPE_DATA, ONLY : DX, DY, DT, CI
  USE CN_DATA
  IMPLICIT NONE
  INTEGER :: I, J
  REAL (8) :: DX2, DY2, DXX, DYY
  COMPLEX (8) :: CDT
!
  DX2 = DX*DX   ! Generating the Coefficients for the
  DY2 = DY*DY   ! C-N method (to solve the spatial part)           
  DXX = 1.D0/DX2
  DYY = 1.D0/DY2
  CDT = CI*DT
  CA0 = 1.D0+CDT*DXX
  CA0R = 1.D0-CDT*DXX
  CB0 = 1.D0+CDT*DYY
  CB0R = 1.D0-CDT*DYY
!
  CT0X = -CDT*DXX/2.D0
  CALA(NXX) = 0.D0
  CGAA(NXX) = -1.D0/CA0
  DO I = NXX, 1, -1
     CALA(I-1) = CT0X*CGAA(I)
     CGAA(I-1) = -1.D0/(CA0+CT0X*CALA(I-1))
  END DO
!
  CT0Y = -CDT*DYY/2.D0
  CALB(NYY) = 0.D0
  CGAB(NYY) = -1.D0/CB0
  DO J = NYY, 1, -1
     CALB(J-1) = CT0Y*CGAB(J)
     CGAB(J-1) = -1.D0/(CB0+CT0Y*CALB(J-1))
  END DO
END SUBROUTINE COEF

SUBROUTINE CALCNU(CP, DT) ! Exact solution
! Solves the partial differential equation with the potential and 
! the nonlinear term.
  USE COMM_DATA, ONLY : NX, NY, NXX, NYY
  USE GPE_DATA, ONLY : CI, V, G 
  IMPLICIT NONE
!----------------------------------------------------------------------
  COMPLEX (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(IN) :: DT
  REAL (8), DIMENSION(0:NX,0:NY) :: P, P2, TMP2D, DP
!
 
  P = ABS(CP)
  P2 = P*P
!   
  TMP2D = DT*(V + G*P2) 
  CP = CP*EXP(-CI*TMP2D)
 
END SUBROUTINE CALCNU

SUBROUTINE LUX(CP)
  ! Solves the partial differential equation only with the X-space
  ! derivative term using the Crank-Nicholson method
  USE COMM_DATA, ONLY : NX, NY, NXX
  USE CN_DATA, ONLY : CA0R, CT0X, CALA, CGAA  
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
  INTEGER :: I, J
  COMPLEX (8) :: CBE(0:NXX,0:NY), CXX(0:NXX,0:NY)
!
 DO J =0,NY 
  CBE(NXX,J) = CP(NXX,J)
  DO I = NXX, 1, -1
     CXX(I,J) = -CT0X*CP(I+1,J)+CA0R*CP(I,J)-CT0X*CP(I-1,J)
  END DO
 END DO

 DO J =0,NY 
  DO I = NXX, 1, -1
     CBE(I-1,J) = CGAA(I)*(CT0X*CBE(I,J)-CXX(I,J))
  END DO
!-----------------------------
! Boundary condition reflecting 
  CP(0,J) = 0.D0
  DO I = 0, NXX
     CP(I+1,J) = CALA(I)*CP(I,J) + CBE(I,J)
  END DO
  CP(NX,J) = 0.D0
 END DO
!-----------------------------
END SUBROUTINE LUX

SUBROUTINE LUY(CP)
  ! Solves the partial differential equation only with the Y-space
  ! derivative term using the Crank-Nicholson method
  USE COMM_DATA, ONLY : NX, NY, NYY
  USE CN_DATA, ONLY :   CB0R, CT0Y, CALB, CGAB   
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
  INTEGER :: I, J
  COMPLEX (8) :: CBE(0:NX,0:NYY), CYY(0:NX,0:NYY)
!
 DO I=0,NX
  CBE(I,NYY) = CP(I,NYY)
  DO J = NYY, 1, -1
     CYY(I,J) = -CT0Y*CP(I,J+1)+CB0R*CP(I,J)-CT0Y*CP(I,J-1)
  END DO
 END DO

 DO I=0,NX
  DO J = NYY, 1, -1
     CBE(I,J-1) = CGAB(J)*(CT0Y*CBE(I,J)-CYY(I,J))
  END DO
!-----------------------------
! Boundary condition reflecting: 
  CP(I,0) = 0.D0
  DO J = 0, NYY
     CP(I,J+1) = CALB(J)*CP(I,J) + CBE(I,J)
  END DO
  CP(I,NY) = 0.D0
 END DO
!-----------------------------
END SUBROUTINE LUY

SUBROUTINE NORM(CP, ZNORM)
! Calculates the normalization of the wave function and sets it to unity.
  USE GPE_DATA, ONLY : DX, DY
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(OUT) :: ZNORM
!----------------------------------------------------------------------
  INTERFACE
    FUNCTION INTEGRATE(P2, DX, DY) RESULT(RES)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:, 0:), INTENT(IN) :: P2
      REAL (8), INTENT(IN) :: DX, DY
      REAL (8) :: RES
    END FUNCTION INTEGRATE
  END INTERFACE
!----------------------------------------------------------------------  
  REAL (8), DIMENSION(0:SIZE(CP,1)-1,0:SIZE(CP,2)-1) :: P, TMP2D
!  
  P = ABS(CP)
  TMP2D = P*P
  ZNORM = INTEGRATE(TMP2D, DX, DY)
  !CP = CP/SQRT(ZNORM)
END SUBROUTINE NORM

SUBROUTINE RAD(P, RMS) 
! Calculates the root mean square size RMS
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : DX, DY, R2, X2, Y2
  IMPLICIT NONE
  REAL (8), DIMENSION(0:,0:), INTENT(IN) :: P
  REAL (8), DIMENSION(:), INTENT(OUT) :: RMS
!----------------------------------------------------------------------
  INTERFACE
    FUNCTION INTEGRATE(P2, DX, DY) RESULT(RES)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:, 0:), INTENT(IN) :: P2
      REAL (8), INTENT(IN) :: DX, DY
      REAL (8) :: RES
    END FUNCTION INTEGRATE
  END INTERFACE
!----------------------------------------------------------------------
  INTEGER :: I,J
  REAL (8), DIMENSION(0:SIZE(P,1)-1,0:SIZE(P,2)-1) :: TMP2D
!  
  FORALL(J=0:NY) TMP2D(:,J) = X2*P(:,J)
    RMS(1) = SQRT(INTEGRATE(TMP2D, DX, DY))
  FORALL(I=0:NX) TMP2D(I,:) = Y2*P(I,:)
    RMS(2) = SQRT(INTEGRATE(TMP2D, DX, DY))
  TMP2D = R2*P
    RMS(3) = SQRT(INTEGRATE(TMP2D, DX, DY))
END SUBROUTINE RAD

SUBROUTINE CHEM(CP, MU, EN)
!  Calculates the chemical potential MU and energy EN.  CP is the wave
!  function, V is the potential and G is the nonlinearity.
  USE COMM_DATA, ONLY : NX, NY,NXX,NYY
  USE GPE_DATA, ONLY : DX, DY, V, G 
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:,0:), INTENT(IN) :: CP
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
    FUNCTION INTEGRATE(P2, DX, DY) RESULT(RES)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:, 0:), INTENT(IN) :: P2
      REAL (8), INTENT(IN) :: DX, DY
      REAL (8) :: RES
    END FUNCTION INTEGRATE
  END INTERFACE
 !-------------------------------------------------
  INTEGER :: I
  REAL (8), DIMENSION(0:NX,0:NY) :: P, DP, DPX, DPY, P2, GP2 
  REAL (8), DIMENSION(0:NX,0:NY) :: DP2, TMP2D, EMP2D 
!
  P = ABS(CP)
   DO I = 0, NX
    DPY(I,:) = DIFF(P(I,:), DY)
   END DO
   DO I = 0, NY
    DPX(:,I) = DIFF(P(:,I), DX)
   END DO
!
  P2 = P * P
  GP2 = G*P2
  DP2 = DPX*DPX + DPY*DPY  
! 
  TMP2D = (V + GP2)*P2 + DP2      
  EMP2D = (V + GP2/2.D0)*P2 + DP2    
!     
  MU = INTEGRATE(TMP2D, DX, DY)
  EN = INTEGRATE(EMP2D, DX, DY)
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
!
  N = SIZE(F) - 1
  F1 = F(1) + F(N-1) ! N EVEN
  F2 = F(2) 
  DO I = 3, N-3, 2
     F1 = F1 + F(I)
     F2 = F2 + F(I+1)
  END DO
  SIMP = DX*(F(0) + 4.D0*F1 + 2.D0*F2 + F(N))/3.D0
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
!
  N = SIZE(P) - 1
  DP(0) = 0.D0
  DP(1) = (P(2) - P(0))/(2.D0*DX)
  FORALL(I=2:N-2)
    DP(I) = (P(I-2)-8.D0*P(I-1)+8.D0*P(I+1)-P(I+2))/(12.D0*DX)
  END FORALL
  DP(N-1) = (P(N) - P(N-2))/(2.D0*DX)
  DP(N) = 0.D0
END FUNCTION DIFF
 
SUBROUTINE WRITE_DENSITY(FUNIT, U2)
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : X, Y
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: FUNIT
  REAL (8), DIMENSION(0:, 0:), INTENT(IN) :: U2
  INTEGER :: I, J
!
  DO I = 0, NX
     DO J = 0, NY
        WRITE(FUNIT, 999) X(I), Y(J), U2(I,J)
     END DO
     WRITE(FUNIT, *)
  END DO
  999 FORMAT(2F12.4, E17.5E3)
END SUBROUTINE WRITE_DENSITY

SUBROUTINE DENSITY_1DX(FUNIT, U2, DE1DX)
  USE COMM_DATA, ONLY : NX
  USE GPE_DATA, ONLY : X, DY
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: FUNIT
  REAL (8), DIMENSION(0:, 0:), INTENT(IN) :: U2
  REAL (8), DIMENSION(0:), INTENT(OUT) :: DE1DX
!-------------------------------------------------
  INTERFACE
    PURE FUNCTION SIMP(F, DX) RESULT (VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
!-------------------------------------------------
  INTEGER :: I
!
      DO I = 0, NX 
	 DE1DX(I) = SIMP(U2(I,:), DY)
         WRITE(FUNIT, 1001) X(I), DE1DX(I)
      END DO
  1001 FORMAT(F10.2, E17.5E3)
END SUBROUTINE DENSITY_1DX

SUBROUTINE DENSITY_1DY(FUNIT, U2, DE1DY)
  USE COMM_DATA, ONLY : NY
  USE GPE_DATA, ONLY : Y, DX
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: FUNIT
  REAL (8), DIMENSION(0:, 0:), INTENT(IN) :: U2
  REAL (8), DIMENSION(0:), INTENT(OUT) :: DE1DY
!-------------------------------------------------
  INTERFACE
    PURE FUNCTION SIMP(F, DX) RESULT (VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
!-------------------------------------------------
  INTEGER :: J
!
      DO J = 0, NY 
	 DE1DY(J) = SIMP(U2(:,J), DX)
         WRITE(FUNIT, 1001) Y(J), DE1DY(J)
      END DO
  1001 FORMAT(F10.2, E17.5E3)
END SUBROUTINE DENSITY_1DY

FUNCTION INTEGRATE(U, DX, DY) RESULT(RES)
  IMPLICIT NONE
  REAL (8), DIMENSION(0:, 0:), INTENT(IN) :: U
  REAL (8), INTENT (IN) :: DX, DY
  REAL (8) :: RES
!-------------------------------------------------
  INTERFACE
    PURE FUNCTION SIMP(F, DX)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: SIMP
    END FUNCTION SIMP
  END INTERFACE
!-------------------------------------------------  
  REAL (8), DIMENSION(0:SIZE(U,1)-1) :: TMP1D
  INTEGER :: I, NX
!
  NX = SIZE(U,1)-1 
  FORALL (I = 0:NX) TMP1D(I) = SIMP(U(I,0:), DY)
  RES = SIMP(TMP1D, DX)
END FUNCTION INTEGRATE

!# File name : real2d.f90
