!#  File name : real1d.f90	!# Last modified : 14 Mar 2016
!#  Fortran program for Gross-Pitaevskii equation in one-dimensional 
!#  harmonic trap by real time propagation (Fortran 90/95 Version)
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
! N : Number of space mesh points
  INTEGER, PARAMETER :: N = 3000, N2 = N/2, NX = N-1
! NSTP : Number of iterations to introduce the nonlinearity. 
! NSTP=0 reads the wave function, /=0 calculates the wave function.
! NPAS : Number of subsequent iterations with fixed nonlinearity.
! NRUN : Number of final time iterations with fixed nonlinearity.
  INTEGER, PARAMETER :: NSTP = 1000000, NPAS = 100, NRUN = 40000
  REAL (8), PARAMETER :: PI = 3.14159265358979D0
END MODULE COMM_DATA
!
 MODULE GPE_DATA
  USE COMM_DATA, ONLY : N, PI
  REAL (8), PARAMETER :: AHO = 1.D-6     		! Unit of length ( l = 1 MICRON)            
  REAL (8), PARAMETER :: Bohr_a0 =  5.2917720859D-11/AHO  	  ! Bohr radius (scaled with AHO)
  COMPLEX (8), PARAMETER :: CI = (0.D0,1.D0) 	 	! Complex i
!
  REAL (8), PARAMETER :: DX = 0.01D0, DT = 0.0001D0	! DX : Space step and DT : Time step
  INTEGER, PARAMETER  :: NATOMS = 2000			! Number of Atoms
  REAL (8), PARAMETER :: AS = 74.1032482D0*Bohr_a0	! Scattering length (in units of Bohr_a0)
  REAL (8), PARAMETER :: GAMMA = 1.D0			! Parameter of Trap
  REAL (8), PARAMETER :: DRHO = 0.5D0, DRHO2 = DRHO * DRHO	  ! DRHO = AH0/SQRT(LAMBDA*KAPPA) : Transverse trap parameter 
! 
  REAL (8), PARAMETER :: G_1D = 4.D0*PI*AS*NATOMS/(2.D0*PI*DRHO2) ! G_1D : Nonlinearity in the one-dimensional GP equation
  REAL (8), PARAMETER :: G_3D = G_1D*2.D0*PI*DRHO2	! G_3D : Three-dimensional nonlinearity 
  REAL (8), PARAMETER :: GPAR = 0.5D0 	                ! Change for dynamics
!
! OPTION   decides which equation to be solved.
! OPTION=1 Solves -psi_xx+V(x)psi+G_1D|psi|^2 psi =i psi_t
! OPTION=2 Solves [-psi_xx+V(x)psi]/2+G_1D|psi|^2 psi =i psi_t
  INTEGER, PARAMETER :: OPTION = 2 
! X(0:N):Space mesh, V(0:N):Potential, CP(0:N):Wave function  
  REAL (8), DIMENSION(0:N) :: X, X2, V
  COMPLEX (8), DIMENSION(0:N) :: CP
  REAL (8) :: G,  GSTP,   XOP
END MODULE GPE_DATA

MODULE CN_DATA
  USE COMM_DATA, ONLY : N
  COMPLEX (8), DIMENSION(0:N) :: CAL, CGA
  COMPLEX (8) :: CAIPM
END MODULE CN_DATA 

PROGRAM GROSS_PITAEVSKII_SSCN_1D
   USE COMM_DATA
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
       COMPLEX (8), DIMENSION(0:), INTENT(INOUT) :: CP
     END SUBROUTINE INITIALIZE
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
      COMPLEX (8), DIMENSION(0:), INTENT(INOUT) :: CP
      REAL (8), INTENT(IN) :: DT
    END SUBROUTINE CALCNU
  END INTERFACE
!------------------------
  INTERFACE 
    SUBROUTINE LU(CP)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:), INTENT(INOUT) :: CP
    END SUBROUTINE LU
  END INTERFACE
!------------------------
  INTERFACE
    SUBROUTINE NORM(CP, ZNORM)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:), INTENT(INOUT) :: CP
      REAL (8), INTENT(OUT) :: ZNORM
    END  SUBROUTINE NORM
  END INTERFACE
!------------------------
  INTERFACE
    SUBROUTINE RAD(CP2, RMS)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: CP2
      REAL (8), INTENT(OUT) :: RMS
    END  SUBROUTINE RAD
  END INTERFACE
!------------------------
  INTERFACE
    SUBROUTINE CHEM(CP, MU, EN)
      IMPLICIT NONE
      COMPLEX (8), DIMENSION(0:), INTENT(IN) :: CP
      REAL (8), INTENT(OUT) :: MU, EN
    END  SUBROUTINE CHEM
  END INTERFACE
!------------------------
  INTERFACE
    SUBROUTINE WRITE_DENSITY(FUNIT, U2)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FUNIT
      REAL (8), DIMENSION(0:), INTENT(IN) :: U2
    END SUBROUTINE WRITE_DENSITY
  END INTERFACE
!------------------------ INTERFACE BLOCKS -----------------------
  INTEGER :: K
  REAL (8) :: ZNORM, MU, EN, RMS, T, T1, T2
  REAL (8), DIMENSION(0:N) :: CP2
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
! CALCULATE_TRAP() initializes the harmonic oscillator potential V.
  CALL CALCULATE_TRAP()
! COEF() defines the coefficients of Crank-Nicholson Scheme.
  CALL COEF()
!
  OPEN(17, FILE = 'real1d-out.txt')
!
  WRITE(17,900) OPTION 
  WRITE(17,*) 
  WRITE(17,901) NATOMS, AHO
  WRITE(17,902) AS/Bohr_a0 
  WRITE(17,910) G_3D 
  WRITE(17,903) G_1D 
  WRITE(17,904) GAMMA
  WRITE(17,*)
  WRITE(17,905) N
  WRITE(17,906) DX
  WRITE(17,907) NSTP, NPAS, NRUN  
  WRITE(17,908) DT 
  WRITE(17,*)
  WRITE(17, 1012) GPAR 
  WRITE(17,*)
  WRITE(17, 1001)
  WRITE(17, 1002)
  WRITE(17, 1001)
  900 FORMAT(' Real time propagation 1d,   OPTION = ',I3 )
  901 FORMAT('  Number of Atoms N =',I10,', Unit of length AHO =',F12.8,' m')
  902 FORMAT('  Scattering length a = ',F9.2,'*a0')
  910 FORMAT('  Nonlinearity G_3D =',F16.7 )
  903 FORMAT('  Nonlinearity G_1D =',F16.7 )
  904 FORMAT('  Parameter of trap: GAMMA =',F7.2)
  905 FORMAT(' # Space Stp:  N = ', I8)
  906 FORMAT('  Space Step: DX = ', F10.6)
  907 FORMAT(' # Time Stp : NSTP = ',I9,', NPAS = ',I9,', NRUN = ',I9)
  908 FORMAT('   Time Step:   DT = ', F10.6)
  1012 FORMAT(' * Change for dynamics: GPAR = ',F11.3, ' *')
  1001 FORMAT (19X,'-----------------------------------------------------')
  1002 FORMAT (20X, 'Norm', 6X, 'Chem', 8X, 'Ener/N', 6X, '<x>', 6X, '|Psi(0)|^2') 
!
!   OPEN(1, FILE = 'real1d-den-init.txt')
!    CALL WRITE_DENSITY(1, CP2)
!   CLOSE(1)
!
  IF (NSTP /= 0) THEN
    GSTP  = XOP*G_1D/DFLOAT(NSTP)
    G = 0.D0
    CALL NORM(CP,ZNORM)		! NORM() calculates norm and restores normalization.
    CALL CHEM(CP, MU, EN)	! CHEM() calculates the chemical potential MU and energy EN.
    CP2 = ABS(CP)**2
    CALL RAD(CP2, RMS)		! RAD() calculates the r.m.s radius RMS.
    WRITE (17, 1003) ZNORM, MU/XOP, EN/XOP, RMS, CP2(N2)
! 
    DO K = 1, NSTP		! NSTP iterations to introduce the nonlinearity
      G = G + GSTP
      CALL CALCNU(CP, DT)
      CALL LU(CP)
    END DO
    CALL NORM(CP, ZNORM)
    CALL CHEM(CP, MU, EN)
    CP2 = ABS(CP)**2
    CALL RAD(CP2, RMS)  
    WRITE (17, 1004) ZNORM, MU/XOP, EN/XOP, RMS, CP2(N2)
    1004 FORMAT('After NSTP iter.:', F11.4, F12.4, F12.4, 2F11.3)
!
!   OPEN(2, FILE = 'real1d-den-nstp.txt')
!    CALL WRITE_DENSITY(2, CP2)
!   CLOSE(2)
!
  ELSE
    G = XOP*G_1D
    CALL NORM(CP,ZNORM)
    CALL CHEM(CP, MU, EN)
    CP2 = ABS(CP)**2
    CALL RAD(CP2, RMS)
    WRITE (17, 1003) ZNORM, MU/XOP, EN/XOP, RMS, CP2(N2)
  END IF
!
  T = 0.D0  
  IF((NPAS /= 0).OR.(NRUN /= 0)) OPEN(8, FILE = 'real1d-dyna.txt')
!
  DO K = 1, NPAS ! NPAS iterations transient
     T = T + DT
     CALL CALCNU(CP, DT)
     CALL LU(CP)
      IF (MOD(K, 50) == 0) THEN
        CP2 = ABS(CP)**2
        CALL RAD(CP2, RMS)
        WRITE(8,1006) T*XOP, RMS
      END IF
  END DO
  IF(NPAS /= 0)THEN  
    CALL NORM(CP, ZNORM)
    CALL CHEM(CP, MU, EN)
    CP2 = ABS(CP)**2
    CALL RAD(CP2, RMS)
    WRITE (17, 1005) ZNORM, MU/XOP, EN/XOP, RMS, CP2(N2)
    1005 FORMAT('After NPAS iter.:',F11.4, F12.4, F12.4, 2F11.3)
  END IF
!
!  OPEN(1, FILE = 'real1d-den.txt')
!   CALL WRITE_DENSITY(1, CP2)
!  CLOSE(1)
!
!  The following line defines a problem which is studied in the time evolution.
  G = GPAR*G
!  
  DO K = 1, NRUN ! NRUN iterations to study a nonstationary problem
     T = T + DT
     CALL CALCNU(CP, DT)
     CALL LU(CP)
      IF (MOD(K, 50) == 0) THEN
        CP2 = ABS(CP)**2
        CALL RAD(CP2, RMS)
        WRITE(8,1006) T*XOP, RMS
      END IF
  END DO
  CLOSE (8)
  IF(NRUN /= 0)THEN 
    CALL NORM(CP, ZNORM)
    CALL CHEM(CP, MU, EN)
    CP2 = ABS(CP)**2
    CALL RAD(CP2, RMS)
    WRITE (17, 1007) ZNORM, MU/XOP, EN/XOP, RMS, CP2(N2)
    1007 FORMAT('After NRUN iter.:',F11.4, F12.4, F12.4, 2F11.3)
!
   OPEN(1, FILE = 'real1d-den.txt')
    CALL WRITE_DENSITY(1, CP2)
   CLOSE(1)
!
  END IF
! 
  CALL SYSTEM_CLOCK (CLCK_COUNTS_END, CLCK_RATE)
  CALL CPU_TIME(T2)
  WRITE (17, 1001)
  WRITE (17,*)
  WRITE (17,'(A,I7,A)') ' Clock Time: ', (CLCK_COUNTS_END - CLCK_COUNTS_BEG)/INT (CLCK_RATE,8), ' seconds'
  WRITE (17,'(A,I7,A)') '   CPU Time: ', INT(T2-T1), ' seconds' 
  CLOSE (17)
  1003 FORMAT ('Initial :    ', 4X, F11.4, F12.4, F12.4, 2F11.3)
  1006 FORMAT(2F12.4)
!
END PROGRAM GROSS_PITAEVSKII_SSCN_1D

SUBROUTINE INITIALIZE(CP)
! Routine that initializes the constant and variables and calculates the initial wave function CP.
  USE COMM_DATA, ONLY : N, N2, PI, NSTP
  USE GPE_DATA, ONLY : DX, X, X2, GAMMA 
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:), INTENT(INOUT) :: CP
  REAL (8) ::  TX, TMP, PI4
  INTEGER :: I
! 
  FORALL (I=0:N) X(I) = (I-N2)*DX
  X2 = X*X
  PI4 = SQRT(SQRT(PI/GAMMA))
!         
  IF (NSTP == 0) THEN
    WRITE(*,'(a)') "Run the program using the input file to read. e.g.: ./real1d < imag1d-den.txt"
     DO I = 0,N
        READ (*,*)TX, TMP	  !READ OUTPUT of IMAGTIME PROGRAM WITH SAME PARAMETERS 
        CP(I) = SQRT(TMP)
     END DO
  ELSE 
   CP = EXP(-X2*GAMMA/2.D0)/PI4
  END IF
END SUBROUTINE INITIALIZE

SUBROUTINE CALCULATE_TRAP()
! Routine that  initializes the harmonic oscillator potential V.
  USE GPE_DATA, ONLY : XOP, X2, V, GAMMA
  IMPLICIT NONE
!
  V = X2*GAMMA**2
END SUBROUTINE CALCULATE_TRAP

SUBROUTINE COEF()
!   Calculates the coefficients needed in subroutine LU.  
  USE COMM_DATA, ONLY :  NX
  USE GPE_DATA, ONLY : DT, DX, CI
  USE CN_DATA
  IMPLICIT NONE
  INTEGER :: J
  REAL (8) :: DX2
  COMPLEX (8) :: CAI0, CDT
!
  CDT = DT/CI       ! Generating the Coefficients for the
  DX2 = DX*DX       ! C-N method (to solve the spatial part)
  CAIPM = CDT/2.D0/DX2
  CAI0 = 1.D0-CDT/DX2
  CAL(NX) = 0.D0
  CGA(NX) = -1.D0/CAI0
  DO J = NX, 1, -1
     CAL(J-1) = CGA(J)*CAIPM
     CGA(J-1) = -1.D0/(CAI0+CAIPM*CAL(J-1))
  END DO
END SUBROUTINE COEF

SUBROUTINE CALCNU(CP, DT) ! Exact solution
! Solves the partial differential equation with the potential and the 
! nonlinear term.
   USE COMM_DATA, ONLY : N, NX
   USE GPE_DATA, ONLY : X2, V, G,  CI 
   IMPLICIT NONE
!-------------------------------------------------
   COMPLEX (8), DIMENSION(0:), INTENT(INOUT) :: CP
   REAL (8), INTENT(IN) :: DT
   REAL (8), DIMENSION(0:N) :: P, P2, TMP1D, DP
!
   P = ABS(CP)
   P2 = P*P
! 
  TMP1D = DT*(V + G*P2)
  CP = CP*EXP(-CI*TMP1D)
END SUBROUTINE CALCNU

 SUBROUTINE LU(CP)
!  Solves the partial differential equation only with the space 
!  derivative term using the Crank-Nicholson method
   USE COMM_DATA, ONLY : N, NX
   USE CN_DATA
   IMPLICIT NONE
   COMPLEX (8), DIMENSION(0:), INTENT(INOUT) :: CP
   COMPLEX (8), DIMENSION(0:NX) :: CBE
   COMPLEX (8) :: CXX
   INTEGER :: I
!   
   CBE(NX) = CP(N)
   DO I = NX, 1, -1
      CXX = CP(I)-(CP(I+1)-2.D0*CP(I)+CP(I-1))*CAIPM
      CBE(I-1) = CGA(I)*(CAIPM*CBE(I)-CXX)
   END DO
!-----------------------------
! Boundary condition reflecting:
   CP(0) = (0.D0, 0.D0)
   DO I = 0, NX-1
      CP(I+1) = CAL(I)*CP(I) + CBE(I)
   END DO
   CP(N) = (0.D0, 0.D0)
!-----------------------------
END SUBROUTINE LU

SUBROUTINE NORM(CP, ZNORM)
! Calculates the normalization of the wave function and sets it to unity.
  USE COMM_DATA, ONLY : N
  USE GPE_DATA, ONLY : DX
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(OUT) :: ZNORM
  !----------------------------------------------------------------------
  INTERFACE
    PURE FUNCTION SIMP(F, DX) RESULT (VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN)  :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
  !----------------------------------------------------------------------
  REAL (8), DIMENSION(0:N) :: P, P2
!
  P = ABS(CP)
  P2 = P*P
  ZNORM = SQRT(SIMP(P2, DX))
!  CP = CP/ZNORM
END SUBROUTINE NORM
 
SUBROUTINE RAD(CP2, RMS)
! Calculates the root mean square size RMS.
  USE COMM_DATA, ONLY : N
  USE GPE_DATA, ONLY : DX, X2
  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(IN) :: CP2
  REAL (8), INTENT(OUT) :: RMS
!----------------------------------------------------------------------
   INTERFACE
    FUNCTION SIMP(F, DX) RESULT (VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
!----------------------------------------------------------------------
 REAL (8), DIMENSION(0:SIZE(CP2)-1) :: P2
!
  P2 = X2 * CP2
  RMS = SQRT(SIMP(P2, DX))
END SUBROUTINE RAD

SUBROUTINE CHEM(CP, MU, EN)
!  Calculates the chemical potential MU and energy EN.  CP is the wave
!  function, V is the potential and G is the nonlinearity. 
  USE COMM_DATA, ONLY : N,NX  
  USE GPE_DATA, ONLY : DX, V, G 
  IMPLICIT NONE
  COMPLEX (8), DIMENSION(0:), INTENT(IN) :: CP
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
    PURE FUNCTION SIMP(F, DX) RESULT(VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
!----------------------------------------------------------------------
  REAL (8), DIMENSION(0:SIZE(CP)-1) :: TMP1D, P2, DP2, GP2
  REAL (8), DIMENSION(0:SIZE(CP)-1) :: DCP, DP,  EMP1D
!
  TMP1D = ABS(CP)
  DCP = DIFF(TMP1D, DX)
  DP2 = DCP*DCP
  P2 = TMP1D * TMP1D
  GP2 = G * P2
!  
  TMP1D = (V + GP2)*P2 + DP2
  EMP1D = (V + GP2/2.D0)*P2 + DP2
!  
  MU = SIMP(TMP1D, DX)
  EN = SIMP(EMP1D, DX)
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
!
  N = SIZE(F) - 1
  F1 = F(1) + F(N-1) ! N EVEN 
  F2 = F(2) 
  DO I = 3, N-3, 2
     F1 = F1 + F(I)
     F2 = F2 + F(I+1)
  END DO
  VALUE = DX*(F(0) + 4.D0*F1 + 2.D0*F2 + F(N))/3.D0
END FUNCTION SIMP

PURE FUNCTION DIFF(P,DX)  RESULT (DP)
! Computes the first derivative DP of P using
! Richardson extrapolation formula. The derivative at the  
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
  USE COMM_DATA, ONLY : N
  USE GPE_DATA, ONLY : X
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: FUNIT
  REAL (8), DIMENSION(0:), INTENT(IN) :: U2
  INTEGER :: I
!
  DO I = 0, N
    WRITE(FUNIT, 999) X(I), U2(I)
  END DO
  999 FORMAT(F12.4, E17.5E3)
END SUBROUTINE WRITE_DENSITY
 
!# File name : real1d.f90
