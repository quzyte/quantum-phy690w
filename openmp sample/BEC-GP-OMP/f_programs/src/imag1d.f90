!#  File name : imag1d.f90	!# Last modified : 14 Mar 2016
!#  Fortran program for Gross-Pitaevskii equation in one-dimensional 
!#  harmonic trap by imaginary time propagation (Fortran 90/95 Version).
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
  INTEGER, PARAMETER :: N = 8000, N2 = N/2, NX = N-1
! NSTP : Number of iterations to introduce the nonlinearity. 
! NPAS : Number of subsequent iterations with fixed nonlinearity.
! NRUN : Number of final time iterations with fixed nonlinearity.
  INTEGER, PARAMETER :: NSTP = 0, NPAS = 200000, NRUN = 20000
  REAL (8), PARAMETER :: PI = 3.14159265358979D0
END MODULE COMM_DATA

MODULE GPE_DATA
  USE COMM_DATA, ONLY : N, PI
  REAL (8), PARAMETER :: AHO = 1.D-6     		! Unit of length (l = 1 MICRON)            
  REAL (8), PARAMETER :: Bohr_a0 =  5.2917720859D-11/AHO	  ! Bohr radius (scaled with AHO)
!
  REAL (8), PARAMETER :: DX = 0.0025D0, DT = 0.00002D0	! DX : Space step and DT : Time step
  INTEGER, PARAMETER  :: NATOMS = 2000			! Number of ATOMS
  REAL (8), PARAMETER :: AS = 74.103248D0*Bohr_a0	! Scattering length (in units of Bohr_a0)
  REAL (8), PARAMETER :: GAMMA = 1.D0			! Parameter of trap
  REAL (8), PARAMETER :: DRHO = 0.5D0, DRHO2 = DRHO * DRHO	  ! DRHO = AH0/SQRT(KAPPA*LAMBDA) : Transverse trap parameter 
!
  REAL (8), PARAMETER :: G_1D = 4.D0*PI*AS*NATOMS/(2.D0*PI*DRHO2) ! G_1D : Nonlinearity in the one-dimensional GP equation
  REAL (8), PARAMETER :: G_3D = G_1D*2.D0*PI*DRHO2	! G_3D : Three-dimensional nonlinearity 
!
! OPTION   decides which equation to be solved.
! OPTION=1 Solves -psi_xx+V(x)psi+G_1D|psi|^2 psi=i psi_t
! OPTION=2 Solves [-psi_xx+V(x)psi]/2+G_1D|psi|^2 psi=i psi_t
  INTEGER, PARAMETER :: OPTION = 2
! X(0:N):Space mesh, V(0:N):Potential, CP(0:N):Wave function 
  REAL (8), DIMENSION(0:N) :: X, X2, V, CP
  REAL (8) :: XOP, G
END MODULE GPE_DATA

MODULE CN_DATA
  USE COMM_DATA, ONLY : N
  REAL (8), DIMENSION(0:N) :: CAL, CGA
  REAL (8) :: CAIPM
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
!------------------------ Interface Blocks -----------------------
 INTERFACE 
    SUBROUTINE INITIALIZE()
      IMPLICIT NONE
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
      REAL (8), DIMENSION(0:), INTENT(INOUT) :: CP
      REAL (8), INTENT(IN) :: DT
    END SUBROUTINE CALCNU
  END INTERFACE
!------------------------
  INTERFACE 
    SUBROUTINE LU(CP)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(INOUT) :: CP
    END SUBROUTINE LU
  END INTERFACE
!------------------------
  INTERFACE
    SUBROUTINE NORM(CP, ZNORM)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(INOUT) :: CP
      REAL (8), INTENT(OUT) :: ZNORM
    END  SUBROUTINE NORM
  END INTERFACE
!------------------------
  INTERFACE
    SUBROUTINE RAD(CP2, RMS)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(INOUT) :: CP2
      REAL (8), INTENT(OUT) :: RMS
    END  SUBROUTINE RAD
  END INTERFACE
!------------------------
  INTERFACE
    SUBROUTINE CHEM(CP, MU, EN)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: CP
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
  REAL (8) :: ZNORM, MU, RMS, EN, T1, T2
  REAL (8), DIMENSION(0:N) :: CP2
  REAL (8) :: GSTP  
  INTEGER :: CLCK_COUNTS_BEG, CLCK_COUNTS_END, CLCK_RATE
!
  CALL SYSTEM_CLOCK ( CLCK_COUNTS_BEG, CLCK_RATE )
  CALL CPU_TIME (T1)
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
  CALL INITIALIZE()
!
  OPEN(7, FILE = 'imag1d-out.txt')
  WRITE(7,900) OPTION 
  WRITE(7,*) 
  WRITE(7,901) NATOMS, AHO
  WRITE(7,902) AS/Bohr_a0 
  WRITE(7,910) G_3D
  WRITE(7,903) G_1D
  WRITE(7,904) GAMMA
  WRITE(7,909) DRHO
  WRITE(7,*)
  WRITE(7,905) N
  WRITE(7,906) DX
  WRITE(7,907) NSTP, NPAS, NRUN  
  WRITE(7,908) DT 
  WRITE(7,*)
  900 FORMAT(' Imaginary time propagation 1d,   OPTION = ',I3 )
  901 FORMAT('  Number of Atoms N =',I10,', Unit of length AHO =',F12.8,' m')
  902 FORMAT('  Scattering length a = ',F9.2,'*a0')
  910 FORMAT('  Nonlinearity G_3D =',F16.7)
  903 FORMAT('  Nonlinearity G_1D =',F16.7)
  904 FORMAT('  Parameter of trap: GAMMA =',F7.2)
  909 FORMAT('  Transverse trap parameter =',F7.3)
  905 FORMAT(' # Space Stp:  N = ', I8)
  906 FORMAT('  Space Step: DX = ', F10.6)
  907 FORMAT(' # Time Stp : NSTP = ',I9,', NPAS = ',I9,', NRUN = ',I9)
  908 FORMAT('   Time Step:   DT = ', F10.6)
!
! CALCULATE_TRAP() initializes the harmonic oscillator potential V.
  CALL CALCULATE_TRAP()
! COEF() defines the coefficients of Crank-Nicholson Scheme.
  CALL COEF()
! NORM() calculates and restore normalization.
  CALL NORM(CP, ZNORM)
! CHEM() calculates the chemical potential MU and energy EN.
  CALL CHEM(CP, MU, EN)
  CP2 = CP * CP
! RAD() calculates the r.m.s radius RMS.
  CALL RAD(CP2, RMS)
!
  WRITE (7, 1001)
  WRITE (7, 1002)
  WRITE (7, 1001)  
  WRITE (7, 1003) ZNORM, MU/XOP, EN/XOP, RMS, CP2(N2) 
  1001 FORMAT (19X,'-----------------------------------------------------')
  1002 FORMAT (20X, 'Norm', 6X, 'Chem', 8X, 'Ener/N', 6X, '<x>', 6X, '|Psi(0)|^2')   
  1003 FORMAT ('Initial : ', 4X, F11.4, 2F12.6, 2F11.5)
!
!   OPEN(1, FILE = 'imag1d-den-init.txt')
!    CALL WRITE_DENSITY(1, CP2)
!   CLOSE(1)
!
  IF (NSTP /= 0) THEN
    GSTP =  XOP * G_1D / DFLOAT(NSTP)    
    G = 0.D0
    DO K = 1, NSTP ! NSTP iterations to introduce the nonlinearity
      G = G + GSTP
      CALL CALCNU(CP, DT)	! CALCNU() performs time propagation with non-derivative parts.
      CALL LU(CP)		! LU() performs the time iteration with space derivative alone.
      CALL NORM(CP, ZNORM)
    END DO
    CALL CHEM(CP, MU, EN)
    CP2 = CP * CP
    CALL RAD(CP2, RMS)
    WRITE (7, 1005) ZNORM, MU/XOP, EN/XOP, RMS, CP2(N2)
    1005 FORMAT('After NSTP iter.:',F8.4, 2F12.6, 2F11.5)
!
!    OPEN(2, FILE = 'imag1d-den-nstp.txt')
!     CALL WRITE_DENSITY(2, CP2)
!    CLOSE(2)
!
  ELSE
    G = XOP * G_1D
  END IF
!
  DO K = 1, NPAS ! NPAS iterations transient
    CALL CALCNU(CP, DT)
    CALL LU(CP)
    CALL NORM(CP, ZNORM)
  END DO
  IF(NPAS /= 0)THEN    
    CALL CHEM(CP, MU, EN)
    CP2 = CP * CP
    CALL RAD(CP2, RMS)
    WRITE (7, 1008) ZNORM, MU/XOP, EN/XOP, RMS, CP2(N2)
    1008 FORMAT('After NPAS iter.:',F8.4, 2F12.6, 2F11.5)
!
!    OPEN(3, FILE = 'imag1d-den-npas.txt')
!     CALL WRITE_DENSITY(3, CP2)
!    CLOSE(3)
!
  END IF
!
  DO K = 1, NRUN ! NRUN iterations to check convergence
    CALL CALCNU(CP, DT)
    CALL LU(CP)
    CALL NORM(CP, ZNORM)
  END DO
  IF(NRUN /= 0)THEN    
    CALL CHEM(CP, MU, EN)
    CP2 = CP * CP 
    CALL RAD(CP2, RMS)
    WRITE (7, 1006) ZNORM, MU/XOP, EN/XOP, RMS, CP2(N2)
    1006 FORMAT('After NRUN iter.:',F8.4, 2F12.6, 2F11.5)
  END IF
!
    OPEN(4, FILE = 'imag1d-den.txt')
      CALL WRITE_DENSITY(4, CP2)
    CLOSE(4)
!  
  CALL SYSTEM_CLOCK (CLCK_COUNTS_END, CLCK_RATE)
  CALL CPU_TIME(T2)
  WRITE (7, 1001)  
  WRITE (7,*)
  WRITE (7,'(A,I7,A)') ' Clock Time:', (CLCK_COUNTS_END - CLCK_COUNTS_BEG)/INT (CLCK_RATE,8), ' seconds'
  WRITE (7,'(A,I7,A)') '   CPU Time:', INT(T2-T1), ' seconds' 
  CLOSE (7)
!
END PROGRAM GROSS_PITAEVSKII_SSCN_1D

SUBROUTINE INITIALIZE()
! Routine that initializes the constant and variables
! Calculates the initial wave function CP
  USE COMM_DATA, ONLY : N, N2, PI
  USE GPE_DATA, ONLY : GAMMA, DX, X, X2, CP 
  IMPLICIT NONE
  REAL (8) :: PI4
  INTEGER :: J 
!
  PI4 = SQRT(SQRT(PI/GAMMA))
  FORALL (J=0:N) X(J) = (J-N2)*DX
  X2 = X*X
  CP = EXP(-X2*GAMMA/2.D0)/PI4
END SUBROUTINE INITIALIZE

SUBROUTINE CALCULATE_TRAP()
! Routine that  initializes the harmonic oscillator potential V.
  USE GPE_DATA, ONLY : XOP, X2, V, GAMMA
  IMPLICIT NONE
!
    V = GAMMA**2*X2
END SUBROUTINE CALCULATE_TRAP

SUBROUTINE COEF() 
! Calculates the coefficients needed in subroutine lu.
  USE COMM_DATA, ONLY :  NX
  USE GPE_DATA, ONLY : DX, DT   
  USE CN_DATA
  IMPLICIT NONE
  INTEGER :: J
  REAL (8) :: DX2, CAI0
!
  DX2 = DX*DX
  CAIPM = -DT/(2.D0*DX2)
  CAI0 = 1.D0 + DT/DX2
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
  USE GPE_DATA, ONLY : V, G
  IMPLICIT NONE
!-------------------------------------------------
  REAL (8), DIMENSION(0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(IN) :: DT
  REAL (8), DIMENSION(0:SIZE(CP)-1) :: P2, TMP1D, DP
!
  P2 = CP * CP
  TMP1D = DT*(V + G*P2)
!  
  CP = CP*EXP(-TMP1D)
END SUBROUTINE CALCNU

SUBROUTINE LU(CP)
! Solves the partial differential equation only with the space
! derivative term using the Crank-Nicholson method
  USE COMM_DATA, ONLY : N, NX
  USE CN_DATA
  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(INOUT) :: CP
  REAL (8), DIMENSION(0:SIZE(CP)-2) :: CBE
  REAL (8) :: CXX
  INTEGER :: I
!
  CBE(NX) = CP(N)
  DO I = NX, 1, -1
     CXX = CP(I)-(CP(I+1)-2.D0*CP(I)+CP(I-1))*CAIPM
     CBE(I-1) = CGA(I)*(CAIPM*CBE(I)-CXX)
  END DO
!-----------------------------
! Boundary condition reflecting:
   CP(0) = 0.D0
   DO I = 0, NX-1
      CP(I+1) = CAL(I)*CP(I) + CBE(I)
   END DO
   CP(N) = 0.D0
!-----------------------------
END SUBROUTINE LU

SUBROUTINE NORM(CP, ZNORM)
! Calculates the normalization of the wave function and sets it to unity.
  USE COMM_DATA, ONLY : N
  USE GPE_DATA, ONLY : DX
  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(OUT) :: ZNORM
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
  REAL (8), DIMENSION(0:SIZE(CP)-1) :: P2
!
  P2 = CP * CP
  ZNORM = SQRT(SIMP(P2, DX))
  CP = CP/ZNORM
END SUBROUTINE NORM

SUBROUTINE RAD(CP2, RMS)
!  Calculates the root mean square size RMS 
  USE COMM_DATA, ONLY : N 
  USE GPE_DATA, ONLY : DX, X2
  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(INOUT) :: CP2
  REAL (8), INTENT(OUT) :: RMS
!----------------------------------------------------------------------      
   INTERFACE
    FUNCTION SIMP(F, DX) RESULT(VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
!----------------------------------------------------------------------      
  REAL (8), DIMENSION(0:SIZE(CP2)-1) :: P2
!
  P2 = X2*CP2
  RMS = SQRT(SIMP(P2, DX))
END SUBROUTINE RAD

SUBROUTINE CHEM(CP, MU, EN)
! Calculates the chemical potential MU and energy EN.  CP is the wave
! function, V is the potential and G is the nonlinearity. 
  USE COMM_DATA, ONLY : N,NX
  USE GPE_DATA, ONLY : DX, V, G 
  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(IN) :: CP
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
  DCP = DIFF(CP, DX)
  DP2 = DCP*DCP
  P2 = CP * CP
  GP2 = G * P2
!   
      TMP1D = (V + GP2 )*P2 + DP2
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
  999 FORMAT(F12.4, E17.7E3)
END SUBROUTINE WRITE_DENSITY

!# File name : imag1d.f90
