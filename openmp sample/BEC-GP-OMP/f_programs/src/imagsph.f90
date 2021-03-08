!#  File name : imagsph.f90	!# Last modified : 14 Mar 2016
!#  Fortran program for Gross-Pitaevskii equation in three-dimensional 
!#  spherically-symmetric trap by imaginary time propagation (Fortran 90/95 Version)
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
  INTEGER, PARAMETER :: N = 3000, NX = N-1
! NPAS : Number of subsequent iterations with fixed nonlinearity, and 
! NRUN : Number of final time iterations with fixed nonlinearity
  INTEGER, PARAMETER :: NSTP=0, NPAS = 200000, NRUN = 20000
END MODULE COMM_DATA
!
MODULE GPE_DATA
  USE COMM_DATA, ONLY : N
  REAL (8), PARAMETER :: PI = 3.14159265358979D0
  REAL (8), PARAMETER :: AHO = 1.D-6     		 ! unit of length (= 1 MICRON)            
  REAL (8), PARAMETER :: Bohr_a0 =  5.2917720859D-11/AHO ! Bohr radius (scaled with AHO)
!
  REAL (8), PARAMETER :: DX = 0.0025D0, DT = 0.00002D0	 ! DX : Space step and DT : Time step
  INTEGER, PARAMETER  :: NATOMS = 2000			 ! Number of ATOMS
  REAL (8), PARAMETER :: AS = 94.351186D0*Bohr_a0	 ! Scattering length (in units of Bohr_a0) 
!   
  REAL (8), PARAMETER :: G0 = 4.D0*PI*AS*NATOMS		 ! G0 : 3-dimensional nonlinearity
!   
! OPTION and XOP decides which equation to be solved. 
! OPTION=1 -> Solves -psi_xx+V(x)psi+G0|psi|^2 psi=i psi_t
! OPTION=2 -> Solves [-psi_xx+V(x)psi]/2+G0|psi|^2 psi=i psi_t
  INTEGER, PARAMETER :: OPTION = 2
! X(0:N) : Space mesh, V(0:N) : Potential
! CP(0:N) : Wave function (COMPLEX)
  REAL (8), DIMENSION(0:N) :: X, X2, V, CP
  REAL (8) :: G, XOP
END MODULE GPE_DATA
!
MODULE CN_DATA
  USE COMM_DATA, ONLY : N, NX
  REAL (8), DIMENSION(N) :: CAL, CGA
  REAL (8) :: CAIPM
END MODULE CN_DATA
!
PROGRAM GROSS_PITAEVSKII_SSCN_SPH
  USE COMM_DATA
  USE GPE_DATA
  IMPLICIT NONE
! Subroutine INITIALIZE() used to initialize the space mesh X(I), 
! potential V(I) and the initial wave function. Subroutine COEF() used to
! generate the coefficients for the Crank-Nicholson Scheme. The routine
! CALCNU() performs time progation for the non-derivative part and LU() performs
! time propagation of derivative part. NORM() calculates the norm and 
! normalizes the wave function, CHEM() and RAD() are used to calculate the 
! chemical potential, energy and the rms radius, respectively. The functions 
! DIFF(), used to calculate the space derivatives of the wave function used 
! in CHEM(). The function SIMP() does the integration by Simpson's rule.
!------------------------ INTERFACE BLOCKS -----------------------
  INTERFACE 
    SUBROUTINE INITIALIZE()
      IMPLICIT NONE
    END SUBROUTINE INITIALIZE
  END INTERFACE
!
  INTERFACE 
    SUBROUTINE CALCNU(CP, DT)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(INOUT) :: CP
      REAL (8), INTENT(IN) :: DT
    END SUBROUTINE CALCNU
  END INTERFACE
!
  INTERFACE 
    SUBROUTINE LU(CP)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(INOUT) :: CP
    END SUBROUTINE LU
  END INTERFACE
!
  INTERFACE
    SUBROUTINE NORM(CP, ZNORM)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(INOUT) :: CP
      REAL (8), INTENT(OUT) :: ZNORM
    END  SUBROUTINE NORM
  END INTERFACE
!
  INTERFACE
    SUBROUTINE RAD(CP, RMS)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(INOUT) :: CP
      REAL (8), INTENT(OUT) :: RMS
    END  SUBROUTINE RAD
  END INTERFACE
!
   INTERFACE
     SUBROUTINE CHEM(CP, MU, EN)
       IMPLICIT NONE
       REAL (8), DIMENSION(0:), INTENT(IN) :: CP
       REAL (8), INTENT(OUT) :: MU, EN
     END  SUBROUTINE CHEM
   END INTERFACE
!------------------------ END INTERFACE BLOCKS -------------------
  INTEGER :: I, K
  REAL (8) :: GSTP, ZNORM, MU, EN, RMS,T1,T2
 INTEGER :: CLCK_COUNTS_BEG, CLCK_COUNTS_END, CLCK_RATE
!
  CALL SYSTEM_CLOCK ( CLCK_COUNTS_BEG, CLCK_RATE )
  CALL CPU_TIME(T1)
  SELECT CASE (OPTION)
    CASE (1)
      XOP = 1.0D0
    CASE (2)
      XOP = 2.0D0
    CASE (:0,4:)
      PRINT *, 'ERROR: Wrong option', option
      STOP
  END SELECT
  G = 0.0D0
  OPEN(7, FILE = 'imagsph-out.txt')
  WRITE(7,900) OPTION
  WRITE(7,*)  
  WRITE(7,906) NATOMS, AHO  
  WRITE(7,907) AS/Bohr_a0   
  WRITE(7,903) G0   
  WRITE(7,*)    
  WRITE(7,901) N
  WRITE(7,908) DX   
  WRITE(7,902) NSTP, NPAS, NRUN  
  WRITE(7,904) DT 
  WRITE(7,*)
  900 FORMAT(' Imaginary time propagation spherically-symmetric trap,   OPTION = ',I3 )
  906 FORMAT('  Number of Atoms N =',I10,', Unit of length AHO =',F12.8,' m')
  907 FORMAT('  Scattering length a = ',F9.2,'*a0')
  903 FORMAT('  Nonlinearity G = ', F16.7)  
  901 FORMAT(' # Space Stp:  N = ', I8)
  908 FORMAT('  Space Step: DX = ', F10.6)
  902 FORMAT(' # Time Stp : NSTP = ',I9,', NPAS = ',I9,', NRUN = ',I9)
  904 FORMAT('   Time Step:   DT = ', F10.6)
! INITIALIZE() initializes the starting normalized wave function 'CP' and
! harmonic oscillator potential V
  CALL INITIALIZE()
! COEF() defines the coefficients of Crank-Nicholson Scheme.
  CALL COEF()
  
!   OPEN(1, FILE = 'imagsph-den-init.txt')
!     DO I = 0, N ! write the initial wave function
!     IF(I.NE.0)   WRITE(1,905) X(I), (CP(I)/X(I))**2
!     IF(I.EQ.0)   WRITE(1,905) X(I), (CP(1)/X(1))**2
!     END DO
!   CLOSE(1)
!
! NORM() calculates and restore normalization
  CALL NORM(CP, ZNORM)
! CHEM() calculates the chemical potential MU and energy EN.
  CALL CHEM(CP, MU, EN)
! RAD() calculates the r.m.s radius RMS
  CALL RAD(CP, RMS)
!   
  WRITE (7, 1001)
  WRITE (7, 1002)
  WRITE (7, 1001)
  WRITE (7, 1003) ZNORM, MU/XOP, EN/XOP, RMS, (CP(1)/X(1))**2
  1001 FORMAT (19X,'-----------------------------------------------------')
  1002 FORMAT (20X, 'Norm', 6X, 'Chem', 8X, 'Ener', 6X, '<r>', 6X, '|Psi(0)|^2')
  1003 FORMAT ('Initial : ', 4X, F11.4, 2F12.6, 3F11.5)

 IF (NSTP /= 0) THEN
    GSTP =  XOP * G0 / DFLOAT(NSTP)    
    G = 0.D0
    DO K = 1, NSTP ! NSTP iterations to introduce the nonlinearity
      G = G + GSTP
      CALL CALCNU(CP, DT)	! CALCNU() performs time propagation with non-derivative parts.
      CALL LU(CP)		! LU() performs the time iteration with space derivative alone.
      CALL NORM(CP, ZNORM)
    END DO
    CALL CHEM(CP, MU, EN)
    CALL RAD(CP, RMS)
    WRITE (7, 1005) ZNORM, MU/XOP, EN/XOP, RMS, (CP(1)/X(1))**2
    1005 FORMAT('After NSTP iter.:',F8.4, 2F12.6, 2F11.5)
!     
!   OPEN(2, FILE = 'imagsph-den-nstp.txt')
!     DO I = 0, N ! write the  wave function after NSTP iterations
!     IF(I.NE.0)   WRITE(2,905) X(I), (CP(I)/X(I))**2
!     IF(I.EQ.0)   WRITE(2,905) X(I), (CP(1)/X(1))**2
!     END DO
!   CLOSE(2)    
! 
  ELSE
    G = XOP * G0
  END IF 
! 
  DO K = 1, NPAS ! NPAS iterations transient
     CALL CALCNU(CP, DT)
     CALL LU(CP)
     CALL NORM(CP, ZNORM)
  END DO
  IF(NPAS /= 0)THEN  
    CALL CHEM(CP, MU, EN)
    CALL RAD(CP, RMS)
    WRITE (7, 1008) ZNORM, MU/XOP, EN/XOP, RMS, (CP(1)/X(1))**2
    1008 FORMAT('After NPAS iter.:', F8.4, 2F12.6, 2F11.5)
!   
!    OPEN(3, FILE = 'imagsph-den-npas.txt')
!      DO I = 0, N ! write the wave function after NPAS iterations
!      IF(I.NE.0)   WRITE(3,905) X(I), (CP(I)/X(I))**2
!      IF(I.EQ.0)   WRITE(3,905) X(I), (CP(1)/X(1))**2
!      END DO
!    CLOSE(3)
! 
  END IF  

  DO K = 1, NRUN ! NRUN iterations to check convergence
     CALL CALCNU(CP, DT)
     CALL LU(CP)
     CALL NORM(CP, ZNORM)
  END DO
  IF(NRUN /= 0)THEN   
    CALL CHEM(CP, MU, EN)
    CALL RAD(CP, RMS)
    WRITE (7, 1006) ZNORM, MU/XOP, EN/XOP, RMS,  (CP(1)/X(1))**2
    1006 FORMAT('After NRUN iter.:',F8.4, 2F12.6, 2F11.5)
  END IF 
!     
    OPEN(4, FILE = 'imagsph-den.txt')
      DO I = 0, N ! write the final wave function
      IF(I.NE.0)   WRITE(4,905) X(I), (CP(I)/X(I))**2
      IF(I.EQ.0)   WRITE(4,905) X(I), (CP(1)/X(1))**2
      END DO
    CLOSE(4)
!       
  CALL SYSTEM_CLOCK (CLCK_COUNTS_END, CLCK_RATE)
  CALL CPU_TIME(T2)
  WRITE (7, 1001)
  WRITE (7,*)
  WRITE (7,'(A,I7,A)') ' Clock Time: ', (CLCK_COUNTS_END - CLCK_COUNTS_BEG)/INT (CLCK_RATE,8), ' seconds'
  WRITE (7,'(A,I7,A)') '   CPU Time: ', INT(T2-T1), ' seconds' 
  905 FORMAT(F12.6, F16.8) 
  CLOSE (7)
!   
END PROGRAM GROSS_PITAEVSKII_SSCN_SPH

SUBROUTINE INITIALIZE()
! Routine that initializes the constant and variables.
! Calculates the potential term V and the initial wave function CP
  USE COMM_DATA
  USE GPE_DATA
  IMPLICIT NONE
  REAL (8) :: PI4
  INTEGER :: J
! 
  PI4 = SQRT(PI*SQRT(PI)) ! (PI)^(3/4)
  FORALL (J=0:N) X(J) = J*DX
  X2 = X*X
!   
  V = X2
  CP = X*EXP(-X2/2.0D0)/PI4     
END SUBROUTINE INITIALIZE

SUBROUTINE COEF()
! Calculates the coefficients needed in subroutine LU.
  USE COMM_DATA
  USE GPE_DATA
  USE CN_DATA
  IMPLICIT NONE
  INTEGER :: J
  REAL (8) :: DX2, CAI0
!   
  DX2 = DX*DX
  CAIPM = -DT/(2.0D0*DX2)
  CAI0 = 1.0D0 + DT/DX2
  CAL(NX) = 0.0D0
  CGA(NX) = -1.0D0/CAI0
  DO J = NX, 1, -1
     CAL(J-1) = CGA(J)*CAIPM
     CGA(J-1) = -1.0D0/(CAI0+CAIPM*CAL(J-1))      
  END DO   
END SUBROUTINE COEF

SUBROUTINE CALCNU(CP, DT)
! Solves the partial differential equation with the potential and the
! nonlinear term.      
  USE GPE_DATA, ONLY : X, V, G
  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(IN) :: DT
  REAL (8), DIMENSION(1:SIZE(CP)-1) :: TMP1D, P, P2
!   
  P = CP(1:)/X(1:)
  P2 = P*P
!   
  TMP1D = DT*(V(1:)+G*P2)
  CP(0) = 0.0D0
  CP(1:) = CP(1:)*EXP(-TMP1D)
END SUBROUTINE CALCNU

SUBROUTINE LU(CP)
! Solves the partial differential equation only with the space
! derivative term using the Crank-Nicholson method
  USE COMM_DATA, ONLY : N, NX
  USE CN_DATA
  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(INOUT) :: CP
  REAL (8), DIMENSION(0:N) :: CBE
  REAL (8) :: CXX
  INTEGER :: I
!   
  CBE(NX) = CP(N)
  DO I = NX, 1, -1
     CXX = CP(I)-(CP(I+1)-2.0D0*CP(I)+CP(I-1))*CAIPM
     CBE(I-1) = CGA(I)*(CAIPM*CBE(I)-CXX)
  END DO
  CP(0) = 0.0D0
  DO I = 0, NX
     CP(I+1) = CAL(I)*CP(I) + CBE(I)
  END DO
END SUBROUTINE LU

SUBROUTINE NORM(CP, ZNORM)
! Calculates the normalization of the wave function and sets it to
! unity.
  USE COMM_DATA, ONLY : N
  USE GPE_DATA, ONLY : DX, PI
  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(OUT) :: ZNORM
!-----------------------------------------------------------  
  INTERFACE
    FUNCTION SIMP(F, DX) RESULT (VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
!-----------------------------------------------------------  
  REAL (8), DIMENSION(0:N) :: P2
!   
  P2 = CP*CP
  ZNORM = SQRT(4.0D0*PI*SIMP(P2, DX))
  CP = CP/ZNORM
END SUBROUTINE NORM

SUBROUTINE RAD(CP, RMS)
! Calculates the root mean square radius RMS.
  USE COMM_DATA, ONLY : N
  USE GPE_DATA, ONLY : DX, X2, PI
  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(OUT) :: RMS
!-----------------------------------------------------------  
  INTERFACE
    FUNCTION SIMP(F, DX) RESULT (VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
!-----------------------------------------------------------  
  REAL (8), DIMENSION(0:SIZE(CP)-1) :: CP2
!   
  CP2 = X2*CP*CP
  RMS = SQRT(4.0D0*PI*SIMP(CP2, DX))
END SUBROUTINE RAD

SUBROUTINE CHEM(CP, MU, EN)
!  Calculates the chemical potential MU and energy EN.  CP is the wave
!  function defined at N+1 points at steps of DX, V is the potential and
!  G is the nonlinearity.
  USE COMM_DATA, ONLY : N
  USE GPE_DATA, ONLY : DX, X2, V, G, PI
  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(IN) :: CP
  REAL (8), INTENT(OUT) :: MU, EN
!-----------------------------------------------------------  
  INTERFACE
    PURE FUNCTION DIFF(P, DX) RESULT(DP)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: P
      REAL (8), INTENT(IN) :: DX
      REAL (8), DIMENSION(0:SIZE(P)-1) :: DP
    END FUNCTION DIFF
  END INTERFACE
!-----------------------------------------------------------  
  INTERFACE
    FUNCTION SIMP(F, DX) RESULT (VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
!-----------------------------------------------------------  
  REAL (8), DIMENSION(0:N) :: DCP, TMP1D
  REAL (8), DIMENSION(1:N) :: DP2, P2, GP2
!   
  DCP = DIFF(CP, DX)
  P2 = CP(1:)*CP(1:)
  DP2 = DCP(1:)*DCP(1:)
  GP2 = G*P2
  TMP1D(0) = DP2(1)
  TMP1D(1:) = (V(1:) + GP2/X2(1:))*P2 + DP2
  MU = 4.0D0*PI*SIMP(TMP1D, DX)
  TMP1D(1:) = (V(1:) + GP2/X2(1:)/2.0D0)*P2 + DP2
  EN = 4.0D0*PI*SIMP(TMP1D, DX)
END SUBROUTINE CHEM

FUNCTION SIMP(F, DX) RESULT (VALUE)
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
  VALUE = DX*(F(0) + 4.0D0*F1 + 2.0D0*F2 + F(N))/3.0D0
END FUNCTION SIMP
!
PURE FUNCTION DIFF(P, DX) RESULT (DP)
! Computes the first derivative DP of P using
! Richardsonextrapolation formula. The derivative at the
! boundaries are assumed to be zero
  IMPLICIT NONE
  REAL (8), DIMENSION(0:), INTENT(IN) :: P
  REAL (8), INTENT(IN) :: DX
  REAL (8), DIMENSION(0:SIZE(P)-1) :: DP
  INTEGER :: I, N
  N = SIZE(P) - 1
  DP(0) = 0.0D0
  DP(1) = (P(2) - P(0))/(2.0D0*DX)
  FORALL(I =2:N-2)
    DP(I) = (P(I-2)-8.0D0*P(I-1)+8.0D0*P(I+1)-P(I+2))/(12.0D0*DX)
  END FORALL
  DP(N-1) = (P(N) - P(N-2))/(2.0D0*DX)
  DP(N) = 0.0D0
END FUNCTION DIFF
 
!# File name : imagsph.f90
