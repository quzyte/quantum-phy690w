!#  File name : imag3d.f90	!# Last modified : 14 Mar 2016
!#  Fortran program for Gross-Pitaevskii equation in three-dimensional 
!#  anisotropic trap by imaginary time propagation (Fortran 90/95 Version)
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
! gfortran -O3 -fopenmp -fno-automatic
!
!Important note: if you encounter segmentation fault erros when running Fortran
!programs on Linux, you should increase the stack size. Running the following
!command should resolve the issue:

!   ulimit -s unlimited          [if you are using bash shell]
!   limit stacksize unlimited    [if you are using csh or tcsh shell]


!
MODULE COMM_DATA
! NX, NY, NZ : Number of space mesh points (X, Y and Z)  
  INTEGER, PARAMETER :: NX = 240, NXX = NX-1, NX2 = NX/2
  INTEGER, PARAMETER :: NY = 200, NYY = NY-1, NY2 = NY/2
  INTEGER, PARAMETER :: NZ = 160, NZZ = NZ-1, NZ2 = NZ/2
! NSTP : Number of iterations to introduce the nonlinearity.
! NPAS : Number of subsequent iterations with fixed nonlinearity.
! NRUN : Number of final time iterations with fixed nonlinearity.
  INTEGER, PARAMETER :: NSTP = 0, NPAS = 5000, NRUN = 500
  INTEGER, PARAMETER :: NO_OF_THREADS = 4 ! CPU: Quad Core, Single Socket
  INTEGER, PARAMETER :: NNN = NX*NY*NZ
  REAL (8), PARAMETER :: PI = 3.14159265358979D0
END MODULE COMM_DATA

MODULE GPE_DATA
  USE COMM_DATA, ONLY : NX, NY, NZ, PI
  REAL (8), PARAMETER :: AHO = 1.D-6     		 ! Unit of length (= 1 MICRON)            
  REAL (8), PARAMETER :: Bohr_a0 =  5.2917720859D-11/AHO ! Bohr radius (scaled with AHO)
!  
  REAL (8), PARAMETER :: DX = 0.05D0, DY = DX, DZ = DX	 ! DX, DY, DZ : SPACE STEPS
  REAL (8), PARAMETER :: DT = 0.0004D0			 ! DT : TIME STEP  
  INTEGER, PARAMETER :: NATOMS = 1000			 ! Number of Atoms
  REAL (8), PARAMETER :: AS = 67.530979D0*Bohr_a0	 ! Scattering length (in units of Bohr_a0)
  REAL (8), PARAMETER :: GAMMA = 1.D0, ANU = 1.414214D0, LAMBDA = 2.D0	! GAMMA, NU, LAMBDA : Parameteres of Trap in x, y and z directions
!   
  REAL (8), PARAMETER :: G0 = 4.D0*PI*AS*NATOMS! 44.907D0		! Three-dimensional nonlinearity 
!
! OPTION decides which equation to be solved.
! OPTION=1 Solves -psi_xx-psi_yy-psi_zz+V(x,y,z)psi+G0|psi|^2 psi =i psi_t
! OPTION=2 Solves [-psi_xx-psi_yy-psi_zz+V(x,y,z)psi]/2+G0|psi|^2 psi =i psi_t
  INTEGER, PARAMETER :: OPTION = 2 
! X(0:NX), Y(0:NY), Z(0:NZ) : Space mesh, V(0:NX,0:NY,0:NZ) : Potential, 
! CP(0:NX,0:NY,0:NZ) : Wave function  
  REAL (8), DIMENSION(0:NX) :: X, X2
  REAL (8), DIMENSION(0:NY) :: Y, Y2
  REAL (8), DIMENSION(0:NZ) :: Z, Z2 
  REAL (8), DIMENSION(0:NX, 0:NY, 0:NZ) :: V, R2, CP
  REAL (8) :: G,   XOP
END MODULE GPE_DATA

MODULE CN_DATA
  USE COMM_DATA, ONLY : NX, NY, NZ
  REAL (8), DIMENSION(0:NX) :: CALA, CGAA
  REAL (8), DIMENSION(0:NY) :: CALB, CGAB
  REAL (8), DIMENSION(0:NZ) :: CALC, CGAC
  REAL (8) :: CT0X, CT0Y, CT0Z
  REAL (8) :: CA0, CB0, CC0, CA0R, CB0R, CC0R
END MODULE CN_DATA 

PROGRAM GROSS_PITAEVSKII_SSCN_3D
  USE COMM_DATA, ONLY : NX, NY, NZ, NX2, NY2, NZ2, NSTP, NPAS, NRUN
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
!------------------------ interface blocks -----------------------       
  INTERFACE 
    SUBROUTINE INITIALIZE()
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
    END SUBROUTINE COEF
  END INTERFACE
!------------------------
  INTERFACE 
    SUBROUTINE CALCNU(CP, DT)
      REAL (8), DIMENSION(0:, 0:, 0:), INTENT(INOUT) :: CP
      REAL (8), INTENT(IN) :: DT
    END SUBROUTINE CALCNU
  END INTERFACE
!------------------------
  INTERFACE 
    SUBROUTINE LUX(CP)
      REAL (8), DIMENSION(0:, 0:, 0:), INTENT(INOUT) :: CP
    END SUBROUTINE LUX
  END INTERFACE
!------------------------
  INTERFACE 
    SUBROUTINE LUY(CP)
      REAL (8), DIMENSION(0:, 0:, 0:), INTENT(INOUT) :: CP
    END SUBROUTINE LUY
  END INTERFACE
!------------------------ 
  INTERFACE 
    SUBROUTINE LUZ(CP)
      REAL (8), DIMENSION(0:, 0:, 0:), INTENT(INOUT) :: CP
    END SUBROUTINE LUZ
  END INTERFACE
!------------------------ 
  INTERFACE
    SUBROUTINE NORM(CP, ZNORM)
      REAL (8), DIMENSION(0:, 0:, 0:), INTENT(INOUT) :: CP
      REAL (8), INTENT(OUT) :: ZNORM
    END  SUBROUTINE NORM
  END INTERFACE
!------------------------ 
  INTERFACE
    SUBROUTINE RAD(CP2, RMS)
      REAL (8), DIMENSION(0:, 0:, 0:), INTENT(IN) :: CP2
      REAL (8), DIMENSION(:), INTENT(OUT) :: RMS
    END  SUBROUTINE RAD
  END INTERFACE
!------------------------ 
  INTERFACE
    SUBROUTINE CHEM(CP, MU, EN)
      REAL (8), DIMENSION(0:, 0:, 0:), INTENT(IN) :: CP
      REAL (8), INTENT(OUT) :: MU, EN
    END SUBROUTINE CHEM
  END INTERFACE
!------------------------ 
  INTERFACE
    SUBROUTINE WRITE_2DXY(FUNIT, U2)
      IMPLICIT NONE
      INTEGER :: FUNIT
      REAL (8), DIMENSION(0:,0:,0:), INTENT(IN) :: U2
    END  SUBROUTINE WRITE_2DXY
  END INTERFACE
!----------------------------
  INTERFACE
    SUBROUTINE WRITE_2DXZ(FUNIT, U2)
      IMPLICIT NONE
      INTEGER :: FUNIT
      REAL (8), DIMENSION(0:,0:,0:), INTENT(IN) :: U2
    END  SUBROUTINE WRITE_2DXZ
  END INTERFACE
!----------------------------
  INTERFACE
    SUBROUTINE DEN2DYZ(FUNIT, U2, DE2DYZ)
      IMPLICIT NONE
      INTEGER :: FUNIT
      REAL (8), DIMENSION(0:,0:,0:), INTENT(IN) :: U2
      REAL (8), DIMENSION(0:,0:), INTENT(OUT) :: DE2DYZ
    END  SUBROUTINE DEN2DYZ
  END INTERFACE
!------------------------ 
  INTERFACE
    SUBROUTINE DEN1DY(FUNIT, U2, DE1DY)
      IMPLICIT NONE
      INTEGER :: FUNIT
      REAL (8), DIMENSION(0:,0:,0:), INTENT(IN) :: U2
      REAL (8), DIMENSION(0:), INTENT(OUT) :: DE1DY
    END  SUBROUTINE DEN1DY
  END INTERFACE
!------------------------
  INTERFACE
    SUBROUTINE WRITE_3D(FUNIT, U2)
      IMPLICIT NONE
      INTEGER :: FUNIT
      REAL (8), DIMENSION(0:,0:,0:), INTENT(IN) :: U2
    END  SUBROUTINE WRITE_3D
  END INTERFACE
!------------------------ 
  INTERFACE
    SUBROUTINE DEN1DZ(FUNIT, U2, DE1DZ)
      IMPLICIT NONE
      INTEGER :: FUNIT
      REAL (8), DIMENSION(0:,0:,0:), INTENT(IN) :: U2
      REAL (8), DIMENSION(0:), INTENT(OUT) :: DE1DZ
    END  SUBROUTINE DEN1DZ
  END INTERFACE
!------------------------ 
  INTERFACE
    SUBROUTINE DEN1DX(FUNIT, U2, DE1DX)
      IMPLICIT NONE
      INTEGER :: FUNIT
      REAL (8), DIMENSION(0:,0:,0:), INTENT(IN) :: U2
      REAL (8), DIMENSION(0:), INTENT(OUT) :: DE1DX
    END  SUBROUTINE DEN1DX
  END INTERFACE
!------------------------ 
  INTERFACE
    SUBROUTINE DEN2DXZ(FUNIT, U2, DE2DXZ)
      IMPLICIT NONE
      INTEGER :: FUNIT
      REAL (8), DIMENSION(0:,0:,0:), INTENT(IN) :: U2
      REAL (8), DIMENSION(0:,0:), INTENT(OUT) :: DE2DXZ
    END  SUBROUTINE DEN2DXZ
  END INTERFACE
!------------------------ 
  INTERFACE
    SUBROUTINE DEN2DXY(FUNIT, U2, DE2DXY)
      IMPLICIT NONE
      INTEGER :: FUNIT
      REAL (8), DIMENSION(0:,0:,0:), INTENT(IN) :: U2
      REAL (8), DIMENSION(0:,0:), INTENT(OUT) :: DE2DXY
    END  SUBROUTINE DEN2DXY
  END INTERFACE
!-------------------- END INTERFACE BLOCKS -----------------------
  INTEGER :: K
  REAL (8), DIMENSION(0:NX, 0:NY, 0:NZ) :: CP2
!  <r> = RMS(1), <x> = RMS(2), <y> = RMS(3), <z> = RMS(4)
  REAL (8), DIMENSION(4) :: RMS
  REAL (8) :: ZNORM, MU, EN, T1, T2 
  REAL (8) :: GSTP  
  REAL (8), DIMENSION(0:NZ) :: DEN1Z 
  REAL (8), DIMENSION(0:NX, 0:NY) :: DEN2XY
  REAL (8), DIMENSION(0:NX) :: DEN1X 
  REAL (8), DIMENSION(0:NX, 0:NZ) :: DEN2XZ
  REAL (8), DIMENSION(0:NY, 0:NZ) :: DEN2YZ
  REAL (8), DIMENSION(0:NY) :: DEN1Y
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
  OPEN(7, FILE = 'imag3d-out.txt')
  OPEN(4, FILE = 'imag3d-rms.txt')
!
  WRITE(7,900) OPTION
  WRITE(4,900) OPTION
  WRITE(7,*) 
  WRITE(4,*) 
  WRITE(7,901) NATOMS, AHO
  WRITE(7,902) AS/Bohr_a0 
  WRITE(7,903) G0 
  WRITE(7,904) GAMMA, ANU, LAMBDA
  WRITE(7,*)
  WRITE(7,905) NX, NY, NZ
  WRITE(7,906) DX, DY, DZ
  WRITE(7,907) NSTP, NPAS, NRUN  
  WRITE(7,908) DT 
  WRITE(7,*)
  900 FORMAT(' Imaginary time propagation 3D, OPTION =',I3)
  901 FORMAT('  Number of Atoms N =',I10,', Unit of length AHO =',F12.8,' m')
  902 FORMAT('  Scattering length a = ',F9.2,'*a0')
  903 FORMAT('  Nonlinearity G_3D =',F16.7 )
  904 FORMAT('  Parameters of trap: GAMMA =',F7.4, ', NU =',F7.4, ', LAMBDA =',F7.4)
  905 FORMAT(' # Space Stp: NX = ', I8, ', NY = ', I8, ', NZ = ', I8)
  906 FORMAT('  Space Step: DX = ', F10.6, ', DY = ', F10.6, ', DZ = ', F10.6)
  907 FORMAT(' # Time Stp : NSTP = ',I9,', NPAS = ',I9,', NRUN = ',I9)
  908 FORMAT('   Time Step:   DT = ', F10.6)
!
! INITIALIZE() initializes the starting normalized wave function 'CP'.
  CALL INITIALIZE()
! CALCULATE_TRAP() initializes the harmonic oscillator potential V.
  CALL CALCULATE_TRAP()
! COEF() defines the coefficients of Crank-Nicholson Scheme.
  CALL COEF()
! NORM() calculates norm and restores normalization.
  CALL NORM(CP, ZNORM)
! CHEM() calculates the chemical potential MU and energy EN.
  CALL CHEM(CP, MU, EN)
  CP2 = CP*CP
! RAD() calculates the r.m.s radius RMS.
  CALL RAD(CP2, RMS)
!
  WRITE (7, 1001)
  WRITE (7, 1002)
  WRITE (7, 1001)
  WRITE (7, 1003) ZNORM, MU/XOP, EN/XOP, RMS(1), CP2(NX2, NY2, NZ2) 
  WRITE (4, 1001)
  WRITE (4, 1012)
  WRITE (4, 1001)
  WRITE (4, 1013) RMS(1:4) 
  1001 FORMAT (19X,'-----------------------------------------------------')
  1002 FORMAT (20X, 'Norm', 6X, 'Chem', 8X, 'Ener/N', 5X, '<r>', 5X, '|Psi(0,0,0)|^2')
  1003 FORMAT ('Initial : ', 4X, F11.4, 2F12.5, 2F11.5)
  1012 FORMAT ('Values of rms size:', 5X, '<r>', 10X, '<x>', 10X, '<y>', 10X, '<z>')
  1013 FORMAT (11X, 'Initial:', F12.5, 3F13.5)
! 
!    OPEN(31, FILE = 'imag3d-den-init.txt')
!    CALL WRITE_3D(31, CP2)
!    CLOSE(31)
!  
!    OPEN(11, FILE = 'imag3d-den-init1d_z.txt')
!    OPEN(12, FILE = 'imag3d-den-init1d_x.txt')
!    OPEN(13, FILE = 'imag3d-den-init1d_y.txt')
!    CALL DEN1DZ(11, CP2, DEN1Z)
!    CALL DEN1DX(12, CP2, DEN1X)
!    CALL DEN1DY(13, CP2, DEN1Y)
!    CLOSE(11)
!    CLOSE(12)
!    CLOSE(13)
! 
  !OPEN(21, FILE = 'imag3d-den-init2d_xz.txt')
  !OPEN(22, FILE = 'imag3d-den-init2d_xy.txt')
  !OPEN(23, FILE = 'imag3d-den-init2d_yz.txt')
  !CALL DEN2DXZ(21, CP2, DEN2XZ)
  !CALL DEN2DXY(22, CP2, DEN2XY)
  !CALL DEN2DYZ(23, CP2, DEN2YZ)
  !CLOSE(21)
  !CLOSE(22)
  !CLOSE(23)

  !OPEN(32, FILE = 'imag3d-den-init3d_x0z.txt')
  !OPEN(33, FILE = 'imag3d-den-init3d_xy0.txt')
  !CALL WRITE_2DXZ(32, CP2)
  !CALL WRITE_2DXY(33, CP2)
  !CLOSE(32)
  !CLOSE(33)
!
  IF(NSTP /= 0)THEN
    GSTP =  XOP*G0/DFLOAT(NSTP)
    G = 0.D0
   DO K = 1, NSTP ! NSTP iterations to introduce the nonlinearity
     G = G + GSTP
     CALL CALCNU(CP, DT)	! CALCNU() performs time propagation with non-derivative parts.
     CALL LUX(CP)		! LU() performs the time iteration with space derivative alone.
     CALL LUY(CP)
     CALL LUZ(CP)
     CALL NORM(CP, ZNORM)
   END DO
    CALL CHEM(CP, MU, EN)
    CP2 = CP*CP
    CALL RAD(CP2, RMS)
    WRITE (7, 1005) ZNORM, MU/XOP, EN/XOP, RMS(1), CP2(NX2, NY2, NZ2)
    WRITE (4, 1015) RMS(1:4)
    1005 FORMAT('After NSTP iter.:', F8.4, 2F12.5, 2F11.5)
    1015 FORMAT (2X, 'After NSTP iter.:', F12.5, 3F13.5)
!
!   OPEN(31, FILE = 'imag3d-den-nstp.txt')
!   CALL WRITE_3D(31, CP2)
!   CLOSE(31)
! 
!   OPEN(11, FILE = 'imag3d-den-nstp1d_z.txt')
!   OPEN(12, FILE = 'imag3d-den-nstp1d_x.txt')
!   OPEN(13, FILE = 'imag3d-den-nstp1d_y.txt')
!   CALL DEN1DZ(11, CP2, DEN1Z)
!   CALL DEN1DX(12, CP2, DEN1X)
!   CALL DEN1DY(13, CP2, DEN1Y)
!   CLOSE(11)
!   CLOSE(12)
!   CLOSE(13)
! 
  !OPEN(21, FILE = 'imag3d-den-nstp2d_xz.txt')
  !OPEN(22, FILE = 'imag3d-den-nstp2d_xy.txt')
  !OPEN(23, FILE = 'imag3d-den-nstp2d_yz.txt')
  !CALL DEN2DXZ(21, CP2, DEN2XZ)
  !CALL DEN2DXY(22, CP2, DEN2XY)
  !CALL DEN2DYZ(23, CP2, DEN2YZ)
  !CLOSE(21)
  !CLOSE(22)
  !CLOSE(23)

  !OPEN(32, FILE = 'imag3d-den-nstp3d_x0z.txt')
  !OPEN(33, FILE = 'imag3d-den-nstp3d_xy0.txt')
  !CALL WRITE_2DXZ(32, CP2)
  !CALL WRITE_2DXY(33, CP2)
  !CLOSE(32)
  !CLOSE(33)
!
  ELSE
     G= XOP*G0
  END IF
!
  DO K = 1, NPAS ! NPAS iterations transient
     CALL CALCNU(CP, DT)
     CALL LUX(CP)
     CALL LUY(CP)
     CALL LUZ(CP)
     CALL NORM(CP, ZNORM)
  END DO
  IF(NPAS /= 0)THEN  
    CALL CHEM(CP, MU, EN)
    CP2 = CP*CP
    CALL RAD(CP2, RMS)
    WRITE (7, 1006) ZNORM, MU/XOP, EN/XOP, RMS(1), CP2(NX2, NY2, NZ2)
    WRITE (4, 1016) RMS(1:4)
    1006 FORMAT('After NPAS iter.:',F8.4, 2F12.5, 2F11.5)
    1016 FORMAT (2X, 'After NPAS iter.:', F12.5, 3F13.5)
!
!    OPEN(31, FILE = 'imag3d-den-npas.txt')
!    CALL WRITE_3D(31, CP2)
!    CLOSE(31)
! 
!    OPEN(11, FILE = 'imag3d-den-npas1d_z.txt')
!    OPEN(12, FILE = 'imag3d-den-npas1d_x.txt')
!    OPEN(13, FILE = 'imag3d-den-npas1d_y.txt')
!    CALL DEN1DZ(11, CP2, DEN1Z)
!    CALL DEN1DX(12, CP2, DEN1X)
!    CALL DEN1DY(13, CP2, DEN1Y)
!    CLOSE(11)
!    CLOSE(12)
!    CLOSE(13)
! 
   !OPEN(21, FILE = 'imag3d-den-npas2d_xz.txt')
   !OPEN(22, FILE = 'imag3d-den-npas2d_xy.txt')
   !OPEN(23, FILE = 'imag3d-den-npas2d_yz.txt')
   !CALL DEN2DXZ(21, CP2, DEN2XZ)
   !CALL DEN2DXY(22, CP2, DEN2XY)
   !CALL DEN2DYZ(23, CP2, DEN2YZ)
   !CLOSE(21)
   !CLOSE(22)
   !CLOSE(23)

   !OPEN(32, FILE = 'imag3d-den-npas3d_x0z.txt')
   !OPEN(33, FILE = 'imag3d-den-npas3d_xy0.txt')
   !CALL WRITE_2DXZ(32, CP2)
   !CALL WRITE_2DXY(33, CP2)
   !CLOSE(32)
   !CLOSE(33)
!   
  END IF  
!
  DO K = 1, NRUN ! NRUN iterations to check convergence
     CALL CALCNU(CP, DT)
     CALL LUX(CP)
     CALL LUY(CP)
     CALL LUZ(CP)
     CALL NORM(CP, ZNORM)
  END DO
  IF(NRUN /= 0)THEN    
  CALL CHEM(CP, MU, EN)
  CP2 = CP*CP
  CALL RAD(CP2, RMS)
  WRITE (7, 1007) ZNORM, MU/XOP, EN/XOP, RMS(1), CP2(NX2, NY2, NZ2)
  WRITE (4, 1017) RMS(1:4)
  1007 FORMAT('After NRUN iter.:',F8.4, 2F12.5, 2F11.5)
  1017 FORMAT (2X, 'After NRUN iter.:', F12.5, 3F13.5)
  END IF 
! 
  OPEN(31, FILE = 'imag3d-den.txt')
  CALL WRITE_3D(31, CP2)
  CLOSE(31)

  OPEN(11, FILE = 'imag3d-den1d_z.txt')
  OPEN(12, FILE = 'imag3d-den1d_x.txt')
  OPEN(13, FILE = 'imag3d-den1d_y.txt')
  CALL DEN1DZ(11, CP2, DEN1Z)
  CALL DEN1DX(12, CP2, DEN1X)
  CALL DEN1DY(13, CP2, DEN1Y)
  CLOSE(11)
  CLOSE(12)
  CLOSE(13)
!
  !OPEN(21, FILE = 'imag3d-den2d_xz.txt')
  !OPEN(22, FILE = 'imag3d-den2d_xy.txt')
  !OPEN(23, FILE = 'imag3d-den2d_yz.txt')
  !CALL DEN2DXZ(21, CP2, DEN2XZ)
  !CALL DEN2DXY(22, CP2, DEN2XY)
  !CALL DEN2DYZ(23, CP2, DEN2YZ)
  !CLOSE(21)
  !CLOSE(22)
  !CLOSE(23)

  !OPEN(32, FILE = 'imag3d-den3d_x0z.txt')
  !OPEN(33, FILE = 'imag3d-den3d_xy0.txt')
  !CALL WRITE_2DXZ(32, CP2)
  !CALL WRITE_2DXY(33, CP2)
  !CLOSE(32)
  !CLOSE(33)
!  
  CALL SYSTEM_CLOCK (CLCK_COUNTS_END, CLCK_RATE)
  CALL CPU_TIME(T2)
  WRITE (7, 1001)  
  WRITE (4, 1001)
  CLOSE (4)
  WRITE (7,*)
  WRITE (7,'(A,I7,A)') ' Clock Time: ', (CLCK_COUNTS_END - CLCK_COUNTS_BEG)/INT (CLCK_RATE,8), ' seconds'
  WRITE (7,'(A,I7,A)') '   CPU Time: ', INT(T2-T1), ' seconds' 
  CLOSE (7)
!
END PROGRAM GROSS_PITAEVSKII_SSCN_3D

SUBROUTINE INITIALIZE()
!   Routine that initializes the constant and variables
!   Calculates the initial wave function CP
  USE COMM_DATA, ONLY : NX, NY, NZ, NX2, NY2, NZ2, PI
  USE GPE_DATA, ONLY : DX, DY, DZ, ANU, LAMBDA, GAMMA, &
        X, Y, Z, X2, Y2, Z2, R2, CP
  IMPLICIT NONE
  REAL (8), DIMENSION(0:NX, 0:NY, 0:NZ) :: TMP3D
  REAL (8) :: SP, PI34	 
  INTEGER :: I, J, K
!
  PI34 = SQRT(PI*SQRT(PI)) ! PI^(3/4)
  SP = PI34 / SQRT(SQRT(ANU * LAMBDA * GAMMA))
!
  FORALL (I=0:NX) X(I) = (I-NX2)*DX
  FORALL (J=0:NY) Y(J) = (J-NY2)*DY
  FORALL (K=0:NZ) Z(K) = (K-NZ2)*DZ
  X2 = X*X
  Y2 = Y*Y
  Z2 = Z*Z
!
  FORALL (J=0:NY, K=0:NZ)
    TMP3D(0:NX,J,K) =  (GAMMA * X2 + ANU * Y2(J) + LAMBDA * Z2(K))
    R2(0:NX,J,K) = X2 + Y2(J) + Z2(K)
  END FORALL
  CP = EXP(-TMP3D/2.D0)/SP
END SUBROUTINE INITIALIZE

SUBROUTINE CALCULATE_TRAP()
! Routine that  initializes the harmonic oscillator potential V.
  USE COMM_DATA, ONLY : NX, NY, NZ
  USE GPE_DATA, ONLY : XOP, GAMMA, ANU, LAMBDA, V, X2, Y2, Z2
  IMPLICIT NONE
  INTEGER :: J, K
  REAL (8) :: ANU2, LAM2, GAM2
!
  GAM2 = GAMMA * GAMMA
  ANU2 = ANU * ANU
  LAM2 = LAMBDA * LAMBDA 
  FORALL (J=0:NY, K=0:NZ)
    V(0:NX,J,K) = (GAM2 * X2 + ANU2 * Y2(J) + LAM2 * Z2(K))
  END FORALL
END SUBROUTINE CALCULATE_TRAP

SUBROUTINE COEF()
!  Calculates the coefficients needed in subroutines LUX, LUY, LUZ.      
  USE COMM_DATA, ONLY : NXX, NYY, NZZ
  USE GPE_DATA, ONLY : DX, DY, DZ, DT
  USE CN_DATA
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL (8) :: DX2, DY2, DZ2
  REAL (8) :: DXX, DYY, DZZ
!
  DX2 = DX*DX   ! generating the coefficients for the
  DY2 = DY*DY   ! c-n method (to solve the spatial part)
  DZ2 = DZ*DZ
!     
  DXX = 1.D0/DX2
  DYY = 1.D0/DY2
  DZZ = 1.D0/DZ2
!
  CA0 = 1.D0 + DT*DXX
  CA0R = 1.D0 - DT*DXX
  CB0 = 1.D0 + DT*DYY
  CB0R = 1.D0 - DT*DYY
  CC0 = 1.D0 + DT*DZZ
  CC0R = 1.D0 - DT*DZZ
!
  CT0X = -DT*DXX/2.D0
  CALA(NXX) = 0.D0
  CGAA(NXX) = -1.D0/CA0
  DO I = NXX, 1, -1
     CALA(I-1) = CT0X*CGAA(I)
     CGAA(I-1) = -1.D0/(CA0+CT0X*CALA(I-1))
  END DO
!
  CT0Y = -DT*DYY/2.D0
  CALB(NYY) = 0.D0
  CGAB(NYY) = -1.D0/CB0
  DO J = NYY, 1, -1
     CALB(J-1) = CT0Y*CGAB(J)
     CGAB(J-1) = -1.D0/(CB0+CT0Y*CALB(J-1))
  END DO
!
  CT0Z = -DT*DZZ/2.D0
  CALC(NZZ) = 0.D0
  CGAC(NZZ) = -1.D0/CC0
  DO K = NZZ, 1, -1
     CALC(K-1) = CT0Z*CGAC(K)
     CGAC(K-1) = -1.D0/(CC0+CT0Z*CALC(K-1))
  END DO
END SUBROUTINE COEF

SUBROUTINE CALCNU(CP, DT) ! Exact solution
!  Solves the partial differential equation with the potential and the
!  nonlinear term.      
  USE COMM_DATA, ONLY : NX, NY, NZ, NXX, NYY, NZZ
  USE GPE_DATA, ONLY : V, G
  IMPLICIT NONE
!-------------------------------------------------
  REAL (8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(IN) :: DT
  REAL (8), DIMENSION(0:NX,0:NY,0:NZ) :: P2, TMP3D, DP
!
  P2 = CP * CP
  TMP3D = V + G*P2
!
  CP = CP*EXP(-DT*TMP3D)
END SUBROUTINE CALCNU

SUBROUTINE LUX(CP)
!  Solves the partial differential equation only with the X-space
!  derivative term using the Crank-Nicholson method 
  USE COMM_DATA
  USE CN_DATA, ONLY : CA0R, CT0X, CALA, CGAA  
  IMPLICIT NONE
  REAL (8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: CP
  INTEGER :: I,J,K
  REAL (8), DIMENSION(0:NX,0:NY,0:NZ) :: CBE,CXX
!  REAL (8), DIMENSION(0:NY,0:NZ) :: CXX
  
  DO J=0,NY
   DO K=0,NZ 
    CBE(NXX,J,K) = CP(NXX,J,K)
    DO I = NXX, 1, -1
       CXX(I,J,K) = -CT0X*CP(I+1,J,K)+CA0R*CP(I,J,K)-CT0X*CP(I-1,J,K)
    END DO
   END DO
  END DO

 DO J=0,NY
  DO K=0,NZ 
   DO I = NXX, 1, -1
     CBE(I-1,J,K) = CGAA(I)*(CT0X*CBE(I,J,K)-CXX(I,J,K))
   END DO
!-----------------------------
! Boundary condition reflecting:
   CP(0,J,K) = 0.D0
   DO I = 0, NXX
      CP(I+1,J,K) = CALA(I)*CP(I,J,K) + CBE(I,J,K)
   END DO
   CP(NX,J,K) = 0.D0
  END DO
 END DO
!-----------------------------
END SUBROUTINE LUX

SUBROUTINE LUY(CP)
!  Solves the partial differential equation only with the Y-space
!  derivative term using the Crank-Nicholson method
  USE COMM_DATA
  USE CN_DATA, ONLY :  CB0R, CT0Y, CALB, CGAB 
  IMPLICIT NONE
  REAL (8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: CP
  INTEGER :: J,I,K
  REAL (8) :: CBE(0:NX,0:NY,0:NZ), CYY(0:NX,0:NY,0:NZ)

  DO I=0,NX
   DO K=0,NZ 
    CBE(I,NYY,K) = CP(I,NYY,K)
    DO J = NYY, 1, -1
       CYY(I,J,K) = -CT0Y*CP(I,J+1,K)+CB0R*CP(I,J,K)-CT0Y*CP(I,J-1,K)
    END DO
   END DO
  END DO

  DO I=0,NX
   DO K=0,NZ 
    DO J = NYY, 1, -1
       CBE(I,J-1,K) = CGAB(J)*(CT0Y*CBE(I,J,K)-CYY(I,J,K))
    END DO
!-----------------------------
! Boundary condition reflecting:
    CP(I,0,K) = 0.D0
    DO J = 0, NYY
       CP(I,J+1,K) = CALB(J)*CP(I,J,K) + CBE(I,J,K)
    END DO
    CP(I,NY,K) = 0.D0
   END DO
  END DO
!-----------------------------
END SUBROUTINE LUY

SUBROUTINE LUZ(CP)
!  Solves the partial differential equation only with the Y-space
!  derivative term using the Crank-Nicholson method
  USE COMM_DATA
  USE CN_DATA, ONLY :   CC0R, CT0Z, CALC, CGAC 
  IMPLICIT NONE
  REAL (8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: CP
  INTEGER :: K,I,J
  REAL (8) :: CBE(0:NX,0:NY,0:NZ), CZZ(0:NX,0:NY,0:NZ)

  DO I=0,NX
   DO  J=0,NY 
    CBE(I,J,NZZ) = CP(I,J,NZZ)
    DO K = NZZ, 1, -1
       CZZ(I,J,K) = -CT0Z*CP(I,J,K+1)+CC0R*CP(I,J,K)-CT0Z*CP(I,J,K-1)
    END DO
   END DO
  END DO

  DO I=0,NX
   DO J=0,NY 
    DO K = NZZ, 1, -1
       CBE(I,J,K-1) = CGAC(K)*(CT0Z*CBE(I,J,K)-CZZ(I,J,K))
    END DO
!-----------------------------
! Boundary condition reflecting:
    CP(I,J,0) = 0.D0
    DO K = 0, NZZ
       CP(I,J,K+1) = CALC(K)*CP(I,J,K) + CBE(I,J,K)
    END DO
    CP(I,J,NZ) = 0.D0
   END DO
  END DO
!-----------------------------
END SUBROUTINE LUZ

SUBROUTINE NORM(CP, ZNORM)
!  Calculates the normalization of the wave function and sets it to unity.
  USE COMM_DATA, ONLY : NX, NY, NZ
  USE GPE_DATA, ONLY : DX, DY, DZ
  IMPLICIT NONE
  REAL (8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(OUT) :: ZNORM
!-------------------------------------------------
  INTERFACE
    FUNCTION INTEGRATE(U, DX, DY, DZ) RESULT(VALUE)
      REAL (8), DIMENSION(0:, 0:, 0:), INTENT(IN) :: U
      REAL (8), INTENT (IN) :: DX, DY, DZ
      REAL (8) :: VALUE
    END FUNCTION INTEGRATE
  END INTERFACE
!-------------------------------------------------
  REAL (8), DIMENSION(0:NX, 0:NY, 0:NZ) :: P2
!
  P2 = CP*CP
  ZNORM = SQRT(INTEGRATE(P2, DX, DY, DZ))
  CP = CP/ZNORM
END SUBROUTINE NORM

SUBROUTINE RAD(CP2, R)
!  Calculates the root mean square size RMS
  USE COMM_DATA, ONLY : NX, NY, NZ
  USE GPE_DATA, ONLY : DX, DY, DZ, X2, Y2, Z2, R2
  IMPLICIT NONE
  REAL (8), DIMENSION(0:,0:,0:), INTENT(IN) :: CP2
  REAL (8), DIMENSION(:), INTENT(OUT) :: R
!-------------------------------------------------
  INTERFACE
    FUNCTION INTEGRATE(U, DX, DY, DZ) RESULT(VALUE)
      REAL (8), DIMENSION(0:, 0:, 0:), INTENT(IN) :: U
      REAL (8), INTENT (IN) :: DX, DY, DZ
      REAL (8) :: VALUE
    END FUNCTION INTEGRATE
  END INTERFACE
!-------------------------------------------------
  INTEGER :: I, J, K
  REAL (8), DIMENSION(0:NX, 0:NY, 0:NZ) :: TMP3D
!
  TMP3D = R2*CP2
   R(1) = SQRT(INTEGRATE(TMP3D, DX, DY, DZ))
  FORALL(J=0:NY, K=0:NZ) TMP3D(:,J,K) = X2*CP2(:,J,K)
   R(2) = SQRT(INTEGRATE(TMP3D, DX, DY, DZ))
  FORALL(I=0:NX, K=0:NZ) TMP3D(I,:,K) = Y2*CP2(I,:,K)
   R(3) = SQRT(INTEGRATE(TMP3D, DX, DY, DZ))
  FORALL(I=0:NX, J=0:NY) TMP3D(I,J,:) = Z2*CP2(I,J,:)
   R(4) = SQRT(INTEGRATE(TMP3D, DX, DY, DZ))
END SUBROUTINE RAD

 SUBROUTINE CHEM(CP, MU, EN)
!  Calculates the chemical potential MU and energy EN.  CP is the wave
!  function, V is the potential and G is the nonlinearity.
   USE COMM_DATA, ONLY : NX, NY, NZ, NXX, NYY, NZZ
   USE GPE_DATA, ONLY : DX, DY, DZ, X2, V, G   
   IMPLICIT NONE
   REAL (8), DIMENSION(0:, 0:, 0:), INTENT(IN) :: CP
   REAL (8), INTENT(OUT) :: MU, EN
!-------------------------------------------------
  INTERFACE
    PURE FUNCTION DIFF(P, DX) RESULT (DP)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: P
      REAL (8), INTENT(IN) :: DX
      REAL (8), DIMENSION(0:SIZE(P)-1) :: DP
    END FUNCTION DIFF
  END INTERFACE
!-------------------------------------------------
  INTERFACE
    FUNCTION INTEGRATE(U, DX, DY, DZ) RESULT(VALUE)
      REAL (8), DIMENSION(0:, 0:, 0:), INTENT(IN) :: U
      REAL (8), INTENT (IN) :: DX, DY, DZ
      REAL (8) :: VALUE
    END FUNCTION INTEGRATE
  END INTERFACE
!-------------------------------------------------
   INTEGER :: I, J, K
   REAL (8), DIMENSION(0:NX, 0:NY, 0:NZ) :: DP, DPX, DPY, DPZ
   REAL (8), DIMENSION(0:NX, 0:NY, 0:NZ) :: P2, GP2,  DP2, TMP3D, EMP3D
!
   DO I = 0, NX; DO J = 0,NY
      DPZ(I,J,0:NZ) = DIFF(CP(I,J,0:NZ), DZ)
   END DO; END DO
   DO I = 0, NX; DO K = 0,NZ
      DPY(I,0:NY,K) = DIFF(CP(I,0:NY,K), DY)
   END DO; END DO
   DO J = 0, NY; DO K = 0,NZ
      DPX(0:NX,J,K) = DIFF(CP(0:NX,J,K), DX)
   END DO; END DO
!
   P2 =  CP * CP
   GP2 = G * P2
   DP2 = DPX*DPX + DPY*DPY + DPZ*DPZ
!
   TMP3D = (V + GP2)*P2 + DP2
   EMP3D = (V + GP2/2.D0)*P2 + DP2
!
   MU = INTEGRATE(TMP3D, DX, DY, DZ)
   EN = INTEGRATE(EMP3D, DX, DY, DZ)
END SUBROUTINE CHEM

FUNCTION INTEGRATE(U, DX, DY, DZ) RESULT(RES)
  IMPLICIT NONE
  REAL (8), DIMENSION(0:, 0:, 0:), INTENT(IN) :: U
  REAL (8), INTENT (IN) :: DX, DY, DZ
  REAL (8) :: RES
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
  REAL (8), DIMENSION(0:SIZE(U,1)-1, 0:SIZE(U,2)-1) :: TMP2D
  REAL (8), DIMENSION(0:SIZE(U,1)-1) :: TMP1D
  INTEGER :: I, J, NX, NY 
!
  NX = SIZE(U,1)-1
  NY = SIZE(U,2)-1
  FORALL (I = 0:NX, J = 0:NY) TMP2D(I,J) = SIMP(U(I,J,0:), DZ)
  FORALL (I = 0:NX) TMP1D(I) = SIMP(TMP2D(I,0:), DY)
  RES = SIMP(TMP1D, DX)
END FUNCTION INTEGRATE

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

PURE FUNCTION DIFF(P,DX) RESULT (DP)
! Computes the first derivative DP of P using Richardsonextrapolation formula. 
! The derivative at the boundaries are assumed to be zero.
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
 
SUBROUTINE WRITE_3D(FUNIT, U2)
  USE COMM_DATA, ONLY : NX, NY, NZ
  INTEGER :: FUNIT
  REAL (8), DIMENSION(0:,0:,0:) :: U2
  INTEGER :: I, J, K
!
  DO I = 0, NX
     DO J = 0, NY
        DO K = 0, NZ
           WRITE(FUNIT, 1000) U2(I,J,K)
        END DO
     END DO
  END DO
   1000 FORMAT(E17.6E3)  
END SUBROUTINE WRITE_3D
 
SUBROUTINE WRITE_2DXY(FUNIT, U2)
  USE COMM_DATA, ONLY : NX, NY, NZ2
  USE GPE_DATA, ONLY : X, Y
  INTEGER :: FUNIT
  REAL (8), DIMENSION(0:,0:,0:) :: U2
  INTEGER :: I, J, K
!
  K = NZ2
  DO I = 0, NX
     DO J = 0, NY
        WRITE(FUNIT, 1000) X(I), Y(J), U2(I,J,K)   
     END DO
     WRITE(FUNIT, *)
  END DO
  1000 FORMAT(2F10.2, E17.6E3)
END SUBROUTINE WRITE_2DXY
 
SUBROUTINE WRITE_2DXZ(FUNIT, U2)
  USE COMM_DATA, ONLY : NX, NZ, NY2
  USE GPE_DATA, ONLY : X, Z
  INTEGER :: FUNIT
  REAL (8), DIMENSION(0:,0:,0:) :: U2
  INTEGER :: I, J, K
!
  J = NY2
  DO I = 0, NX
     DO K = 0, NZ
        WRITE(FUNIT, 1000) X(I), Z(K), U2(I,J,K)
     END DO
     WRITE(FUNIT, *)
  END DO
  1000 FORMAT(2F10.2, E17.6E3)
END SUBROUTINE WRITE_2DXZ

SUBROUTINE DEN2DXY(FUNIT, U2, DE2DXY)
  USE COMM_DATA, ONLY : NX, NY  
  USE GPE_DATA, ONLY : X,  DZ, Y
  IMPLICIT NONE
  REAL (8), DIMENSION(0:, 0:, 0:), INTENT(IN) :: U2
  REAL (8), DIMENSION(0:, 0:), INTENT(OUT) :: DE2DXY
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
  INTEGER :: I, J, FUNIT	
! 
  FORALL (I = 0:NX, J = 0:NY) DE2DXY(I,J) = SIMP(U2(I,J,0:), DZ)
  DO I = 0, NX
   DO J = 0,NY
     WRITE(FUNIT, 1000) X(I),Y(J), DE2DXY(I,J)
   END DO
    WRITE(FUNIT, 1000)
  END DO 
  1000 FORMAT(2F10.2, E17.6E3)
END SUBROUTINE DEN2DXY

SUBROUTINE DEN2DXZ(FUNIT, U2, DE2DXZ)
  USE COMM_DATA, ONLY : NX, NZ  
  USE GPE_DATA, ONLY : X, DY, Z
  IMPLICIT NONE
  REAL (8), DIMENSION(0:, 0:, 0:), INTENT(IN) :: U2
  REAL (8), DIMENSION(0:, 0:), INTENT(OUT) :: DE2DXZ
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
  INTEGER :: I, K, FUNIT
!	 
  FORALL (I = 0:NX, K = 0:NZ) DE2DXZ(I,K) = SIMP(U2(I,0:,K), DY)
  DO I = 0, NX
    DO K = 0,NZ
      WRITE(FUNIT, 1000) X(I),Z(K), DE2DXZ(I,K)
    END DO
     WRITE(FUNIT, 1000)
  END DO
  1000 FORMAT(2F10.2, E17.6E3) 
END SUBROUTINE DEN2DXZ

SUBROUTINE DEN2DYZ(FUNIT, U2, DE2DYZ)
  USE COMM_DATA, ONLY : NY, NZ  
  USE GPE_DATA, ONLY : Y, DX, Z
  IMPLICIT NONE
  REAL (8), DIMENSION(0:, 0:, 0:), INTENT(IN) :: U2
  REAL (8), DIMENSION(0:, 0:), INTENT(OUT) :: DE2DYZ
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
  INTEGER :: J, K, FUNIT
!
  FORALL (J = 0:NY, K = 0:NZ) DE2DYZ(J,K) = SIMP(U2(0:,J,K), DX)
  DO J = 0, NY
  DO K = 0,NZ
        WRITE(FUNIT, 1000) Y(J),Z(K), DE2DYZ(J,K)
  END DO
      WRITE(FUNIT, 1000)
  END DO
  1000 FORMAT(2F10.2, E17.6E3) 
END SUBROUTINE DEN2DYZ

SUBROUTINE DEN1DX(FUNIT, U2, DE1DX)
  USE COMM_DATA, ONLY : NX, NZ
  USE GPE_DATA, ONLY : X, DY, DZ 
  IMPLICIT NONE
  REAL (8), DIMENSION(0:, 0:, 0:), INTENT(IN) :: U2
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
  REAL (8), DIMENSION(0:NX, 0:NZ) :: TMP2D
  INTEGER :: I, K, FUNIT	 
!
  FORALL (I = 0:NX, K = 0:NZ) TMP2D(I,K) = SIMP(U2(I,0:,K), DY)
      DO I = 0, NX 
	 DE1DX(I) = SIMP(TMP2D(I,0:), DZ)
         WRITE(FUNIT, 1001) X(I), DE1DX(I)
      END DO
  1001 FORMAT(F10.2, E17.6E3)
END SUBROUTINE DEN1DX

SUBROUTINE DEN1DY(FUNIT, U2, DE1DY)
  USE COMM_DATA, ONLY : NY, NZ
  USE GPE_DATA, ONLY : Y, DX, DZ  
  IMPLICIT NONE
  REAL (8), DIMENSION(0:, 0:, 0:), INTENT(IN) :: U2
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
  REAL (8), DIMENSION(0:NY, 0:NZ) :: TMP2D
  INTEGER :: J, K, FUNIT
!
  FORALL (J = 0:NY, K = 0:NZ) TMP2D(J,K) = SIMP(U2(0:,J,K), DX)
      DO J = 0, NY 
	 DE1DY(J) = SIMP(TMP2D(J,0:), DZ)
         WRITE(FUNIT, 1001) Y(J), DE1DY(J)
      END DO
  1001 FORMAT(F10.2, E17.6E3) 
END SUBROUTINE DEN1DY

SUBROUTINE DEN1DZ(FUNIT, U2, DE1DZ)
  USE COMM_DATA, ONLY : NY, NZ
  USE GPE_DATA, ONLY : Z, DX, DY   
  IMPLICIT NONE
  REAL (8), DIMENSION(0:, 0:, 0:), INTENT(IN) :: U2
  REAL (8), DIMENSION(0:), INTENT(OUT) :: DE1DZ
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
  REAL (8), DIMENSION(0:NY, 0:NZ) :: TMP2D
  INTEGER :: J, K, FUNIT
!
  FORALL (J = 0:NY, K = 0:NZ) TMP2D(J,K) = SIMP(U2(0:,J,K), DX)
      DO K = 0, NZ 
	 DE1DZ(K) = SIMP(TMP2D(0:,K), DY)
         WRITE(FUNIT, 1001) Z(K), DE1DZ(K)
      END DO
  1001 FORMAT(F10.2, E17.6E3)
END SUBROUTINE DEN1DZ

!# File name : imag3d.f90
