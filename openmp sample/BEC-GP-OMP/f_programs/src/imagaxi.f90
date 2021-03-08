!#  File name : imagaxi.f90	!# Last modified : 14 Mar 2016
!#  Fortran program for Gross-Pitaevskii equation in three-dimensional 
!#  axially-symmetric trap by imaginary time propagation (Fortran 90/95 Version)
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
  INTEGER, PARAMETER :: NX = 300, NXX = NX-1 ! No. of X steps
  INTEGER, PARAMETER :: NY = 300, NYY = NY-1, NY2 = NY/2 ! No. of Y steps
! NSTP : Number of iterations to introduce the nonlinearity. 
! NPAS : Number of subsequent iterations with fixed nonlinearity.
! NRUN : Number of final time iterations with fixed nonlinearity.
  INTEGER, PARAMETER :: NSTP = 0, NPAS = 100000 , NRUN = 20000
END MODULE COMM_DATA

MODULE GPE_DATA
  USE COMM_DATA, ONLY : NX, NY
  REAL (8), PARAMETER :: PI = 3.14159265358979D0
  REAL (8), PARAMETER :: AHO = 1.D-6     		 ! Unit of length (l= 1 MICRON)            
  REAL (8), PARAMETER :: Bohr_a0 =  5.2917720859D-11/AHO ! Bohr radius (scaled with AHO)
! 
  REAL (8), PARAMETER :: DX = 0.02D0, DY = 0.02D0, DT = 0.00004D0	! DX, DY : Space step and DT : Time step
  INTEGER, PARAMETER  :: NATOMS = 400			 ! Number of Atoms
  REAL (8), PARAMETER :: AS = 70.71602D0*Bohr_a0	 ! Scattering length (in units of Bohr_a0)
  REAL (8), PARAMETER :: NU = 1.0D0, LAMBDA = 4.0D0	 ! NU, LAMBDA : Parameteres of Trap
!   
  REAL (8), PARAMETER :: G0 = 4.D0*PI*AS*NATOMS !18.81D0 ! Three-dimensional nonlinearity 
!
! OPTION and XOP decides which equation to be solved. 
  ! OPTION = 1 Solves -psi_xx-psi_yy+V(x,y)psi+G0|psi|^2 psi=i psi_t
  ! OPTION = 2 Solves [-psi_xx-psi_yy+V(x,y)psi]/2+G0|psi|^2 psi=i psi_t
  ! OPTION = 3 Solves -psi_xx-psi_yy+V(x,y)psi/4+G0|psi|^2 psi=i psi_t
  INTEGER, PARAMETER :: OPTION = 2
  REAL (8) :: G, XOP, TPI
! X(0:N), Y(0:N) : Space mesh, V(0:N) : Potential, and 
! CP(0:N) : Wave function 
  REAL (8), DIMENSION(0:NX, 0:NY) :: V, CP
  REAL (8), DIMENSION(0:NX) :: X, X2
  REAL (8), DIMENSION(0:NY) :: Y, Y2
END MODULE GPE_DATA

MODULE CN_DATA
  USE COMM_DATA, ONLY : NX, NY
  REAL (8), DIMENSION(0:NX) :: CALX, CGAX, CAIPX
  REAL (8), DIMENSION(1:NX) :: CPX
  REAL (8), DIMENSION(0:NY) :: CALY, CGAY
  REAL (8) :: CAIPMX, CAIPMY
  REAL (8) :: CAI0, CAIM
END MODULE CN_DATA 

PROGRAM GROSS_PITAEVSKII_SSCN_AXI
  USE COMM_DATA, ONLY : NX, NY, NY2, NSTP, NPAS, NRUN
  USE GPE_DATA
  IMPLICIT NONE
! Subroutine INITIALIZE() used to initialize the space mesh X(I), 
! potential V(I) and the initial wave function. Subroutine COEF() used to
! generate the coefficients for the Crank-Nicholson Scheme. The routine
! CALCNU() performs time progation for the non-derivative part and LU() performs
! time propagation of derivative part. NORM() calculates the norm and 
! normalizes the wave function, CHEM() and RAD() are used to calculate the 
! chemical potential, energy and the rms radius, respectively. The function 
! DIFF() used to calculate the space derivatives of the wave function used 
! in CHEM() and SIMP() does the integration by Simpson's rule.
!------------------------ INTERFACE BLOCKS -----------------------
  INTERFACE 
    SUBROUTINE INITIALIZE()
      IMPLICIT NONE
    END SUBROUTINE INITIALIZE
  END INTERFACE
!
  INTERFACE 
    SUBROUTINE COEF()
      IMPLICIT NONE
    END SUBROUTINE COEF
  END INTERFACE
!
  INTERFACE 
    SUBROUTINE CALCNU(CP, DT)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
      REAL (8), INTENT(IN) :: DT
    END SUBROUTINE CALCNU
  END INTERFACE
!
  INTERFACE 
    SUBROUTINE LUX(CP)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
    END SUBROUTINE LUX
  END INTERFACE
!
  INTERFACE 
    SUBROUTINE LUY(CP)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
    END SUBROUTINE LUY
  END INTERFACE
!
  INTERFACE
    SUBROUTINE NORM(CP, ZNORM)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
      REAL (8), INTENT(OUT) :: ZNORM
    END  SUBROUTINE NORM
  END INTERFACE
!
  INTERFACE
    SUBROUTINE RAD(CP, RHORMS, ZRMS, RRMS)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
      REAL (8), INTENT(OUT) :: RHORMS, ZRMS, RRMS
    END  SUBROUTINE RAD
  END INTERFACE
!
  INTERFACE
    SUBROUTINE CHEM(CP, MU, EN)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:,0:), INTENT(IN) :: CP
      REAL (8), INTENT(OUT) :: MU, EN
    END  SUBROUTINE CHEM
  END INTERFACE
!------------------------ END INTERFACE BLOCKS -------------------
  INTEGER :: I, J, K
  REAL (8) :: ZNORM, RHORMS, ZRMS, RRMS, MU, EN, TMp
  REAL (8) :: GSTP, T1, T2
 REAL (8), DIMENSION(0:NX,0:NY):: CP2
  INTEGER :: CLCK_COUNTS_BEG, CLCK_COUNTS_END, CLCK_RATE

!-----------------------------------------------------------
  INTERFACE
    PURE FUNCTION SIMP(F, DX) RESULT(VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
!-----------------------------------------------------------
!  
  CALL SYSTEM_CLOCK ( CLCK_COUNTS_BEG, CLCK_RATE )
  CALL CPU_TIME (T1)
!
  SELECT CASE (OPTION)
    CASE (1)
      XOP = 1.0D0
    CASE (2)
      XOP = 2.0D0
    CASE (:0,3:)
      PRINT *, 'ERROR: Wrong option', OPTION
      STOP
  END SELECT
!
! INITIALIZE() initializes the starting normalized wave function 'CP'.
  CALL INITIALIZE()
!
  OPEN(7, FILE = 'imagaxi-out.txt')
  OPEN(4, FILE = 'imagaxi-rms.txt')
! 
  WRITE(7,900) OPTION 
  WRITE(4,900) OPTION
  WRITE(7,*) 
  WRITE(4,*)
  WRITE(7,901) NATOMS, AHO
  WRITE(7,902) AS/Bohr_a0 
  WRITE(7,903) G0  
  WRITE(7,904) NU, LAMBDA
  WRITE(7,*)
  WRITE(7,905) NX, NY
  WRITE(7,906) DX, DY 
  WRITE(7,907) NSTP, NPAS, NRUN  
  WRITE(7,908) DT 
  WRITE(7,*)
  900 FORMAT(' Imaginary time propagation axially-symmetric trap,   OPTION = ',I3 )
  901 FORMAT('  Number of Atoms N =',I10,', Unit of length AHO =',F12.8,' m')
  902 FORMAT('  Scattering length a = ',F9.2,'*a0 ')
  903 FORMAT('  Nonlinearity G_3D =',F16.7 )
  904 FORMAT('  Parameters of trap: NU =',F7.2, ', LAMBDA =',F7.2)
  905 FORMAT(' # Space Stp: NX = ', I8, ', NY = ', I8)
  906 FORMAT('  Space Step: DX = ', F10.6, ', DY = ', F10.6)
  907 FORMAT(' # Time Stp : NSTP = ',I9,', NPAS = ',I9,', NRUN = ',I9)
  908 FORMAT('   Time Step:   DT = ', F10.6)
! 
! COEF() defines the coefficients of Crank-Nicholson Scheme.
  CALL COEF()
! NORM() calculates norm and restores normalization.
  CALL NORM(CP, ZNORM)
! CHEM() calculates the chemical potential MU and energy EN.
  CALL CHEM(CP, MU, EN)
! RAD() calculates the r.m.s radius RMS.
  CALL RAD(CP, RHORMS, ZRMS, RRMS)
! 
  WRITE (7, 1001)
  WRITE (7, 1002)
  WRITE (7, 1001)
  WRITE (7, 1003) ZNORM, MU/XOP, EN/XOP, RRMS, CP(0,NY2+1)**2
  WRITE (4, 1001)
  WRITE (4, 1012)
  WRITE (4, 1001)
  WRITE (4, 1013) RRMS, RHORMS, ZRMS
  1001 FORMAT(19X,53('-'))
  1002 FORMAT(20X,'Norm',6X,'Chem',8X,'Ener/N',6X,'<r>',3X,'|Psi(0,0)|^2')
  1003 FORMAT ('Initial : ', 4X, F11.4, 2F12.5, 2F11.5)
  1012 FORMAT ('Values of rms size:', 8X, '<r>', 13X, '<rho>', 13X, '<z>')
  1013 FORMAT (11X, 'Initial:', F14.5, 2F16.5)
! 
!   OPEN(1, FILE = 'imagaxi-den-init.txt')
!     DO I = 0, NX, 10 ! Writes initial density
!      DO J = 0, NY, 10
!         WRITE (1,999) X(I), Y(J), (CP(I,J))**2
!      END DO
!      WRITE (1,*)
!     END DO
!   CLOSE(1)
!   CP2 = CP * CP
!   OPEN(12, FILE = 'imagaxi-den-init1d_rho.txt')
!    DO I = 0, NX 
! 	 TMP = SIMP(CP2(I,:), DY)
!          WRITE(12, 99) X(I), TMP
!    END DO
!   CLOSE(12)
! 
!   OPEN(13, FILE = 'imagaxi-den-init1d_z.txt')
!    DO I = 0, NY 
! 	 TMP =2.d0*PI* SIMP(CP2(:,I)*X(:), DX)
!          WRITE(13, 99) Y(I), TMP
!    END DO
!   CLOSE(13)
!
  IF(NSTP /= 0)THEN
    GSTP =  XOP*G0/DFLOAT(NSTP)
    G = 0.D0
    DO K = 1, NSTP ! NSTP iterations to introduce the nonlinearity
      G = G + GSTP
      CALL CALCNU(CP, DT)	! CALCNU() performs time propagation with non-derivative parts.
      CALL LUX(CP)		! LU() performs the time iteration with space derivative alone.
      CALL LUY(CP)
      CALL NORM(CP,ZNORM)
    END DO
    CALL CHEM(CP, MU, EN)
    CALL RAD(CP, RHORMS, ZRMS, RRMS)
    WRITE (7, 1004) ZNORM, MU/XOP, EN/XOP, RRMS, CP(0,NY2+1)**2
    WRITE (4, 1014) RRMS, RHORMS, ZRMS
    1004 FORMAT('After NSTP iter.:', F8.4, 2F12.5, 2F11.5)
    1014 FORMAT (2X, 'After NSTP iter.:', F14.5, 2F16.5)
!
!   OPEN(2, FILE = 'imagaxi-den-nstp.txt')
!     DO I = 0, NX, 10 ! Writes density after NSTP iterations 
!      DO J = 0, NY, 10
!         WRITE (2,999) X(I), Y(J), (CP(I,J))**2
!      END DO
!      WRITE (2,*)
!     END DO
!   CLOSE(2)  
!   CP2 = CP * CP
!   OPEN(22, FILE = 'imagaxi-den-nstp1d_rho.txt')
!    DO I = 0, NX 
! 	 TMP = SIMP(CP2(I,:), DY)
!          WRITE(22, 99) X(I), TMP
!    END DO
!   CLOSE(22)
! 
!   OPEN(23, FILE = 'imagaxi-den-nstp1d_z.txt')
!    DO I = 0, NY 
! 	 TMP =2.d0*PI* SIMP(CP2(:,I)*X(:), DX)
!          WRITE(23, 99) Y(I), TMP
!    END DO
!   CLOSE(23)
!
  ELSE
     G = XOP*G0
  END IF
!
  DO K = 1, NPAS ! NPAS Iterations transient
     CALL CALCNU(CP, DT)
     CALL LUX(CP)
     CALL LUY(CP)
     CALL NORM(CP,ZNORM)
  END DO
  IF(NPAS /= 0)THEN 
    CALL CHEM(CP, MU, EN)
    CALL RAD(CP, RHORMS, ZRMS, RRMS)
    WRITE(7,1005) ZNORM, MU/XOP, EN/XOP, RRMS, CP(0,NY2+1)**2
    WRITE (4, 1015)  RRMS, RHORMS, ZRMS
    1005 FORMAT('After NPAS iter.:',F8.4, 2F12.5, 2F11.5)
    1015 FORMAT (2X, 'After NPAS iter.:', F14.5, 2F16.5)
! 
!     OPEN(3, FILE = 'imagaxi-den-npas.txt')
!       DO I = 0, NX, 10 ! Writes intermediate wave funtion after NPAS iterations
! 	DO J = 0, NY, 10
!         WRITE (3,999) X(I), Y(J), CP(I,J)
! 	END DO
! 	WRITE (3,*)
!       END DO
!     CLOSE(3)
!   CP2 = CP * CP
!   OPEN(32, FILE = 'imagaxi-den-npas1d_rho.txt')
!    DO I = 0, NX 
! 	 TMP = SIMP(CP2(I,:), DY)
!          WRITE(32, 99) X(I), TMP
!    END DO
!   CLOSE(32)
! 
!   OPEN(33, FILE = 'imagaxi-den-npas1d_z.txt')
!    DO I = 0, NY 
! 	 TMP =2.d0*PI* SIMP(CP2(:,I)*X(:), DX)
!          WRITE(33, 99) Y(I), TMP
!    END DO
!   CLOSE(33)

!
  END IF

  DO K = 1, NRUN 
     CALL CALCNU(CP, DT)
     CALL LUX(CP)
     CALL LUY(CP)
     CALL NORM(CP, ZNORM)
  END DO
  IF(NRUN /= 0)THEN 
    CALL CHEM(CP, MU, EN)
    CALL RAD(CP, RHORMS, ZRMS, RRMS)
    WRITE(7,1006) ZNORM, MU/XOP, EN/XOP, RRMS, CP(0,NY2+1)**2
    WRITE(4, 1016) RRMS, RHORMS, ZRMS
    1006 FORMAT('After NRUN iter.:',F8.4, 2F12.5, 2F11.5)
    1016 FORMAT (2X, 'After NRUN iter.:', F14.5, 2F16.5)
  END IF
!      
  OPEN(5, FILE = 'imagaxi-den.txt') ! Writes final wave funtion
      DO I = 0, NX 
	DO J = 0, NY
	  WRITE(5,999) X(I), Y(J), (CP(I,J))**2
	END DO
	WRITE(5,*)
      END DO 
  CLOSE(5)

  CP2 = CP * CP
  OPEN(55, FILE = 'imagaxi-den1d_rho.txt')
   DO I = 0, NX 
	 TMP = SIMP(CP2(I,:), DY)
         WRITE(55, 99) X(I), TMP
   END DO
  CLOSE(55)

  OPEN(65, FILE = 'imagaxi-den1d_z.txt')
   DO I = 0, NY 
	 TMP =2.d0*PI* SIMP(CP2(:,I)*X(:), DX)
         WRITE(65, 99) Y(I), TMP
   END DO
  CLOSE(65)
! 
   99 FORMAT(F10.3, E16.8)
  CALL SYSTEM_CLOCK (CLCK_COUNTS_END, CLCK_RATE)
  CALL CPU_TIME(T2)
  WRITE (7, 1001)
  WRITE (4, 1001)
  CLOSE (4)
  WRITE (7,*)
  WRITE (7,'(A,I7,A)') ' Clock Time: ', (CLCK_COUNTS_END - CLCK_COUNTS_BEG)/INT (CLCK_RATE,8), ' seconds'
  WRITE (7,'(A,I7,A)') '   CPU Time: ', INT(T2-T1), ' seconds' 
  CLOSE (7)
  999 FORMAT(2F10.3, E16.8)
! 
END PROGRAM GROSS_PITAEVSKII_SSCN_AXI

SUBROUTINE INITIALIZE()
!   Routine that initizlizes the constant and variables.
!   Calculates the potential term V and the initial wave function CP
  USE COMM_DATA, ONLY : NX, NY, NY2
  USE GPE_DATA
  IMPLICIT NONE
  REAL (8), DIMENSION(0:NX,0:NY):: TMP2D
  REAL (8) :: PI3, FAC, NU2, LAM2
  INTEGER :: I
! 
  TPI = 2.0D0*PI
  NU2 = NU * NU
  LAM2 = LAMBDA * LAMBDA
  PI3 = PI**3
  FAC = SQRT(SQRT(LAMBDA*NU2/(PI3)))
  FORALL (I=0:NX) X(I) = I*DX
  FORALL (I=0:NY) Y(I) = (I-NY2)*DY
  X2 = X*X
  Y2 = Y*Y
! 
  FORALL(I=0:NX) V(I,:) = (NU2*X2(I) + LAM2*Y2)
  FORALL(I=0:NX) TMP2D(I,:) = (NU*X2(I) + LAMBDA*Y2)/2.0D0
  CP = FAC*EXP(-TMP2D)
END SUBROUTINE INITIALIZE
!
SUBROUTINE COEF()
!   Routine that initizlizes the constant and variables.
!   Calculates the potential term V and the initial wave function CP
  USE COMM_DATA, ONLY : NX, NY, NXX, NYY
  USE GPE_DATA, ONLY : DX, DY, DT, X
  USE CN_DATA
  IMPLICIT NONE
  INTEGER :: I, J
  REAL (8) :: DX2, DY2
!
  DX2 = DX*DX   ! Generating the Coefficients for the
  DY2 = DY*DY   ! C-N method (to solve the spatial part)           
  CAIPMX = -DT/(2.0D0*DX2)
  CAI0 = 1.D0 + DT/DX2
  CALX(1) = 1.0D0
  CGAX(1) = -1.0D0/(CAI0+CAIPMX+DT/(4.D0*DX*X(1)))
!
  DO J = 1, NXX
     CPX(J) = -DT/(4.0D0*DX*X(J))
     CAIPX(J) = CAIPMX-CPX(J)
     CAIM = CAIPMX+CPX(J)
     CALX(J+1) = CGAX(J)*CAIM
     CGAX(J+1) = -1.0D0/(CAI0+(CAIPMX+DT/(4.D0*DX*X(J+1)))*CALX(J+1))
  END DO
!
  CAIPMY = -DT/(2.0D0*DY2)
  CAI0 = 1.D0 + DT/DY2
  CALY(NYY) = 0.0D0
  CGAY(NYY) = -1.0D0/CAI0
!
  DO J = NYY, 1, -1
     CALY(J-1) = CGAY(J)*CAIPMY
     CGAY(J-1) = -1.0D0/(CAI0+CAIPMY*CALY(J-1))
  END DO
END SUBROUTINE COEF

SUBROUTINE CALCNU(CP, DT) ! Exact solution
! Solves the partial differential equation with the potential and the
! nonlinear term.
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : V, G
  IMPLICIT NONE
  REAL (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(IN) :: DT
  REAL (8), DIMENSION(0:NX,0:NY) :: P2, TMP2D
! 
  P2 = CP*CP
  TMP2D = V + G*P2
  CP = CP*EXP(-DT*TMP2D)
END SUBROUTINE CALCNU

SUBROUTINE LUX(CP)
! Solves the partial differential equation only with the X-space
! derivative term using the Crank-Nicholson method
  USE COMM_DATA, ONLY : NX, NY, NXX
  USE CN_DATA, ONLY : CALX, CGAX, CAIPX, CAIPMX, CPX
  IMPLICIT NONE
  REAL (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
  INTEGER :: I, J
  REAL (8), DIMENSION (0:NX, 0:NY) :: CBE,CXX
!  REAL (8), DIMENSION (0:NY) :: CXX
! 
  DO J = 0, NY
     CBE(1,J) = 0.0D0
     DO I = 1, NXX
        CXX(I,J) = CP(I,J)-CAIPMX*(CP(I+1,J)-2.0D0*CP(I,J)+CP(I-1,J)) &
                - CPX(I)*(CP(I+1,J)-CP(I-1,J))
     END DO
  END DO

  DO J = 0, NY
     DO I = 1, NXX
        CBE(I+1,J) = CGAX(I)*(CAIPX(I)*CBE(I,J)-CXX(I,J))
     END DO
     CP(NX,J) = 0.0D0
     DO I = NX, 1, -1
        CP(I-1,J) = CALX(I)*CP(I,J) + CBE(I,J)
     END DO
  END DO
END SUBROUTINE LUX

SUBROUTINE LUY(CP)
!  Solves the partial differential equation only with the Y-space
!  derivative term using the Crank-Nicholson method
  USE COMM_DATA, ONLY : NX, NY, NYY
  USE CN_DATA, ONLY : CALY, CGAY, CAIPMY
  IMPLICIT NONE
  REAL (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
  INTEGER :: I, J
  REAL (8), DIMENSION (0:NX,0:NY) :: CBE,CYY
!  REAL (8), DIMENSION (NX) :: CYY
! 
  DO I = 0, NX
     CBE(I,NYY) = CP(I, NY)
     DO J = NYY, 1, -1
        CYY(I,J) = CP(I,J)-CAIPMY*(CP(I,J+1)-2.0D0*CP(I,J)+CP(I,J-1))
     END DO
  END DO

  DO I = 0, NX
   DO J = NYY, 1, -1
      CBE(I,J-1) = CGAY(J)*(CAIPMY*CBE(I,J)-CYY(I,J))
   END DO
   DO J = 0, NY-1
      CP(I,J+1) = CALY(J)*CP(I,J)+CBE(I,J)
   END DO
  END DO
END SUBROUTINE LUY

SUBROUTINE NORM(CP, ZNORM)
!  Calculates the normalization of the wave function and sets it to
!  unity.
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : DX, DY, TPI, X
  IMPLICIT NONE
  REAL (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(OUT) :: ZNORM
!-----------------------------------------------------------
  INTERFACE
    PURE FUNCTION SIMP(F, DX) RESULT(VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
!-----------------------------------------------------------
  INTEGER :: I, J
  REAL (8), DIMENSION(0:NX,0:NY) :: TMP2D
  REAL (8), DIMENSION(0:NX) :: TMP1D
! 
  TMP2D = CP*CP
  FORALL (I = 0:NX) TMP1D(I) = X(I)*SIMP(TMP2D(I,:), DY)
  ZNORM = SQRT(TPI*SIMP(TMP1D, DX))
  CP = CP/ZNORM
END SUBROUTINE NORM

SUBROUTINE RAD(CP, RHORMS, ZRMS, RRMS) 
! Calculates the root mean square size RMS
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : DX, DY, TPI, X, X2, Y2
  IMPLICIT NONE
  REAL (8), DIMENSION(0:,0:), INTENT(INOUT) :: CP
  REAL (8), INTENT(OUT) :: RHORMS, ZRMS, RRMS
!-----------------------------------------------------------
  INTERFACE
    PURE FUNCTION SIMP(F, DX) RESULT(VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
!-----------------------------------------------------------
  INTEGER :: I
  REAL (8), DIMENSION(0:NX,0:NY) :: TMP2D, TMPY1, TMPY2
  REAL (8), DIMENSION(0:NX) :: TMP1D
! 
  TMP2D = CP*CP
  FORALL (I = 0:NX) TMPY1(I,:) = X2(I)*TMP2D(I,:)
  FORALL (I = 0:NY) TMPY2(:,I) = Y2(I)*TMP2D(:,I)
  FORALL (I = 0:NX) TMP1D(I) = X(I)*SIMP(TMPY1(I,:), DY)
  RHORMS = SQRT(TPI*SIMP(TMP1D, DX))
  FORALL (I = 0:NX) TMP1D(I) = X(I)*SIMP(TMPY2(I,:), DY)
  ZRMS = SQRT(TPI*SIMP(TMP1D, DX))
  RRMS = SQRT(RHORMS*RHORMS + ZRMS*ZRMS)
END SUBROUTINE RAD

SUBROUTINE CHEM(CP, MU, EN)
  ! Calculates the chemical potential MU and energy EN.  CP is the wave
  ! function, V is the potential and G is the nonlinearity.
  USE COMM_DATA, ONLY : NX, NY
  USE GPE_DATA, ONLY : DX, DY, TPI, X, V, G
  IMPLICIT NONE
  REAL (8), DIMENSION(0:,0:), INTENT(IN) :: CP
  REAL (8), INTENT(OUT) :: MU, EN
!-----------------------------------------------------------
  INTERFACE
    PURE FUNCTION DIFF(P, DX) RESULT (DP)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: P
      REAL (8), INTENT(IN) :: DX
      REAL (8), DIMENSION(0:SIZE(P)-1) :: DP
    END FUNCTION DIFF
  END INTERFACE
!-----------------------------------------------------------
  INTERFACE
    PURE FUNCTION SIMP(F, DX) RESULT(VALUE)
      IMPLICIT NONE
      REAL (8), DIMENSION(0:), INTENT(IN) :: F
      REAL (8), INTENT(IN) :: DX
      REAL (8) :: VALUE
    END FUNCTION SIMP
  END INTERFACE
!-----------------------------------------------------------
  INTEGER :: I
  REAL (8), DIMENSION(0:NX,0:NY) :: DPX, DPY, P2, GP2, DP2, TMP2D
  REAL (8), DIMENSION(0:NX) :: TMP1D
! 
 DO I = 0, NX  
    DPY(I,0:NY) = DIFF(CP(I,0:NY), DY)
 END DO 
 DO I = 0,NY
    DPX(0:NX,I) = DIFF(CP(0:NX,I), DX)
 END DO
! 
  P2 = CP*CP
  GP2 = G*P2
  DP2 = DPX*DPX + DPY*DPY
  TMP2D = (V + GP2)*P2 + DP2
  FORALL (I = 0:NX) TMP1D(I) = X(I)*SIMP(TMP2D(I,:), DY)
  MU = TPI*SIMP(TMP1D, DX)
  TMP2D = (V + GP2/2.0D0)*P2 + DP2
  FORALL (I = 0:NX) TMP1D(I) = X(I)*SIMP(TMP2D(I,:), DY)
  EN = TPI*SIMP(TMP1D, DX)
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
  VALUE = DX*(F(0) + 4.0D0*F1 + 2.0D0*F2 + F(N))/3.0D0
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
  DP(0) = 0.0D0
  DP(1) = (P(2) - P(0))/(2.0D0*DX)
  FORALL(I=2:N-2)
    DP(I) = (P(I-2)-8.0D0*P(I-1)+8.0D0*P(I+1)-P(I+2))/(12.0D0*DX)
  END FORALL
  DP(N-1) = (P(N) - P(N-2))/(2.0D0*DX)
  DP(N) = 0.0D0
END FUNCTION DIFF

!# File name : imagaxi.f90
