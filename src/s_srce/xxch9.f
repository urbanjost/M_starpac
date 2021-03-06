*XXCH9
      SUBROUTINE XXCH9(LDSTAK)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     TEST SUBPROGRAM FOR SIMPLE TEST OF
C     THE NONLINEAR LEAST SQUARES FAMILY OF ROUTINES.
C
C     DATA IS FROM DANIAL AND WOOD [1980], PAGES 428-441.
C
C     WRITTEN BY - JANET R. DONALDSON
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  AUGUST 3, 1987
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   LDSTAK
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  ARRAYS IN COMMON
      DOUBLE PRECISION DSTAK(12)
C
C  LOCAL SCALARS
      INTEGER
     +   IPRT,IXM,M,N,NPAR
C
C  LOCAL ARRAYS
      REAL
     +   PAR(5),RES(10),STP(5),XM(10,2),Y(10)
C
C  EXTERNAL SUBROUTINES
      EXTERNAL DCKLS,DRV1A,DRV1B,IPRINT,MDL1,NLS,NLSD,STPLS
C
C  COMMON BLOCKS
      COMMON /CSTAK/DSTAK
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     EXTERNAL DRV1A, DRV1B
C        THE NAME OF THE ''USER SUPPLIED'' DERIVATIVE ROUTINES.
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER IERR
C        THE INTEGER VALUE DESIGNATING WHETHER ANY ERRORS WERE
C        DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .EQ. 1, ERRORS WERE DETECTED.
C     INTEGER IPRT
C        LOGICAL OUTPUT UNIT.
C     INTEGER IXM
C        THE FIRST DIMENSION OF THE MATRIX X.
C     INTEGER LDSTAK
C        THE LENGTH OF DSTAK IN COMMON /CSTAK/.
C     INTEGER M
C        THE NUMBER OF INDEPENDENT VARIABLES.
C     EXTERNAL MDL1
C        THE NAME OF THE ''USER SUPPLIED'' MODEL ROUTINES.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS IN EACH PROBLEM.
C     INTEGER NPAR
C        THE NUMBER OF PARAMETERS TO BE ESTIMATED.
C     REAL PAR(5)
C        THE ARRAY IN WHICH THE CURRENT ESTIMATES OF THE UNKNOWN
C        PARAMETERS ARE STORED.
C     REAL RES(10)
C        THE RESIDUALS.
C     REAL STP(5)
C        THE STEP SIZES SELECTED FOR GENERATING FINITE DIFFERENCE
C        DERIVATIVES.
C     REAL XM(10,2)
C        THE INDEPENDENT VARIABLE.
C     REAL Y(10)
C        THE DEPENDENT VARIABLE.
C
C
      DATA Y(1), Y(2), Y(3), Y(4), Y(5), Y(6)
     +   /2.138E0, 3.421E0, 3.597E0, 4.340E0, 4.882E0, 5.660E0/
C
      DATA XM(1,1), XM(2,1), XM(3,1), XM(4,1), XM(5,1), XM(6,1)
     +   /1.309E0, 1.471E0, 1.490E0, 1.565E0, 1.611E0, 1.680E0/
C
C     SET PARAMETERS NECESSARY FOR THE COMPUTATIONS
C
      CALL IPRINT(IPRT)
      IXM = 10
      N = 6
      M = 1
      NPAR = 2
C
C     PRINT HEADER
C
      WRITE (IPRT,1000)
C
C     RUN SIMPLE EXAMPLE OF NLS
C
      WRITE (IPRT,1100)
      PAR(1) = 0.725
      PAR(2) = 4.000
      CALL NLS(Y, XM, N, M, IXM, MDL1, PAR, NPAR, RES, LDSTAK)
      WRITE (IPRT,2000) IERR
C
C     RUN SIMPLE EXAMPLE OF NLSD
C
      WRITE (IPRT,1200)
      PAR(1) = 0.725
      PAR(2) = 4.000
      CALL NLSD(Y, XM, N, M, IXM, MDL1, DRV1A, PAR, NPAR, RES, LDSTAK)
      WRITE (IPRT,2000) IERR
C
C     RUN SIMPLE EXAMPLE OF STPLS
C
      WRITE (IPRT,1300)
      PAR(1) = 0.725
      PAR(2) = 4.000
      CALL STPLS(XM, N, M, IXM, MDL1, PAR, NPAR, LDSTAK, STP)
      WRITE (IPRT,2000) IERR
C
C     RUN SIMPLE EXAMPLE OF DCKLS
C
      WRITE (IPRT,1400)
      PAR(1) = 0.000
      PAR(2) = 4.000
      CALL DCKLS(XM, N, M, IXM, MDL1, DRV1B, PAR, NPAR, LDSTAK)
      WRITE (IPRT,2000) IERR
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT ('1*CH9')
 1100 FORMAT (' SIMPLE TEST OF NLS')
 1200 FORMAT ('1SIMPLE TEST OF NLSD')
 1300 FORMAT ('1SIMPLE TEST OF STPLS')
 1400 FORMAT ('1SIMPLE TEST OF DCKLS')
 2000 FORMAT (/' THE VALUE OF IERR IS ', I4)
C
      END
