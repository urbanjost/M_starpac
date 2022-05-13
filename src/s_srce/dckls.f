*DCKLS
      SUBROUTINE DCKLS(XM, N, M, IXM, MDL, DRV, PAR, NPAR, LDSTAK)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE USER CALLABLE SUBROUTINE FOR CHECKING USER SUPPLIED
C     ANALYTIC DERIVATIVES AGAINST NUMERICAL DERIVATIVES
C     FOR THE NONLINEAR LEAST SQUARES ROUTINES (SHORT CALL).
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  OCTOBER 3, 1983
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   IXM,LDSTAK,M,N,NPAR
C
C  ARRAY ARGUMENTS
      REAL
     +   PAR(*),XM(*)
C
C  SUBROUTINE ARGUMENTS
      EXTERNAL DRV,MDL
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
     +   IPRT,LSCALE,NETA,NPRT,NROW,NTAU
C
C  LOCAL ARRAYS
      REAL
     +   SCALE(1)
      CHARACTER
     +   NMSUB(6)*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL DCKDRV,IPRINT
C
C  COMMON BLOCKS
      COMMON /CSTAK/DSTAK
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     EXTERNAL DRV
C        THE NAME OF THE USER SUPPLIED SUBROUTINE WHICH COMPUTES THE
C        ANALYTIC DERIVATIVES (JACOBIAN MATRIX) OF THE MODEL.
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .EQ. 1, ERRORS WERE DETECTED.
C     INTEGER IPRT
C        THE UNIT NUMBER FOR PRINTED OUTPUT.
C     INTEGER IXM
C        THE FIRST DIMENSION OF THE INDEPENDENT VARIABLE ARRAY XM.
C     INTEGER LDSTAK
C        THE LENGTH OF THE ARRAY DSTAK.
C     INTEGER LSCALE
C        THE LENGTH OF VECTOR SCALE.
C     INTEGER M
C        THE NUMBER OF INDEPENDENT VARIABLES.
C     EXTERNAL MDL
C        THE NAME OF THE USER SUPPLIED SUBROUTINE WHICH COMPUTES THE
C        PREDICTED VALUES BASED ON THE CURRENT PARAMETER ESTIMATES.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS OF DATA.
C     INTEGER NETA
C        THE NUMBER OF ACCURATE DIGITS IN THE MODEL.
C     CHARACTER*1 NMSUB(6)
C        THE NAME OF THE SUBROUTINE CALLING THE ERROR CHECKING
C        SUBROUTINES.
C     INTEGER NPAR
C        THE NUMBER OF PARAMETERS IN THE MODEL.
C     INTEGER NPRT
C        THE INDICATOR VARIABLE USED TO SPECIFY WHETHER OR NOT
C        PRINTED OUTPUT IS TO BE PROVIDED, WHERE IF THE VALUE OF
C        NPRT IS ZERO, NO PRINTED OUTPUT IS GIVEN.
C     INTEGER NROW
C        THE NUMBER OF THE ROW OF THE INDEPENDENT VARIABLE ARRAY AT
C        WHICH THE DERIVATIVE IS TO BE CHECKED.
C     INTEGER NTAU
C        THE NUMBER OF DIGITS OF AGREEMENT REQUIRED BETWEEN THE
C        NUMERICAL DERIVATIVES AND THE USER SUPPLIED DERIVATIVES.
C     REAL PAR(NPAR)
C        THE ARRAY IN WHICH THE CURRENT ESTIMATES OF THE
C        PARAMETERS ARE STORED.
C     REAL SCALE(1)
C        A DUMMY ARRAY, INDICATING USE OF DEFAULT VALUES FOR
C        THE TYPICAL SIZE OF THE PARAMETERS.
C     REAL XM(IXM,M)
C        THE INDEPENDENT VARIABLE MATRIX.
C
C     SET UP NAME ARRAYS
C
      DATA NMSUB(1), NMSUB(2), NMSUB(3), NMSUB(4), NMSUB(5), NMSUB(6) /
     +   'D','C','K','L','S',' '/
C
C     SET DEFAULT VALUES
C
      NETA = 0
      NTAU = 0
      SCALE(1) = 0.0E0
      LSCALE = 1
      NPRT = 1
      NROW = 0
C
C     PASS CONTROL OF DERIVATIVE CHECKING TO DCKDRV
C
      CALL DCKDRV(NMSUB, LDSTAK, XM, N, M, IXM, MDL, DRV, PAR, NPAR,
     +   NETA, NTAU, SCALE, LSCALE, NROW, NPRT)
C
      IF (IERR.NE.1) RETURN
C
      CALL IPRINT(IPRT)
      WRITE (IPRT,1000)
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (/42H THE CORRECT FORM OF THE CALL STATEMENT IS//
     +   '       CALL DCKLS (XM, N, M, IXM, NLSMDL, NLSDRV,'/
     +   '      +            PAR, NPAR, LDSTAK)')
      END
