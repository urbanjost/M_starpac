*STPLS
      SUBROUTINE STPLS(XM, N, M, IXM, MDL, PAR, NPAR, LDSTAK, STP)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE USER CALLABLE SUBROUTINE FOR SELECTING STEP SIZES
C     TO BE USED IN COMPUTING FORWARD DIFFERENCE QUOTIENT ESTIMATES
C     OF THE NUMERICAL DERIVATIVES FOR THE NONLINEAR LEAST SQUARES
C     ROUTINES (SHORT CALL).
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
     +   PAR(*),STP(*),XM(*)
C
C  SUBROUTINE ARGUMENTS
      EXTERNAL MDL
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  ARRAYS IN COMMON
      DOUBLE PRECISION DSTAK(12)
C
C  LOCAL SCALARS
      REAL
     +   EXMPT
      INTEGER
     +   IPRT,LSCALE,NETA,NPRT
C
C  LOCAL ARRAYS
      REAL
     +   SCALE(1)
      CHARACTER
     +   NMSUB(6)*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL IPRINT,STPDRV
C
C  COMMON BLOCKS
      COMMON /CSTAK/DSTAK
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     REAL EXMPT
C        THE PROPORTION OF OBSERVATIONS FOR WHICH THE COMPUTED
C        NUMERICAL DERIVATIVES WRT A GIVEN PARAMETER ARE EXCEPTED
C        FROM MEETING THE DERIVATIVE ACCEPTANCE CRITERIA.
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
C        THE NUMBER OF OBSERVATIONS.
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
C     REAL PAR(NPAR)
C        THE ARRAY IN WHICH THE CURRENT ESTIMATES OF THE
C        PARAMETERS ARE STORED.
C     REAL SCALE(1)
C        A DUMMY VECTOR USED TO DESIGNATE USE OF THE DEFAULT VALUES OF
C        THE TYPICAL SIZE OF THE PARAMETERS.
C     REAL STP(NPAR)
C        THE SELECTED STEP SIZES.
C     REAL XM(IXM,M)
C        THE INDEPENDENT VARIABLE ARRAY
C
C     SET UP NAME ARRAYS
C
      DATA NMSUB(1), NMSUB(2), NMSUB(3), NMSUB(4), NMSUB(5), NMSUB(6) /
     +   'S','T','P','L','S',' '/
C
C     SET UP DEFAULT VALUES
C
      EXMPT = 0.1E0
      NETA = 0
      SCALE(1) = 0.0E0
      LSCALE = 1
      NPRT = 1
C
C     PASS CONTROL TO STEP SIZE SELECTION DRIVER
C
      CALL STPDRV(NMSUB, XM, N, M, IXM, MDL, PAR, NPAR, LDSTAK, STP,
     +   NETA, EXMPT, SCALE, LSCALE, NPRT)
C
      IF (IERR.NE.1) RETURN
C
      CALL IPRINT(IPRT)
      WRITE (IPRT,1000)
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (/42H THE CORRECT FORM OF THE CALL STATEMENT IS//
     +   '       CALL STPLS (XM, N, M, IXM, NLSMDL, PAR, NPAR, LDSTAK,',
     +   ' STP)')
      END
