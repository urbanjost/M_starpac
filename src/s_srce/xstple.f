*XSTPLE
      SUBROUTINE XSTPLE(LDSTAK)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     ROUTINE TO TEST DERIVATIVE CHECKING ROUTINES.
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  NOVEMBER 21, 1983
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
      REAL
     +   DELTA,EXMPT
      INTEGER
     +   I,IPRT,IXM,LDSMIN,M,N,NETA,NPAR,NPRT,NTEST
C
C  LOCAL ARRAYS
      REAL
     +   PAR(10),SCALE(10),STP(10),XM(200,2)
C
C  EXTERNAL SUBROUTINES
      EXTERNAL IPRINT,LDSCMP,MDL4,STPLS,STPLS1,STPLS2,STPLSC
C
C  COMMON BLOCKS
      COMMON /CSTAK/DSTAK
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     REAL DELTA
C        *
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     REAL EXMPT
C        THE PROPORTION OF OBSERVATIONS FOR WHICH THE COMPUTED
C        NUMERICAL DERIVATIVES WRT A GIVEN PARAMETER ARE EXEMPTED
C        FROM MEETING THE DERIVATIVE ACCEPTANCE CRITERIA.
C     INTEGER I
C        AN INDEX VARIABLE.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .EQ. 1, ERRORS WERE DETECTED.
C     INTEGER IPRT
C        THE UNIT NUMBER FOR PRINTED OUTPUT.
C     INTEGER IXM
C        THE FIRST DIMENSION OF THE INDEPENDENT VARIABLE ARRAY XM.
C     INTEGER LDSMIN
C        THE MINIMUM LENGTH OF THE ARRAY DSTAK ALLOWED.
C     INTEGER LDSTAK
C        THE LENGTH OF THE ARRAY DSTAK.
C     INTEGER M
C        THE NUMBER OF INDEPENDENT VARIABLES.
C     EXTERNAL MDL4
C        THE NAME OF THE USER SUPPLIED SUBROUTINE WHICH COMPUTES THE
C        PREDICTED VALUES BASED ON THE CURRENT PARAMETER ESTIMATES.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS OF DATA.
C     INTEGER NETA
C        THE NUMBER OF RELIABLE DIGITS IN THE MODEL.
C     INTEGER NPAR
C        THE NUMBER OF UNKNOWN PARAMETERS IN THE MODEL.
C     INTEGER NPRT
C        THE INDICATOR VARIABLE USED TO SPECIFY WHETHER OR NOT
C        PRINTED OUTPUT IS TO BE PROVIDED, WHERE IF THE VALUE OF
C        NPRT IS ZERO, NO PRINTED OUTPUT IS GIVEN.
C     INTEGER NTEST
C        THE NUMBER OF THE CURRENT TEST.
C     REAL PAR(10)
C        THE ARRAY IN WHICH THE CURRENT ESTIMATES OF THE UNKNOWN
C        PARAMETERS ARE STORED.
C     REAL SCALE(10)
C        A DUMMY ARRAY, INDICATING USE OF DEFAULT VALUES FOR
C        THE TYPICAL SIZE OF THE PARAMETERS.
C     REAL STP(10)
C        THE SELECTED STEP SIZES FOR EACH PARAMETER.
C     REAL XM(200,2)
C        THE ARRAY IN WHICH ONE ROW OF THE INDEPENDENT VARIABLE ARRAY
C        IS STORED.
C
      CALL IPRINT(IPRT)
C
C     SET PARAMETER VALUES
C
      CALL STPLS1(N, M, IXM, PAR, NPAR, NETA, EXMPT, SCALE, NPRT)
      CALL STPLS2(NPAR, STP)
      CALL LDSCMP(14, 0, 2*(N+NPAR), 0, 0, 0, 'S', 10*N, LDSMIN)
C
      IF (LDSMIN.LE.LDSTAK) GO TO 5
C
      WRITE (IPRT, 1040) LDSMIN
      RETURN
C
    5 CONTINUE
C
C     CREATE INDEPENDENT VARIABLE
C
      DELTA = 0.0625E0
      XM(1,1) = 0.0E0
      DO 10 I=2,N
         XM(I,1) = XM(I-1,1) + DELTA
   10 CONTINUE
C
      NTEST = 0
C     CHECK ERROR HANDLING
C
C        TEST 1  -  MISCELANEOUS ERROR CHECKING
C
      N = -5
      M = -5
      IXM = -10
      NPAR = -10
C
C
      NTEST = NTEST + 1
      WRITE (IPRT,1090) NTEST
      WRITE (IPRT,1020)
      WRITE (IPRT,1000)
      IERR = -1
      CALL STPLS(XM, N, M, IXM, MDL4, PAR, NPAR, LDSTAK, STP)
      WRITE (IPRT,1050) IERR
      WRITE (IPRT,1010)
      IERR = -1
      CALL STPLSC(XM, N, M, IXM, MDL4, PAR, NPAR, LDSTAK, STP, NETA,
     +   EXMPT, SCALE, NPRT)
      WRITE (IPRT,1050) IERR
C
C        TEST 2  -  MISCELANEOUS ERROR CHECKING (CONTINUED)
C
      CALL STPLS1(N, M, IXM, PAR, NPAR, NETA, EXMPT, SCALE, NPRT)
      SCALE(2) = 0.0E0
C
      NTEST = NTEST + 1
      WRITE (IPRT,1090) NTEST
      WRITE (IPRT,1030)
      WRITE (IPRT,1000)
      IERR = -1
      CALL STPLS(XM, N, M, IXM, MDL4, PAR, NPAR, LDSMIN-1, STP)
      WRITE (IPRT,1050) IERR
      WRITE (IPRT,1010)
      IERR = -1
      CALL STPLSC(XM, N, M, IXM, MDL4, PAR, NPAR, LDSMIN-1, STP, NETA,
     +   EXMPT, SCALE, NPRT)
      WRITE (IPRT,1050) IERR
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (15H TEST OF STPLS )
 1010 FORMAT (15H TEST OF STPLSC)
 1020 FORMAT (32H CHECK ERROR HANDLING  -  TEST 1)
 1030 FORMAT (32H CHECK ERROR HANDLING  -  TEST 2)
 1040 FORMAT (45H1 *** LDSTAK MUST BE GREATER THAN OR EQUAL TO , I6)
 1050 FORMAT (/29H ***** RETURNED RESULTS *****, 5X, 15H (-1 INDICATES ,
     +   39HVALUE NOT CHANGED BY CALLED SUBROUTINE)//9H IERR IS , I3)
 1090 FORMAT (54H1DERIVATIVE STEP SIZE SELECTION SUBROUTINE TEST NUMBER,
     +   I5)
      END
