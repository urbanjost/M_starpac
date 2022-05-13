*PPM
      SUBROUTINE PPM(YM, YMMISS, X, XMISS, N)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE USER CALLABLE ROUTINE WHICH PRODUCES A SIMPLE PAGE
C     PLOT FOR DATA WITH MISSING OBSERVATIONS (SHORT CALL).
C
C     WRITTEN BY - LINDA L. MITCHELL AND JANET R. DONALDSON
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  JANUARY 21, 1982
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      REAL
     +   XMISS
      INTEGER
     +   N
C
C  ARRAY ARGUMENTS
      REAL
     +   X(*),YM(*),YMMISS(1)
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      REAL
     +   XLB,XUB,YLB,YUB
      INTEGER
     +   ILOG,IPRT,ISCHCK,ISIZE,IYM,LISYM,M,NOUT
      LOGICAL
     +   MISS,MULTI
C
C  LOCAL ARRAYS
      INTEGER
     +   ISYM(1)
      CHARACTER
     +   NMSUB(6)*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL IPRINT,PPCNT
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     INTEGER IERR
C        THE COMMON VARIABLE USED TO INDICATE ERRORS,
C        IF =0, THEN NO ERRORS
C     INTEGER ILOG
C        THE TWO DIGIT INTEGER, PQ, USED TO SELECT AXIS SCALE, WHERE
C        P DESIGNATES THE X-AXIS AND Q DESIGNATES THE Y-AXIS.
C        IF P.EQ.0 (Q.EQ.0), THEN THE X-AXIS (Y-AXIS) IS LINEAR.
C        IF P.NE.0 (Q.NE.0), THEN THE X-AXIS (Y-AXIS) IS LOG.
C     INTEGER IPRT
C        OUTPUT LOGICAL UNIT NUMBER
C     INTEGER ISCHCK
C        THE INTEGER VALUE INDICATING HOW THE PLOTTING SYMBOLS
C           WILL BE DESIGNATED, WHERE
C           0 INDICATES THE PLOTTING SYMBOLS HAVE NOT BEEN DESIGNATED IN
C             THE N VECTOR ISYM AND ONLY THE SYMBOL + IS TO BE USED
C           1 INDICATES THE PLOTTING SYMBOLS HAVE BEEN DESIGNATED IN THE
C             N VECTOR ISYM
C           2 INDICATES THAT M SERIES ARE BEING PLOTTED.
C             SYMBOL I+4 WILL BE USED FOR COLUMN I OF YM.
C     INTEGER ISIZE
C        THE TWO DIGIT INTEGER, PQ, USED TO SELECT AXIS SIZE, WHERE
C        P DESIGNATES THE X-AXIS AND Q DESIGNATES THE Y-AXIS.
C        IF P.EQ.0 (Q.EQ.0), THEN THE X-AXIS (Y-AXIS) IS THE MAXIMUM.
C        IF P.NE.0 (Q.NE.0), THEN THE X-AXIS (Y-AXIS) IS HALF THE MAXIMU
C     INTEGER ISYM(1)
C        VECTOR CONTAINING SYMBOLS FOR PLOTTING.
C     INTEGER IYM
C        THE FIRST DIMENSION OF ARRAY YM.
C     INTEGER LISYM
C        THE LENGTH OF ARRAY ISYM.
C     INTEGER M
C        NUMBER OF Y VECTORS
C     LOGICAL MISS
C        INDICATOR VARIABLE USED TO DESIGNATE WHETHER MISSING VALUES
C        MAY BE PRESENT (MISS = .TRUE.) OR NOT (MISS = .FALSE.)
C     LOGICAL MULTI
C        INDICATOR VARIABLE USED TO DESIGNATE WHETHER MULTIPLE Y VALUES
C        ARE TO BE PLOTTED (MULTI = .TRUE.) OR NOT (MULTI = .FALSE.)
C     INTEGER N
C        LENGTH OF VECTORS
C     CHARACTER*1 NMSUB(6)
C        THE CHARACTERS OF THE CALLING ROUTINES NAME.
C     INTEGER NOUT
C        USED TO INDICATE HOW MANY OF THE POINTS OUTSIDE THE BOUNDS
C        OF THE PLOT ARE TO BE LISTED.
C     REAL X(N)
C        VECTOR OF OBSERVATIONS FOR X COORDINATES
C     REAL XLB
C        THE LOWER BOUND FOR THE X-AXIS.  (XLB=XUB INDICATES LIMITS ARE
C        TO BE DETERMINED FROM THE RANGE OF THE DATA.)
C     REAL XMISS
C        THE MISSING VALUE CODE FOR THE X-AXIS.
C     REAL XUB
C        THE UPPER BOUND FOR THE X-AXIS.  (XLB=XUB INDICATES LIMITS ARE
C        TO BE DETERMINED FROM THE RANGE OF THE DATA.)
C     REAL YLB
C        THE LOWER BOUND FOR THE Y-AXIS.  (YLB=YUB INDICATES LIMITS ARE
C        TO BE DETERMINED FROM THE RANGE OF THE DATA.)
C     REAL YM(N,1)
C        VECTOR OF OBSERVATIONS FOR THE Y COORDINATES
C     REAL YMMISS(1)
C        THE MISSING VALUE CODE FOR THE Y-AXIS.
C     REAL YUB
C        THE UPPER BOUND FOR THE Y-AXIS.  (YLB=YUB INDICATES LIMITS ARE
C        TO BE DETERMINED FROM THE RANGE OF THE DATA.)
C
C
C     SET UP NAME ARRAYS
C
      DATA
     +  NMSUB(1),  NMSUB(2),  NMSUB(3),  NMSUB(4),  NMSUB(5),  NMSUB(6)
     +    /  'P',       'P',       'M',       ' ',       ' ',       ' '/
C
C     COMMENCE BODY OF ROUTINE
C
C     SET DEFAULT VALUES
C
      M = 1
      IYM = N
      MULTI = .FALSE.
      ILOG = -1
      YLB = 0.0E0
      YUB = 0.0E0
      XLB = 0.0E0
      XUB = 0.0E0
      ISCHCK = 0
      ISIZE = -1
      NOUT = 0
      MISS = .TRUE.
      LISYM = 1
C
      CALL PPCNT (YM, YMMISS, X, XMISS, N, M, IYM, MULTI, ILOG,
     +  YLB, YUB, XLB, XUB, NMSUB, ISCHCK, ISYM, ISIZE, NOUT, MISS,
     +  LISYM)
C
      IF (IERR.NE.0) THEN
        IERR = 1
        CALL IPRINT(IPRT)
        WRITE (IPRT,1000)
      END IF
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (/42H THE CORRECT FORM OF THE CALL STATEMENT IS//
     +   '       CALL PPM (Y, YMISS, X, XMISS, N)')
      END
