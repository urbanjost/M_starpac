*SPPC
      SUBROUTINE SPPC(YM, X, N, ISYM, ILOG, ISIZE, NOUT, YLB, YUB,
     +  XLB, XUB)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE USER CALLABLE ROUTINE WHICH PRODUCES A SIMPLE PAGE
C     PLOT WITH USER CONTROL OF PLOT SYMBOLS (LONG CALL).
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
      DOUBLE PRECISION
     +   XLB,XUB,YLB,YUB
      INTEGER
     +   ILOG,ISIZE,N,NOUT
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   X(*),YM(*)
      INTEGER
     +   ISYM(*)
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   XMISS
      INTEGER
     +   IPRT,ISCHCK,IYM,LISYM,M
      LOGICAL
     +   MISS,MULTI
C
C  LOCAL ARRAYS
      DOUBLE PRECISION
     +   YMMISS(1)
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
C     INTEGER ISYM(N)
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
C     DOUBLE PRECISION X(N)
C        VECTOR OF OBSERVATIONS FOR X COORDINATES
C     DOUBLE PRECISION XLB
C        THE LOWER BOUND FOR THE X-AXIS.  (XLB=XUB INDICATES LIMITS ARE
C        TO BE DETERMINED FROM THE RANGE OF THE DATA.)
C     DOUBLE PRECISION XMISS
C        THE MISSING VALUE CODE FOR THE X-AXIS.
C     DOUBLE PRECISION XUB
C        THE UPPER BOUND FOR THE X-AXIS.  (XLB=XUB INDICATES LIMITS ARE
C        TO BE DETERMINED FROM THE RANGE OF THE DATA.)
C     DOUBLE PRECISION YLB
C        THE LOWER BOUND FOR THE Y-AXIS.  (YLB=YUB INDICATES LIMITS ARE
C        TO BE DETERMINED FROM THE RANGE OF THE DATA.)
C     DOUBLE PRECISION YM(N,1)
C        VECTOR OF OBSERVATIONS FOR THE Y COORDINATES
C     DOUBLE PRECISION YMMISS(1)
C        THE MISSING VALUE CODE FOR THE Y-AXIS.
C     DOUBLE PRECISION YUB
C        THE UPPER BOUND FOR THE Y-AXIS.  (YLB=YUB INDICATES LIMITS ARE
C        TO BE DETERMINED FROM THE RANGE OF THE DATA.)
C
C
C     SET UP NAME ARRAYS
C
      DATA
     +  NMSUB(1),  NMSUB(2),  NMSUB(3),  NMSUB(4),  NMSUB(5),  NMSUB(6)
     +    /  'S',       'P',       'P',       'C',       ' ',       ' '/
C
C     COMMENCE BODY OF ROUTINE
C
C     SET DEFAULT VALUES
C
      YMMISS(1) = 1.0D0
      XMISS = 1.0D0
      M = 1
      IYM = N
      MULTI = .FALSE.
      ISCHCK = 1
      MISS = .FALSE.
      LISYM = N
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
     +  '       CALL SPPC (Y, X, N, ISYM, ILOG,'/
     +  '      +           ISIZE, NOUT, YLB, YUB, XLB, XUB)')
      END
