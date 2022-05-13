*SVPML
      SUBROUTINE SVPML(YM, YMMISS, N, NS, ISYM, ILOG)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE USER CALLABLE ROUTINE WHICH PRODUCES A VERTICAL
C     PLOT WITH MISSING DATA AND USER CONTROL OF THE PLOT SYMBOL USED
C     (LOG PLOT OPTION).
C
C     WRITTEN BY - JANET R. DONALDSON
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  JANUARY 21, 1982
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   ILOG,N,NS
C
C  ARRAY ARGUMENTS
      REAL
     +   YM(*),YMMISS(*)
      INTEGER
     +   ISYM(*)
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      REAL
     +   XINC,XLB,YLB,YUB
      INTEGER
     +   IBAR,IPRT,IRLIN,ISCHCK,ISIZE,IYM,LISYM,M
      LOGICAL
     +   MISS,MULTI
C
C  LOCAL ARRAYS
      CHARACTER
     +   NMSUB(6)*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL IPRINT,VPCNT
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     INTEGER IBAR
C        THE VARIABLE USED TO DETERMINE IF SINGLE POINTS (IBAR .EQ. 0)
C        OR BARS (IBAR .NE. 0) ARE TO BE PLOTTED.
C     INTEGER IERR
C        A COMMON VARIABLE USED AS A FLAG TO INDICATE WHETHER
C        OR NOT THERE ARE ANY ERRORS, IF =0 THEN NO ERRORS.
C     INTEGER ILOG
C        THE TWO DIGIT INTEGER, PQ, USED TO SELECT AXIS SCALE, WHERE
C        P DESIGNATES THE X-AXIS AND Q DESIGNATES THE Y-AXIS.
C        IF P.EQ.0 (Q.EQ.0), THEN THE X-AXIS (Y-AXIS) IS LINEAR.
C        IF P.NE.0 (Q.NE.0), THEN THE X-AXIS (Y-AXIS) IS LOG.
C     INTEGER IPRT
C        OUTPUT LOGICAL UNIT NUMBER
C     INTEGER IRLIN
C        THE INDICATOR VARIABLE USED TO DESIGNATE WHETHER ZERO OR THE
C        SERIES MEAN IS TO BE PLOTTED AS A REFERENCE LINE, OR WHETHER
C        NO REFERENCE LINE IS TO BE PLOTTED.
C        IF IRLIN .LE. -1, NO REFERENCE LINE IS PLOTTED.
C        IF IRLIN .EQ.  0, ZERO IS PLOTTED AS THE REFERENCE LINE.
C        IF IRLIN .GE.  1, THE SERIES MEAN IS PLOTTED.
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
C        VECTOR CONTAINING SYMBOL DESIGNATIONS FOR PLOTTING
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
C     INTEGER NS
C        THE SAMPLING FREQUENCY,
C        WHERE IF NS .LE. 1, EVERY POINT IS PLOTTED,
C                       = 2, EVERY OTHER POINT IS PLOTTED,
C                       = 3, EVERY THIRD POINT IS PLOTTED, ETC.
C     REAL XINC, XLB
C        INCREMENT AND LOWER BOUNDS FOR X-AXIS.
C     REAL YLB
C        LOWER BOUND FOR Y-AXIS.
C     REAL YM(N,1)
C        MULTIVARIATE OBSERVATIONS FOR THE Y COORDINATES
C     REAL YMMISS(1)
C        THE MISSING VALUE CODE FOR THE Y-AXIS.
C     REAL YUB
C        UPPER BOUND FOR Y-AXIS.
C
C
C     SET UP NAME ARRAYS
C
      DATA
     +  NMSUB(1),  NMSUB(2),  NMSUB(3),  NMSUB(4),  NMSUB(5),  NMSUB(6)
     +    /  'S',       'V',       'P',       'M',       'L',       ' '/
C
C     SET DEFAULT VALUES
C
      M = 1
      IYM = N
      MULTI = .FALSE.
      YLB = 0.0E0
      YUB = 0.0E0
      XLB = 1.0E0
      XINC = 1.0E0
      ISCHCK = 1
      ISIZE = -1
      MISS = .TRUE.
      LISYM = N
      IRLIN = -1
      IBAR = 0
C
C     COMMENCE BODY OF ROUTINE
C
      CALL VPCNT (YM, YMMISS, N, M, IYM, MULTI, ILOG, YLB, YUB,
     +  XLB, XINC, NS, IRLIN, IBAR, NMSUB, ISCHCK, ISYM, ISIZE,
     +  MISS, LISYM)
C
      IF (IERR.NE.0) THEN
        IERR = 1
        CALL IPRINT(IPRT)
        WRITE (IPRT,1000)
      END IF
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (/42H THE CORRECT FORM OF THE CALL STATEMENT IS//
     +  '       CALL SVPML (Y, YMISS, N, NS, ISYM, ILOG)')
      END
