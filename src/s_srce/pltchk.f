*PLTCHK
      SUBROUTINE PLTCHK (YM, YMMISS, X, XMISS, N, M, IYM, MULTI,
     +   ILOG, YLB, YUB, XLB, XUB, NMSUB, MISS, XCHECK)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS PLOT FAMILY ROUTINE CHECKS FOR ERRORS IN THE PARAMETER LISTS
C     OF THE MULTIPLE PLOT ROUTINES
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
     +   XLB,XMISS,XUB,YLB,YUB
      INTEGER
     +   ILOG,IYM,M,N
      LOGICAL
     +   MISS,MULTI,XCHECK
C
C  ARRAY ARGUMENTS
      REAL
     +   X(*),YM(*),YMMISS(*)
      CHARACTER
     +   NMSUB(6)*1
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      INTEGER
     +   NV
      LOGICAL
     +   ERR01,ERR02,ERR03,ERR04,ERR05,ERR06,ERR07,ERR08,ERR09,HEAD
C
C  LOCAL ARRAYS
      INTEGER
     +   ILOGXY(2)
      CHARACTER
     +   LIYM(8)*1,LM(8)*1,LN(8)*1,LONE(8)*1,LX(8)*1,LXLB(8)*1,
     +   LXUB(8)*1,LY(8)*1,LYLB(8)*1,LYM(8)*1,LYUB(8)*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL EISGE,ERAGT,ERAGTM,ERSGT,ERVGT,ERVGTM,PRTCNT
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MAX
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     LOGICAL ERR01, ERR02, ERR03, ERR04, ERR05, ERR06, ERR07, ERR08,
C    1   ERR09
C        VALUES INDICATING WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT
C        (FALSE).
C     LOGICAL HEAD
C        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
C        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
C        OF HEAD WILL BE CHANGED TO FALSE.
C     INTEGER IERR
C        THE COMMON VARIABLE USED AS AN ERROR FLAG
C        IF = 0 THEN NO ERORRS
C     INTEGER ILOG
C        THE TWO DIGIT INTEGER, PQ, USED TO SELECT AXIS SCALE, WHERE
C        P DESIGNATES THE X-AXIS AND Q DESIGNATES THE Y-AXIS.
C        IF P.EQ.0 (Q.EQ.0), THEN THE X-AXIS (Y-AXIS) IS LINEAR.
C        IF P.NE.0 (Q.NE.0), THEN THE X-AXIS (Y-AXIS) IS LOG.
C     INTEGER ILOGXY(2)
C        ...
C     INTEGER IYM
C        ACTUAL ROW DIMENSION OF YM DECLARED IN USERS MAIN PROGRAM
C     CHARACTER*1 LIYM(8), LM(8), LN(8), LONE(8), LX(8), LXLB(8),
C    *  LXUB(8), LY(8), LYLB(8), LYM(8), LYUB(8)
C        THE ARRAY(S) CONTAINING THE NAME(S) OF THE INPUT PARAMETERS(S)
C        CHECKED FOR ERRORS.
C     INTEGER M
C        THE NUMBER OF VECTORS IN YM
C     LOGICAL MISS
C        INDICATOR VARIABLE USED TO DESIGNATE WHETHER MISSING VALUES
C        MAY BE PRESENT (MISS = .TRUE.) OR NOT (MISS = .FALSE.)
C     LOGICAL MULTI
C        AN INDICATOR VARIABLE USED TO DESIGNATE WHETHER THE CALLING
C        ROUTINE HAS AN M PREFIX (TRUE) OR NOT (FALSE).
C     INTEGER N
C        THE LENGTH OF THE VECTORS
C     CHARACTER*1 NMSUB(6)
C        THE CHARACTERS OF THE CALLING ROUTINES NAME.
C     INTEGER NV
C        THE NUMBER OF VIOLATIONS FOUND IN THE X AND Y AXIS ARRAYS.
C     REAL X(N)
C        VECTOR OF OBSERVATIONS FOR X COORDINATES
C     LOGICAL XCHECK
C        INDICATOR VARIABLE USED TO DESIGNATE WHETHER X-AXIS VALUES
C        ARE TO BE CHECKED (XCHECK = .TRUE.) OR NOT (XCHECK = .FALSE.)
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
C     REAL YM(IYM,M)
C        MULTIVARIATE OBSERVATIONS FOR THE Y (VERTICAL) COORDINATES.
C     REAL YMMISS(M)
C        THE MISSING VALUE CODE FOR EACH COLUMN OF YM.
C     REAL YUB
C        THE UPPER BOUND FOR THE Y-AXIS.  (YLB=YUB INDICATES LIMITS ARE
C        TO BE DETERMINED FROM THE RANGE OF THE DATA.)
C
C
C     SET UP NAME ARRAYS
C
      DATA
     +   LIYM(1),   LIYM(2),   LIYM(3),   LIYM(4),   LIYM(5),   LIYM(6)
     + /     'I',       'Y',       'M',       ' ',       ' ',       ' '/
      DATA LIYM(7), LIYM(8)
     +   /     ' ',     ' '/
      DATA
     +     LM(1),     LM(2),     LM(3),     LM(4),     LM(5),     LM(6)
     + /     'M',       ' ',       ' ',       ' ',       ' ',       ' '/
      DATA   LM(7),   LM(8)
     +   /     ' ',     ' '/
      DATA
     +     LN(1),     LN(2),     LN(3),     LN(4),     LN(5),     LN(6)
     + /     'N',       ' ',       ' ',       ' ',       ' ',       ' '/
      DATA   LN(7),   LN(8)
     +   /     ' ',     ' '/
      DATA   LONE(1),  LONE(2),  LONE(3),  LONE(4),  LONE(5),  LONE(6),
     +       LONE(7),  LONE(8)/'O', 'N', 'E', ' ', ' ', ' ', ' ', ' '/
      DATA
     +     LX(1),     LX(2),     LX(3),     LX(4),     LX(5),     LX(6)
     + /     'X',       ' ',       ' ',       ' ',       ' ',       ' '/
      DATA   LX(7),   LX(8)
     +   /     ' ',     ' '/
      DATA
     +   LXLB(1),   LXLB(2),   LXLB(3),   LXLB(4),   LXLB(5),   LXLB(6)
     + /     'X',       'L',       'B',       ' ',       ' ',       ' '/
      DATA LXLB(7), LXLB(8)
     +   /     ' ',     ' '/
      DATA
     +   LXUB(1),   LXUB(2),   LXUB(3),   LXUB(4),   LXUB(5),   LXUB(6)
     + /     'X',       'U',       'B',       ' ',       ' ',       ' '/
      DATA LXUB(7), LXUB(8)
     +   /     ' ',     ' '/
      DATA
     +     LY(1),     LY(2),     LY(3),     LY(4),     LY(5),     LY(6)
     + /     'Y',       ' ',       ' ',       ' ',       ' ',       ' '/
      DATA   LY(7),   LY(8)
     +   /     ' ',     ' '/
      DATA
     +   LYLB(1),   LYLB(2),   LYLB(3),   LYLB(4),   LYLB(5),   LYLB(6)
     + /     'Y',       'L',       'B',       ' ',       ' ',       ' '/
      DATA LYLB(7), LYLB(8)
     +   /     ' ',     ' '/
      DATA
     +    LYM(1),    LYM(2),    LYM(3),    LYM(4),    LYM(5),    LYM(6)
     + /     'Y',       'M',       ' ',       ' ',       ' ',       ' '/
      DATA  LYM(7),  LYM(8)
     +   /     ' ',     ' '/
      DATA
     +   LYUB(1),   LYUB(2),   LYUB(3),   LYUB(4),   LYUB(5),   LYUB(6)
     + /     'Y',       'U',       'B',       ' ',       ' ',       ' '/
      DATA LYUB(7), LYUB(8)
     +   /     ' ',     ' '/
C
C     COMMENCE BODY OF ROUTINE
C
      IERR = 0
      HEAD = .TRUE.
C
C     NUMBER OF POINTS MUST BE AT LEAST 1
C
      CALL EISGE(NMSUB, LN, N, 1, 2, HEAD, ERR01, LONE)
C
C     THERE MUST BE AT LEAST 1 COLUMN OF VECTORS
C
      CALL EISGE(NMSUB, LM, M, 1, 2, HEAD, ERR02, LONE)
C
C     THE ACTUAL LENGTH OF YM MUST EQUAL OR EXCEED THE NUMBER OF
C     OBSERVATIONS
C
      ERR03 = .TRUE.
      IF (.NOT.ERR01)
     +   CALL EISGE(NMSUB, LIYM, IYM, N, 3, HEAD, ERR03, LN)
C
C     IF THIS IS A LOG PLOT CHECK FOR NON-POSITIVE VALUES IN DATA
C
      IF (ERR01 .OR. ERR02 .OR. ERR03) IERR = 1
      IF (ILOG .LE. 0) RETURN
C
      ERR04 = .FALSE.
      ERR05 = .FALSE.
      ERR06 = .FALSE.
      ERR07 = .FALSE.
      ERR08 = .FALSE.
      ERR09 = .FALSE.
C
      CALL PRTCNT (MAX(0,ILOG),2,ILOGXY)
      IF ((ILOGXY(1).NE.0) .AND. XCHECK) THEN
        IF (.NOT.ERR01) THEN
C
C         IF X AXIS IS LOG SCALE, CHECK FOR NEGATIVE X AXIS VALUES
C
          IF (MISS) THEN
            CALL ERVGTM(NMSUB, LX, X, XMISS, N, 0.0E0, 0, HEAD, 1,
     +        NV, ERR04, LX)
          ELSE
            CALL ERVGT(NMSUB, LX, X, N, 0.0E0, 0, HEAD, 1, NV, ERR04,
     +        LX)
          END IF
        END IF
C
        IF (XLB.LT.XUB) THEN
C
C         CHECK FOR NEGATIVE PLOT BOUNDS
C
          CALL ERSGT(NMSUB, LXLB, XLB, 0.0E0, 1, HEAD, ERR05, LXLB)
          CALL ERSGT(NMSUB, LXUB, XUB, 0.0E0, 1, HEAD, ERR06, LXUB)
        END IF
      END IF
      IF (ILOGXY(2).NE.0) THEN
        IF ((.NOT.ERR01) .AND. (.NOT.ERR02) .AND. (.NOT.ERR03)) THEN
C
C         IF Y AYIS IS LOG SCALE, CHECK FOR NEGATIVE Y AYIS VALUES
C
          IF (MISS) THEN
            IF (MULTI) THEN
              CALL ERAGTM(NMSUB, LYM, YM, YMMISS, N, M, IYM, 0.0E0, 0,
     +           HEAD, 1, NV, ERR04, LYM)
            ELSE
              CALL ERVGTM(NMSUB, LY, YM, YMMISS(1), N, 0.0E0, 0, HEAD,
     +           1, NV, ERR04, LY)
            END IF
          ELSE
            IF (MULTI) THEN
              CALL ERAGT(NMSUB, LYM, YM, N, M, IYM, 0.0E0, 0, HEAD,
     +           1, NV, ERR04, LYM)
            ELSE
              CALL ERVGT(NMSUB, LY, YM, N, 0.0E0, 0, HEAD, 1,
     +           NV, ERR04, LY)
            END IF
          END IF
        END IF
C
        IF (YLB.LT.YUB) THEN
C
C         CHECK FOR NEGATIVE PLOT BOUNDS
C
          CALL ERSGT(NMSUB, LYLB, YLB, 0.0E0, 1, HEAD, ERR05, LYLB)
          CALL ERSGT(NMSUB, LYUB, YUB, 0.0E0, 1, HEAD, ERR06, LYUB)
        END IF
      END IF
C
      IF (ERR04 .OR. ERR05 .OR. ERR06 .OR. ERR07 .OR. ERR08 .OR. ERR09)
     +   IERR = 1
C
      RETURN
C
      END
