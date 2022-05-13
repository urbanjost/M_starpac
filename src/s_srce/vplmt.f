*VPLMT
      SUBROUTINE VPLMT (YM, YMMISS, N, M, IYM, YLB, YUB, YMN,
     +  YMX, ERROR, NMSUB, MISS, NSAMPL)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE SETS THE PLOT LIMITS FOR VERTICAL PLOTS
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
      REAL
     +   YLB,YMN,YMX,YUB
      INTEGER
     +   IYM,M,N,NSAMPL
      LOGICAL
     +   ERROR,MISS
C
C  ARRAY ARGUMENTS
      REAL
     +   YM(IYM,M),YMMISS(M)
      CHARACTER
     +   NMSUB(6)*1
C
C  LOCAL SCALARS
      INTEGER
     +   I,II,IPRT,J
      LOGICAL
     +   HEAD,SETLMT
C
C  EXTERNAL FUNCTIONS
      LOGICAL
     +   MVCHK
      EXTERNAL MVCHK
C
C  EXTERNAL SUBROUTINES
      EXTERNAL ADJLMT,EHDR,IPRINT
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MAX,MIN
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     LOGICAL ERROR
C        A VALUE INDICATING WHETHER AN ERROR WAS DETECTED (TRUE)
C        OR NOT (FALSE).
C     LOGICAL HEAD
C        PRINT HEADING (HEAD=TRUE) OR NOT (HEAD=FALSE).
C     INTEGER I, II
C        INDEXING VARIABLES.
C     INTEGER IPRT
C        *
C     INTEGER IYM
C        ACTUAL ROW DIMENSION OF YM DECLARED IN THE USERS MAIN PROGRAM
C     INTEGER J
C        AN INDEX VARIABLE.
C     INTEGER M
C        THE NUMBER OF VECTORS IN YM
C     LOGICAL MISS
C        INDICATOR VARIABLE USED TO DESIGNATE WHETHER MISSING VALUES
C        MAY BE PRESENT (MISS = .TRUE.) OR NOT (MISS = .FALSE.)
C     INTEGER N
C        THE INTEGER NUMBER OF OBSERVATIONS .
C     CHARACTER*1 NMSUB(6)
C        THE CHARACTERS OF THE CALLING ROUTINES NAME.
C     INTEGER NSAMPL
C        *
C     LOGICAL SETLMT
C        AN INDICATOR VARIABLE USED TO DETERMINE IF STARTING VALUES
C        FOR XLB, XINC, YMN, YMX HAVE BEEN FOUND.
C     REAL YLB
C        THE USER SUPPLIED Y-AXIS LOWER BOUND.
C     REAL YM(IYM,M)
C        THE ARRAY CONTAINING THE DEPENDENT VARIABLE(S).
C     REAL YMMISS(M)
C        THE USER SUPPLIED CODE WHICH IS USED TO DETERMINE WHETHER OR
C        NOT AN OBSERVATION IS MISSING.
C        IF YM(I,J) = YMMISS(J), THE VALUE IS ASSUMED MISSING, OTHERWISE
C        IT IS NOT.
C     REAL YMN, YMX
C        THE Y-AXIS LOWER AND UPPER LIMITS ACTUALLY USED.
C     REAL YUB
C        THE USER SUPPLIED Y-AXIS UPPER BOUNDS.
C
      ERROR = .FALSE.
C
      IF (YLB.LT.YUB) THEN
C
C       SET LIMITS TO USER SPECIFIED VALUES
C
        YMN = YLB
        YMX = YUB
C
      ELSE
C
C       SET LIMITS TO RANGE OF VALUES WITHIN ANY USER SPECIFIED VALUES
C
        SETLMT = .FALSE.
        II = 1
C
C       FIND FIRST VALUE TO BE PLOTTED
C
        DO 20 I=1,N,NSAMPL
           DO 10 J=1,M
              IF (MISS .AND. MVCHK(YM(I,J),YMMISS(J))) GO TO 10
              IF ((YLB.LT.YUB) .AND. ((YM(I,J).LT.YLB) .OR.
     +           (YUB.LT.YM(I,J)))) GO TO 10
              IF (SETLMT) THEN
                YMN = MIN(YMN, YM(I,J))
                YMX = MAX(YMX, YM(I,J))
              ELSE
                YMN = YM(I,J)
                YMX = YM(I,J)
                SETLMT = .TRUE.
                II = I + NSAMPL
              END IF
   10      CONTINUE
           IF (SETLMT) GO TO 30
   20   CONTINUE
C
   30   IF (II.LE.1) THEN
C
C         NO VALUES TO BE PLOTTED.  PRINT ERROR MESSAGE
C
          ERROR = .TRUE.
          CALL IPRINT(IPRT)
          HEAD = .TRUE.
          CALL EHDR(NMSUB,HEAD)
          IF (YLB.GE.YUB) THEN
            WRITE (IPRT, 1010)
          ELSE
            WRITE (IPRT, 1020)
          END IF
          WRITE (IPRT, 1030)
C
        ELSE
C
C         FIND LIMITS FROM REMAINING VALUES
C
          IF (II.LE.N) THEN
            DO 50 I=II,N,NSAMPL
               DO 40 J=1,M
                  IF (MISS .AND. MVCHK(YM(I,J),YMMISS(J))) GO TO 40
                  IF ((YLB.LT.YUB) .AND. ((YM(I,J).LT.YLB) .OR.
     +               (YUB.LT.YM(I,J)))) GO TO 40
                  YMN = MIN(YMN, YM(I,J))
                  YMX = MAX(YMX, YM(I,J))
   40          CONTINUE
   50       CONTINUE
          END IF
        END IF
C
C       ADJUST Y AXIS LIMITS IF EQUAL
C
        IF (YMN .GE. YMX) CALL ADJLMT(YMN, YMX)
C
      END IF
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1010 FORMAT (/
     +   44H NO NON-MISSING PLOT COORDINATES WERE FOUND.)
 1020 FORMAT (/
     +   40H NO NON-MISSING VALUES WERE FOUND WITHIN,
     +   26H THE USER SUPPLIED LIMITS.)
 1030 FORMAT (/
     +   30H THE PLOT HAS BEEN SUPPRESSED.)
      END
