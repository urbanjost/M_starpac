*TAPER
      SUBROUTINE TAPER (Y, N, TAPERP, YT)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE USER ROUTINE FOR APPLYING A SPLIT-COSINE-BELL
C     TAPER TO THE (CENTERED) OBSERVED SERIES Y, RETURNING THE TAPERED
C     SERIES IN YT.  THIS ROUTINE IS ADAPTED FROM BLOOMFIELDS
C     ROUTINE TAPER.
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  DECEMBER 7, 1981
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      REAL
     +   TAPERP
      INTEGER
     +   N
C
C  ARRAY ARGUMENTS
      REAL
     +   Y(*),YT(*)
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      REAL
     +   PI,WEIGHT
      INTEGER
     +   I,IPRT,J,M
      LOGICAL
     +   ERR01,HEAD
C
C  LOCAL ARRAYS
      CHARACTER
     +   LN(8)*1,NMSUB(6)*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL CENTER,EISGE,GETPI,IPRINT
C
C  INTRINSIC FUNCTIONS
      INTRINSIC COS,INT
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     LOGICAL ERR01
C        VALUE(S) INDICATING WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT
C        (FALSE).
C     LOGICAL HEAD
C        A VARIABLE USED TO INDICARE WHETHER A HEADING IS NEEDED FOR
C        ERROR MESSAGES (TRUE) OR NOT (FALSE).
C     INTEGER I
C        AN INDEXING VARIABLE.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .EQ. 1, ERRORS HAVE BEEN DETECTED.
C     INTEGER IPRT
C        THE LOGICAL UNIT NUMBER USED FOR OUTPUT.
C     INTEGER J
C        AN INDEXING VARIABLE.
C     CHARACTER*1 LN(8)
C        THE ARRAY(S) CONTAINING THE NAME(S) OF THE PARAMETER(S) CHECKED
C        FOR ERRORS.
C     INTEGER M
C        THE NUMBER OF POINTS AT EACH END OF THE SERIES TO BE
C        TAPERED.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS IN THE SERIES Y.
C     CHARACTER*1 NMSUB(6)
C        THE ARRAY CONTAINING THE NAME OF THIS SUBROUTINE.
C     REAL PI
C        THE VALUE OF PI.
C     REAL TAPERP
C        THE TOTAL PERCENTAGE OF THE DATA TO BE TAPERED.
C     REAL WEIGHT
C        THE ITH TAPER WEIGHT.
C     REAL Y(N)
C        THE VECTOR CONTAINING THE OBSERVED TIME SERIES.
C     REAL YT(N)
C        THE VECTOR IN WHICH THE TAPERED SERIES IS RETURNED.
C
C     SET UP NAME ARRAYS
C
      DATA
     +  NMSUB(1),  NMSUB(2),  NMSUB(3),  NMSUB(4),  NMSUB(5),  NMSUB(6)
     + /     'T',       'A',       'P',       'E',       'R',       ' '/
      DATA
     + LN(1), LN(2), LN(3), LN(4), LN(5), LN(6), LN(7), LN(8)
     + /'N',' ',' ',' ',' ',' ',' ',' '/
C
C     SET UP FOR ERROR CHECKING
C
      IERR = 0
      HEAD = .TRUE.
C
C     CALL ERROR CHECKING ROUTINES
C
      CALL EISGE(NMSUB, LN, N, 17, 1, HEAD, ERR01, LN)
      IF (.NOT. ERR01) GO TO 5
C
      IERR = 1
      CALL IPRINT (IPRT)
      WRITE (IPRT, 1000)
      RETURN
C
    5 CONTINUE
C
      CALL CENTER (Y, N, YT)
C
      IF ((TAPERP .LE. 0.0E0) .OR. (TAPERP .GT. 1.0E0)) RETURN
C
      CALL GETPI(PI)
C
      M = INT(TAPERP * N + 0.5E0) / 2
      IF (M .EQ. 0) RETURN
C
      DO 20 I = 1, M
         WEIGHT = 0.5E0 - 0.5E0 * COS(PI * (I-0.5E0) / M)
         YT(I) = WEIGHT * YT(I)
         J = N + 1 - I
         YT(J) = WEIGHT * YT(J)
   20 CONTINUE
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (/42H THE CORRECT FORM OF THE CALL STATEMENT IS//
     +   36H       CALL TAPER (Y, N, TAPERP, YT))
      END
