*NRANDC
      SUBROUTINE NRANDC(Y, N, ISEED, YMEAN, SIGMA)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS SUBROUTINE GENERATES N NORMALLY DISTRIBUTED PSEUDO-
C     RANDOM NUMBERS WITH MEAN YMEAN AND STANDARD DEVIATION SIGMA.  THE
C     NUMBERS GENERATED ARE DETERMINED BY ISEED.  THEY ARE RETURNED IN Y
C
C     ORIGIN - CONCEIVED BY DR. PETER TRYON TO FACILITATE USE OF
C          EXISTING RANDOM NUMBER GENERATOR
C
C     WRITTEN BY -
C          JOHN E. KOONTZ AND JANET R. DONALDSON
C          STATISTICAL ENGINEERING DIVISION
C          NATIONAL BUREAU OF STANDARDS,
C          BOULDER, COLORADO 80302
C
C     CREATION DATE  -  DECEMBER 2, 1985
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      REAL
     +   SIGMA,YMEAN
      INTEGER
     +   ISEED,N
C
C  ARRAY ARGUMENTS
      REAL
     +   Y(*)
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      INTEGER
     +   I,IPRT,ISEEDU
      LOGICAL
     +   ERR01,ERR02,HEAD
C
C  LOCAL ARRAYS
      CHARACTER
     +   LN(8)*1,LONE(8)*1,LSIGMA(8)*1,LZERO(8)*1,NMSUB(6)*1
C
C  EXTERNAL FUNCTIONS
      REAL
     +   RANDN
      EXTERNAL RANDN
C
C  EXTERNAL SUBROUTINES
      EXTERNAL EISGE,EISRNG,ERSGE,IPRINT
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     LOGICAL ERR01, ERR02
C        VALUE(S) INDICATING WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT
C        (FALSE).
C     LOGICAL HEAD
C        A FLAG INDICATING WHETHER THE HEADING SOULD BE PRINTED
C        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
C        OF HEAD WILL BE CHANGED TO FALSE.
C     INTEGER I
C          THE INDEX OF THE COMPUTING LOOP
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THEIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .EQ. 1, ERRORS HAVE BEEN FOUND.
C     INTEGER IPRT
C        THE STANDARD OUTPUT FILE UNIT NUMBER
C     INTEGER ISEED
C        THE ISEED TO THE RANDOM NUMBER GENERATOR.
C        ISEED MUST LIE BETWEEN 0 AND 2**((MIN(32,I1MACH(8)+1))-1) -1,
C        INCLUSIVE.  IF ISEED IS NOT EQUAL TO 0, ISEED MUST BE ODD.
C     INTEGER ISEEDU
C        THE VALUE OF THE SEED ACTUALLY USED.
C     CHARACTER*1 LN(8), LONE(8), LSIGMA(8), LZERO(8)
C        THE ARRAY(S) CONTAINING THE NAME(S) OF THE VARIABLES(S) CHECKED
C        FOR ERRORS
C     INTEGER N
C        THE LENGTH OF DATA SET GENERATED
C     CHARACTER*1 NMSUB(6)
C        THE NAME OF THIS SUBROUTINE
C     REAL SIGMA
C        THE STANDARD DEVIATION OF THE GENERATED VALUES.
C     REAL Y(N)
C        THE GENERATED RANDOM VALUES.
C     REAL YMEAN
C        THE MEAN OF THE GENERATED VALUES.
C
      DATA  NMSUB(1), NMSUB(2), NMSUB(3), NMSUB(4), NMSUB(5), NMSUB(6)
     +   /       'N',      'R',      'A',      'N',      'D',      'C'/
      DATA     LN(1),    LN(2),    LN(3),    LN(4),    LN(5),    LN(6),
     +         LN(7),    LN(8)/'N', ' ', ' ', ' ', ' ', ' ', ' ', ' '/
      DATA   LONE(1),  LONE(2),  LONE(3),  LONE(4),  LONE(5),  LONE(6),
     +       LONE(7),  LONE(8)/'O', 'N', 'E', ' ', ' ', ' ', ' ', ' '/
      DATA LSIGMA(1),LSIGMA(2),LSIGMA(3),LSIGMA(4),LSIGMA(5),LSIGMA(6),
     +     LSIGMA(7),LSIGMA(8)/'S', 'I', 'G', 'M', 'A', ' ', ' ', ' '/
      DATA  LZERO(1), LZERO(2), LZERO(3), LZERO(4), LZERO(5), LZERO(6),
     +      LZERO(7), LZERO(8)/'Z', 'E', 'R', 'O', ' ', ' ', ' ', ' '/
C
      IERR = 0
C
      HEAD = .TRUE.
C
C     CHECK FOR INPUT ERRORS
C
      CALL EISGE(NMSUB, LN, N, 1, 2, HEAD, ERR01, LONE)
      CALL ERSGE(NMSUB, LSIGMA, SIGMA, 0.0E0, 2, HEAD, ERR02, LZERO)
      CALL EISRNG(NMSUB, ISEED, ISEEDU, HEAD)
C
      IF (ERR01.OR.ERR02) THEN
C
        CALL IPRINT(IPRT)
        WRITE (IPRT,1000)
        IERR = 1
C
      ELSE
C
C     GENERATE THE PSEUDO-RANDOM NUMBERS
C
        Y(1) = RANDN(ISEEDU)
        DO 20 I=1,N
           Y(I) = RANDN(0)*SIGMA + YMEAN
   20   CONTINUE
      END IF
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (/42H THE CORRECT FORM OF THE CALL STATEMENT IS//
     +   '       CALL NRANDC (Y, N, ISEED, YMEAN, SIGMA)')
      END
