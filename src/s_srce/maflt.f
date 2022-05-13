*MAFLT
      SUBROUTINE MAFLT (Y, N, K, YF, NYF)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS SUBROUTINE PERFORMS A SIMPLE MOVING AVERAGE FILTERING
C     OPERATION ON AN INPUT SERIES Y, RETURNING THE FILTERED SERIES
C     IN YF.
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  NOVEMBER 26, 1980
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   K,N,NYF
C
C  ARRAY ARGUMENTS
      REAL
     +   Y(*),YF(*)
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      REAL
     +   HMA
      INTEGER
     +   IPRT
      LOGICAL
     +   ERR01,ERR02,ERR03,HEAD
C
C  LOCAL ARRAYS
      CHARACTER
     +   LK(8)*1,LN(8)*1,LONE(8)*1,NMSUB(6)*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL EISGE,EISII,ERIODD,FLTMA,IPRINT
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     LOGICAL ERR01, ERR02, ERR03
C        VALUE(S) INDICATING WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT
C        (FALSE).
C     LOGICAL HEAD
C        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
C        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
C        OF HEAD WILL BE CHANGED TO FALSE.
C     REAL HMA
C        THE VALUE OF EACH OF THE SIMPLE MOVING AVERAGE LINEAR FILTER
C        COEFFICIENTS.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .EQ. 1, ERRORS HAVE BEEN DETECTED.
C     INTEGER IPRT
C        THE UNIT NUMBER USED FOR OUTPUT.
C     INTEGER K
C        THE NUMBER OF FILTER TERMS.
C     CHARACTER*1 LK(8), LN(8), LONE(8)
C        THE ARRAYS CONTAINING THE NAMES OF THE VARIABLES K AND N.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS IN THE SERIES Y.
C     CHARACTER*1 NMSUB(6)
C        THE ARRAY CONTAINING THE NAME OF THIS SUBROUTINE.
C     INTEGER NYF
C        THE NUMBER OF OBSERVATIONS IN THE FILTERED SERIES YF.
C     REAL Y(N)
C        THE VECTOR CONTAINING THE OBSERVED TIME SERIES.
C     REAL YF(N)
C        THE VECTOR IN WHICH THE FILTERED SERIES IS RETURNED.
C
C     SET UP NAME ARRAYS
C
      DATA
     +  NMSUB(1),  NMSUB(2),  NMSUB(3),  NMSUB(4),  NMSUB(5),  NMSUB(6)
     + /     'M',       'A',       'F',       'L',       'T',       ' '/
      DATA
     +  LK(1), LK(2), LK(3), LK(4), LK(5), LK(6), LK(7), LK(8)
     + /  'K', ' ', ' ', ' ', ' ', ' ', ' ', ' '/
      DATA
     +  LN(1), LN(2), LN(3), LN(4), LN(5), LN(6), LN(7), LN(8)
     + /  'N', ' ', ' ', ' ', ' ', ' ', ' ', ' '/
      DATA
     +  LONE(1), LONE(2), LONE(3), LONE(4), LONE(5), LONE(6), LONE(7),
     +  LONE(8)  /  ' ', ' ', 'O', 'N', 'E', ' ', ' ', ' '/
C
C     SET UP FOR ERROR CHECKING
C
      IERR = 0
      HEAD = .TRUE.
C
C     CALL ERROR CHECKING ROUTINES
C
      CALL EISGE(NMSUB, LN, N, 3, 1, HEAD, ERR01, LN)
C
      CALL EISII(NMSUB, LK, K, 1, N, 1, HEAD, ERR02, LONE, LN)
C
      CALL ERIODD(NMSUB, LK, K, 1, HEAD, ERR03)
C
      IF (ERR01 .OR. ERR02 .OR. ERR03) GO TO 10
      GO TO 20
C
   10 CONTINUE
      IERR = 1
      CALL IPRINT (IPRT)
      WRITE (IPRT, 1000)
      RETURN
C
   20 CONTINUE
C
C     COMPUTE THE SIMPLE MOVING AVERAGE COEFFICIENTS
C
      HMA = K
      HMA = 1.0E0/HMA
C
      CALL FLTMA (Y, N, K, HMA, YF, NYF)
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (/42H THE CORRECT FORM OF THE CALL STATEMENT IS//
     +   '       CALL MAFLT (Y, N, K, YF, NYF)')
      END
