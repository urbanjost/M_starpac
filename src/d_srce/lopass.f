*LOPASS
      SUBROUTINE LOPASS (Y, N, FC, K, HLP, YF, NYF)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS SUBROUTINE CARRIES OUT LOW-PASS FILTERING OF THE
C     SERIES.  THE FILTER IS THE K-TERM
C     LEAST SQUARES APPROXIMATION TO THE CUTOFF FILTER
C     WITH CUTOF FREQUENCY FC.  ITS TRANSFER FUNCTION
C     HAS A TRANSITION BAND OF WIDTH DELTA SURROUNDING FC,
C     WHERE DELTA = 4*PI/K.
C
C     WRITTEN BY  -  PETER BLOOMFIELD
C                    FOURIER ANALYSIS OF TIME SERIES- AN
C                       INTRODUCTION
C                    JOHN WILEY AND SONS, NEW YORK, 1976
C                    PAGE 149
C     ADAPTED FOR STARPAC BY  -  JANET R. DONALDSON
C                                STATISTICAL ENGINEERING DIVISION
C                                NATIONAL BUREAU OF STANDARDS
C                                BOULDER, COLORADO
C
C     CREATION DATE  -  NOVEMBER 26, 1980
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      DOUBLE PRECISION
     +   FC
      INTEGER
     +   K,N,NYF
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   HLP(*),Y(*),YF(*)
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      INTEGER
     +   IPRT
      LOGICAL
     +   ERR01,ERR02,ERR03,ERR04,ERR05,HEAD
C
C  LOCAL ARRAYS
      CHARACTER
     +   LFC(8)*1,LK(8)*1,LN(8)*1,NMSUB(6)*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL EISGE,EISII,ERIODD,ERSII,ERSLFS,FLTSL,IPRINT,LPFLT
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     LOGICAL ERR01, ERR02, ERR03, ERR04, ERR05
C        VALUE(S) INDICATING WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT
C        (FALSE).
C     DOUBLE PRECISION FC
C        THE USER SUPPLIED CUTOFF FREQUENCY.
C     LOGICAL HEAD
C        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
C        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
C        OF HEAD WILL BE CHANGED TO FALSE.
C     DOUBLE PRECISION HLP(K)
C        THE ARRAY IN WHICH THE -IDEAL- HIGH PASS FILTER COEFFICIENTS
C        WILL BE RETURNED.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .EQ. 1, ERRORS HAVE BEEN DETECTED.
C     INTEGER IPRT
C        THE UNIT NUMBER USED FOR OUTPUT.
C     INTEGER K
C        THE NUMBER OF FILTER TERMS TO BE COMPUTED.
C     CHARACTER*1 LFC(8), LK(8), LN(8)
C        THE ARRAY CONTAINING THE NAMES OF THE VARIABLES FC, K AND N.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS IN THE SERIES Y.
C     CHARACTER*1 NMSUB(6)
C        THE ARRAY CONTAINING THE NAME OF THIS SUBROUTINE.
C     INTEGER NYF
C        THE NUMBER OF OBSERVATIONS IN THE FILTERED SERIES YF.
C     DOUBLE PRECISION Y(N)
C        THE VECTOR CONTAINING THE OBSERVED TIME SERIES.
C     DOUBLE PRECISION YF(N)
C        THE VECTOR IN WHICH THE FILTERED SERIES IS RETURNED.
C
C     SET UP NAME ARRAYS
C
      DATA
     +  NMSUB(1),  NMSUB(2),  NMSUB(3),  NMSUB(4),  NMSUB(5),  NMSUB(6)
     + /     'L',       'O',       'P',       'A',       'S',       'S'/
      DATA
     +  LFC(1), LFC(2), LFC(3), LFC(4), LFC(5), LFC(6), LFC(7), LFC(8)
     + /  'F',   'C',   ' ',   ' ',   ' ',   ' ',   ' ',   ' '/
      DATA
     +  LK(1), LK(2), LK(3), LK(4), LK(5), LK(6), LK(7), LK(8)
     + /  'K',   ' ',   ' ',   ' ',   ' ',   ' ',   ' ',   ' '/
      DATA
     +  LN(1), LN(2), LN(3), LN(4), LN(5), LN(6), LN(7), LN(8)
     + /  'N',   ' ',   ' ',   ' ',   ' ',   ' ',   ' ',   ' '/
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
      CALL ERSII(NMSUB, LFC, FC, 0.0D0, 0.5D0, 2, HEAD, ERR02, LFC, LFC)
C
      CALL EISII(NMSUB, LK, K, 1, N, 2, HEAD, ERR03, LK, LK)
C
      CALL ERIODD(NMSUB, LK, K, 1, HEAD, ERR04)
C
      IF (ERR01 .OR. ERR02 .OR. ERR03 .OR. ERR04) GO TO 10
C
      CALL ERSLFS(NMSUB, FC, K, HEAD, ERR05)
C
      IF (.NOT. ERR05) GO TO 20
C
   10 CONTINUE
      IERR = 1
      CALL IPRINT (IPRT)
      WRITE (IPRT, 1000)
      RETURN
C
   20 CONTINUE
C
      CALL LPFLT (FC, K, HLP)
C
      CALL FLTSL (Y, N, K, HLP, YF, NYF)
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (/42H THE CORRECT FORM OF THE CALL STATEMENT IS//
     +   46H       CALL LOPASS (Y, N, FC, K, HLP, YF, NYF))
      END
