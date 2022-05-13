*ARFLT
      SUBROUTINE ARFLT (Y, N, IAR, PHI, YF, NYF)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS SUBROUTINE PERFORMS THE AUTOREGRESSIVE FILTERING
C     OPERATION DEFINED BY PHI, RETURNING THE FILTERED SERIES
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
     +   IAR,N,NYF
C
C  ARRAY ARGUMENTS
      REAL
     +   PHI(*),Y(*),YF(*)
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      REAL
     +   YMEAN
      INTEGER
     +   I,IPRT
      LOGICAL
     +   ERR01,HEAD
C
C  LOCAL ARRAYS
      CHARACTER
     +   LN(8)*1,NMSUB(6)*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL AMEAN,EISGE,FLTAR,IPRINT
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
C        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
C        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
C        OF HEAD WILL BE CHANGED TO FALSE.
C     INTEGER I
C        AN INDEX VARIABLE.
C     INTEGER IAR
C        THE NUMBER OF FILTER COEFFICIENTS.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .EQ. 1, ERRORS HAVE BEEN DETECTED.
C     INTEGER IPRT
C        THE UNIT NUMBER USED FOR OUTPUT.
C     CHARACTER*1 LN(8)
C        THE ARRAY(S) CONTAINING THE NAME(S) OF INPUT PARAMETER(S)
C        CHECKED FOR ERRORS.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS IN THE SERIES Y.
C     CHARACTER*1 NMSUB(6)
C        THE ARRAY CONTAINING THE NAME OF THIS SUBROUTINE.
C     INTEGER NYF
C        THE NUMBER OF OBSERVATIONS IN THE FILTERED SERIES YF.
C     REAL PHI(IAR)
C        THE VECTOR CONTAINING THE FILTER COEFFICIENTS.
C     REAL Y(N)
C        THE VECTOR CONTAINING THE OBSERVED TIME SERIES.
C     REAL YF(N)
C        THE VECTOR IN WHICH THE FILTERED SERIES IS RETURNED.
C     REAL YMEAN
C        THE MEAN OF THE INPUT SERIES Y.
C
C     SET UP NAME ARRAYS
C
      DATA
     +  NMSUB(1),  NMSUB(2),  NMSUB(3),  NMSUB(4),  NMSUB(5),  NMSUB(6)
     + /     'A',       'R',       'F',       'L',       'T',       ' '/
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
C
      IF (.NOT. ERR01) GO TO 10
C
      IERR = 1
      CALL IPRINT (IPRT)
      WRITE (IPRT, 1000)
      RETURN
C
   10 CONTINUE
C
C     COMPUTE ARITHMETIC MEAN
C
      CALL AMEAN(Y, N, YMEAN)
C
      DO 20 I = 1, N
         YF(I) = Y(I) - YMEAN
   20 CONTINUE
C
      CALL FLTAR (YF, N, IAR, PHI, YF, NYF)
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (/42H THE CORRECT FORM OF THE CALL STATEMENT IS//
     +   43H       CALL ARFLT (Y, N, IAR, PHI, YF, NYF))
      END
