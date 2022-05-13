*GAMR
      REAL FUNCTION GAMR (X)
C JULY 1977 EDITION.  W. FULLERTON, C3, LOS ALAMOS SCIENTIFIC LAB.
C THIS ROUTINE, NOT GAMMA(X), SHOULD BE THE FUNDAMENTAL ONE.
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      REAL X
C
C  LOCAL SCALARS
      REAL ALNGX,SGNGX
      INTEGER IROLD
C
C  EXTERNAL FUNCTIONS
      REAL GAMMA
      EXTERNAL GAMMA
C
C  EXTERNAL SUBROUTINES
      EXTERNAL ALGAMS,XERCLR,XGETF,XSETF
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,AINT,EXP
C
      GAMR = 0.0
      IF (X.LE.0.0 .AND. AINT(X).EQ.X) RETURN
C
      CALL XGETF (IROLD)
      CALL XSETF (1)
      IF (ABS(X).GT.10.0) GO TO 10
      GAMR = 1.0/GAMMA(X)
      CALL XERCLR
      CALL XSETF (IROLD)
      RETURN
C
 10   CALL ALGAMS (X, ALNGX, SGNGX)
      CALL XERCLR
      CALL XSETF (IROLD)
      GAMR = SGNGX * EXP(-ALNGX)
      RETURN
C
      END
