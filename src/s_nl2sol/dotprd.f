*DOTPRD
      REAL FUNCTION DOTPRD(P, X, Y)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C  ***  RETURN THE INNER PRODUCT OF THE P-VECTORS X AND Y.  ***
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   P
C
C  ARRAY ARGUMENTS
      REAL
     +   X(*),Y(*)
C
C  LOCAL SCALARS
      REAL
     +   ONE,SQTETA,T,ZERO
      INTEGER
     +   I
C
C  EXTERNAL FUNCTIONS
      REAL
     +   RMDCON
      EXTERNAL RMDCON
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,MAX
C
C     INTEGER P
C     REAL X(*), Y(*)
C
C     INTEGER I
C     REAL ONE, SQTETA, T, ZERO
C/+
C     REAL MAX, ABS
C/
C     EXTERNAL RMDCON
C     REAL RMDCON
C
C  ***  RMDCON(2) RETURNS A MACHINE-DEPENDENT CONSTANT, SQTETA, WHICH
C  ***  IS SLIGHTLY LARGER THAN THE SMALLEST POSITIVE NUMBER THAT
C  ***  CAN BE SQUARED WITHOUT UNDERFLOWING.
C
      DATA ONE/1.0E0/, SQTETA/0.0E0/, ZERO/0.0E0/
C
      DOTPRD = ZERO
      IF (P .LE. 0) GO TO 999
      IF (SQTETA .EQ. ZERO) SQTETA = RMDCON(2)
      DO 20 I = 1, P
         T = MAX(ABS(X(I)), ABS(Y(I)))
         IF (T .GT. ONE) GO TO 10
         IF (T .LT. SQTETA) GO TO 20
         T = (X(I)/SQTETA)*Y(I)
         IF (ABS(T) .LT. SQTETA) GO TO 20
 10      DOTPRD = DOTPRD + X(I)*Y(I)
 20   CONTINUE
C
 999  RETURN
C  ***  LAST CARD OF DOTPRD FOLLOWS  ***
      END
