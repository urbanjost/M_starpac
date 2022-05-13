*V2NORM
      REAL FUNCTION V2NORM(P, X)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C  ***  RETURN THE 2-NORM OF THE P-VECTOR X, TAKING  ***
C  ***  CARE TO AVOID THE MOST LIKELY UNDERFLOWS.    ***
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
     +   X(*)
C
C  LOCAL SCALARS
      REAL
     +   ONE,R,SCALE,SQTETA,T,XI,ZERO
      INTEGER
     +   I,J
C
C  EXTERNAL FUNCTIONS
      REAL
     +   RMDCON
      EXTERNAL RMDCON
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,SQRT
C
C
      DATA ONE/1.0E0/, SQTETA/0.0E0/, ZERO/0.0E0/
C
      IF (P .GT. 0) GO TO 10
         V2NORM = ZERO
         GO TO 999
 10   DO 20 I = 1, P
         IF (X(I) .NE. ZERO) GO TO 30
 20      CONTINUE
      V2NORM = ZERO
      GO TO 999
C
 30   SCALE = ABS(X(I))
      IF (I .LT. P) GO TO 40
         V2NORM = SCALE
         GO TO 999
 40   T = ONE
      IF (SQTETA .EQ. ZERO) SQTETA = RMDCON(2)
C
C     ***  SQTETA IS (SLIGHTLY LARGER THAN) THE SQUARE ROOT OF THE
C     ***  SMALLEST POSITIVE FLOATING POINT NUMBER ON THE MACHINE.
C     ***  THE TESTS INVOLVING SQTETA ARE DONE TO PREVENT UNDERFLOWS.
C
      J = I + 1
      DO 60 I = J, P
         XI = ABS(X(I))
         IF (XI .GT. SCALE) GO TO 50
              R = XI / SCALE
              IF (R .GT. SQTETA) T = T + R*R
              GO TO 60
 50           R = SCALE / XI
              IF (R .LE. SQTETA) R = ZERO
              T = ONE  +  T * R*R
         SCALE = XI
 60      CONTINUE
C
      V2NORM = SCALE * SQRT(T)
 999  RETURN
C  ***  LAST CARD OF V2NORM FOLLOWS  ***
      END
