*DUPDAT
      SUBROUTINE DUPDAT(D, IV, J, N, NN, P, V)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C  ***  UPDATE SCALE VECTOR D FOR NL2ITR (NL2SOL VERSION 2.2)  ***
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   N,NN,P
C
C  ARRAY ARGUMENTS
      REAL
     +   D(P),J(NN,P),V(1)
      INTEGER
     +   IV(1)
C
C  LOCAL SCALARS
      REAL
     +   SII,T,VDFAC,ZERO
      INTEGER
     +   D0,DFAC,DTYPE,I,JTOL0,JTOLI,NITER,S,S1
C
C  EXTERNAL FUNCTIONS
      REAL
     +   V2NORM
      EXTERNAL V2NORM
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MAX,SQRT
C
C  ***  PARAMETER DECLARATIONS  ***
C
C     INTEGER IV(1), N, NN, P
C     REAL D(P), J(NN,P), V(1)
C     DIMENSION IV(*), V(*)
C
C  ***  LOCAL VARIABLES  ***
C
C     INTEGER D0, I, JTOLI, S1
C     REAL SII, T, VDFAC
C
C     ***  CONSTANTS  ***
C     REAL ZERO
C
C/
C  ***  EXTERNAL FUNCTION  ***
C
C     EXTERNAL V2NORM
C     REAL V2NORM
C
C  ***  SUBSCRIPTS FOR IV AND V  ***
C
C     INTEGER DFAC, DTYPE, JTOL0, NITER, S
      DATA DFAC/41/, DTYPE/16/, JTOL0/86/, NITER/31/, S/53/
C
      DATA ZERO/0.0E0/
C
C-----------------------------------------------------------------------
C
      I = IV(DTYPE)
      IF (I .EQ. 1) GO TO 20
         IF (IV(NITER) .GT. 0) GO TO 999
C
 20   VDFAC = V(DFAC)
      D0 = JTOL0 + P
      S1 = IV(S) - 1
      DO 30 I = 1, P
         S1 = S1 + I
         SII = V(S1)
         T = V2NORM(N, J(1,I))
         IF (SII .GT. ZERO) T = SQRT(T*T + SII)
         JTOLI = JTOL0 + I
         D0 = D0 + 1
         IF (T .LT. V(JTOLI)) T = MAX(V(D0), V(JTOLI))
         D(I) = MAX(VDFAC*D(I), T)
 30      CONTINUE
C
 999  RETURN
C  ***  LAST CARD OF DUPDAT FOLLOWS  ***
      END
