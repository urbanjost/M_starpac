*SLUPDT
      SUBROUTINE SLUPDT(A, COSMIN, P, SIZE, STEP, U, W, WCHMTD, WSCALE,
     +                  Y)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C
C  ***  UPDATE SYMMETRIC  A  SO THAT  A * STEP = Y  ***
C  ***  (LOWER TRIANGLE OF  A  STORED ROWWISE       ***
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      DOUBLE PRECISION
     +   COSMIN,SIZE,WSCALE
      INTEGER
     +   P
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   A(1),STEP(P),U(P),W(P),WCHMTD(P),Y(P)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   DENMIN,HALF,ONE,SDOTWM,T,UI,WI,ZERO
      INTEGER
     +   I,J,K
C
C  EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   DOTPRD,V2NORM
      EXTERNAL DOTPRD,V2NORM
C
C  EXTERNAL SUBROUTINES
      EXTERNAL SLVMUL
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,MIN
C
C  ***  PARAMETER DECLARATIONS  ***
C
C     INTEGER P
C     DOUBLE PRECISION A(1), COSMIN, SIZE, STEP(P), U(P), W(P),
C    1                 WCHMTD(P), WSCALE, Y(P)
C     DIMENSION A(P*(P+1)/2)
C
C  ***  LOCAL VARIABLES  ***
C
C     INTEGER I, J, K
C     DOUBLE PRECISION DENMIN, SDOTWM, T, UI, WI
C
C     ***  CONSTANTS  ***
C     DOUBLE PRECISION HALF, ONE, ZERO
C
C/
C  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  ***
C
C     EXTERNAL DOTPRD, SLVMUL, V2NORM
C     DOUBLE PRECISION DOTPRD, V2NORM
C
      DATA HALF/0.5D0/, ONE/1.0D0/, ZERO/0.0D0/
C
C-----------------------------------------------------------------------
C
      SDOTWM = DOTPRD(P, STEP, WCHMTD)
      DENMIN = COSMIN * V2NORM(P,STEP) * V2NORM(P,WCHMTD)
      WSCALE = ONE
      IF (DENMIN .NE. ZERO) WSCALE = MIN(ONE, ABS(SDOTWM/DENMIN))
      T = ZERO
      IF (SDOTWM .NE. ZERO) T = WSCALE / SDOTWM
      DO 10 I = 1, P
 10      W(I) = T * WCHMTD(I)
      CALL SLVMUL(P, U, A, STEP)
      T = HALF * (SIZE * DOTPRD(P, STEP, U)  -  DOTPRD(P, STEP, Y))
      DO 20 I = 1, P
 20      U(I) = T*W(I) + Y(I) - SIZE*U(I)
C
C  ***  SET  A = A + U*(W**T) + W*(U**T)  ***
C
      K = 1
      DO 40 I = 1, P
         UI = U(I)
         WI = W(I)
         DO 30 J = 1, I
              A(K) = SIZE*A(K) + UI*W(J) + WI*U(J)
              K = K + 1
 30           CONTINUE
 40      CONTINUE
C
      RETURN
C  ***  LAST CARD OF SLUPDT FOLLOWS  ***
      END
