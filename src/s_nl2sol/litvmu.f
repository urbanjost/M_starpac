*LITVMU
      SUBROUTINE LITVMU(N, X, L, Y)
C
C  ***  SOLVE  (L**T)*X = Y,  WHERE  L  IS AN  N X N  LOWER TRIANGULAR
C  ***  MATRIX STORED COMPACTLY BY ROWS.  X AND Y MAY OCCUPY THE SAME
C  ***  STORAGE.  ***
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   N
C
C  ARRAY ARGUMENTS
      REAL
     +   L(1),X(N),Y(N)
C
C  LOCAL SCALARS
      REAL
     +   XI,ZERO
      INTEGER
     +   I,I0,II,IJ,IM1,J,NP1
C
      DATA ZERO/0.0E0/
C
      DO 10 I = 1, N
 10      X(I) = Y(I)
      NP1 = N + 1
      I0 = N*(N+1)/2
      DO 30 II = 1, N
         I = NP1 - II
         XI = X(I)/L(I0)
         X(I) = XI
         IF (I .LE. 1) GO TO 999
         I0 = I0 - I
         IF (XI .EQ. ZERO) GO TO 30
         IM1 = I - 1
         DO 20 J = 1, IM1
              IJ = I0 + J
              X(J) = X(J) - XI*L(IJ)
 20           CONTINUE
 30      CONTINUE
 999  RETURN
C  ***  LAST CARD OF LITVMU FOLLOWS  ***
      END
