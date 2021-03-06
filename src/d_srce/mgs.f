*MGS
      SUBROUTINE MGS(A, B, N, NP, X, C, D, R, IR, IA, IER)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE COMPUTES THE SOLUTION  X  TO THE LINEAR SYSTEM OF
C     EQUATIONS  AX=B, USING THE METHOD OF MODIFIED GRAM-SCHMIDT.
C     THE MATRIX A IS DECOMPOSED INTO THREE MATRICES
C        Q  AN ORTHOGONAL MATRIX
C        D  A DIAGONAL MATRIX AND
C        R  AN UPPER TRIANGULAR MATRIX
C     THE SOLUTION VECTOR X IS THE VECTOR WHICH SOLVES THE SYSTEM
C     OF EQUATIONS  RX = C
C     X, A, AND B ARE NOT PRESERVED ON OUTPUT
C
C     ADAPTED FROM OMNITAB II BY -
C                  JANET R. DONALDSON
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  MAY 17, 1982
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   IA,IER,IR,N,NP
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   A(IA,NP),B(N),C(NP),D(NP),R(IR,NP),X(NP)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   SM1,SM2
      INTEGER
     +   I,J,JJ,K,NPJJMJ
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION A(IA,NP)
C           THE COEFFICIENTS MATRIX (NOT PRESERVED ON OUTPUT)
C     DOUBLE PRECISION B(N)
C           THE CONSTANT COLUMN MATRIX OF THE SYSTEM (NOT PRESERVED
C           ON OUTPUT)
C     DOUBLE PRECISION C(NP)
C           THE MATRIX C DESCRIBED ABOVE
C     DOUBLE PRECISION D(NP)
C           THE DIAGONAL ELEMENTS OF THE MATRIX D DESCRIBED ABOVE
C     INTEGER I
C           *
C     INTEGER IA
C           THE ROW DIMENSION OF A.
C     INTEGER IER
C           *
C     INTEGER IR
C           THE ROW DIMENSION OF R.
C     INTEGER J
C           *
C     INTEGER JJ
C           *
C     INTEGER K
C           *
C     INTEGER N
C           THE NUMBER OF OBSERVATIONS
C     INTEGER NP
C           THE NUMBER OF PARAMETERS
C     INTEGER NPJJMJ
C           *
C     DOUBLE PRECISION R(IR,NP)
C           THE UPPER ELEMENTS OF THE MATRIX R DESCRIBED ABOVE
C     DOUBLE PRECISION SM1
C           *
C     DOUBLE PRECISION SM2
C           *
C     DOUBLE PRECISION X(NP)
C           THE SOLUTION MATRIX
C
C
      IER = 0
C
      SM1 = 0.0D0
      SM2 = 0.0D0
      DO 10 I=1,N
         SM1 = A(I,1)*A(I,1) + SM1
         SM2 = A(I,1)*B(I) + SM2
   10 CONTINUE
      IF (SM1.EQ.0.0D0) GO TO 100
      D(1) = SM1
      C(1) = SM2/SM1
      IF (NP.EQ.1) GO TO 70
      DO 60 K=2,NP
         DO 40 J=K,NP
            SM1 = 0.0D0
            DO 20 I=1,N
               SM1 = A(I,K-1)*A(I,J) + SM1
   20       CONTINUE
            R(K-1,J) = SM1/D(K-1)
            DO 30 I=1,N
               A(I,J) = A(I,J) - A(I,K-1)*R(K-1,J)
   30       CONTINUE
   40    CONTINUE
         SM1 = 0.0D0
         SM2 = 0.0D0
         DO 50 I=1,N
            B(I) = B(I) - A(I,K-1)*C(K-1)
            SM1 = A(I,K)*A(I,K) + SM1
            SM2 = A(I,K)*B(I) + SM2
   50    CONTINUE
         IF (SM1.EQ.0.0D0) GO TO 100
         D(K) = SM1
         C(K) = SM2/SM1
   60 CONTINUE
C
C     COMPLETE BACKSOLVE
C
   70 X(NP) = C(NP)
      IF (NP.EQ.1) RETURN
      DO 90 I=2,NP
         K = NP + 1 - I
         JJ = K + 1
         SM1 = 0.0D0
         DO 80 J=JJ,NP
            NPJJMJ = NP + JJ - J
            SM1 = R(K,NPJJMJ)*X(NPJJMJ) + SM1
   80    CONTINUE
         X(K) = C(K) - SM1
   90 CONTINUE
      RETURN
  100 IER = 1
      RETURN
      END
