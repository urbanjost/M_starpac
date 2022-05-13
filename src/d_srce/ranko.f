*RANKO
      SUBROUTINE RANKO(N, Y, H, R, T)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     VERSION  45.0    RANKO    3/ 6/70
C     *****
C     PUTS RANK OF N X"S IN VECTOR R. VECTOR H IS USED FOR STORAGE.
C     X,H AND R MUST BE DIMENSIONED N OR GREATER.
C     STORES CORRECTION FOR TIES IN T = SUM(T-1)*T*(T+1).
C        N.B.  T IS 12 TIMES VALUE COMPUTED BY ORIGINAL OMNITAB ROUTINE.
C     T=0  MEANS NO TIES.
C     WRITTEN BY DAVID HOGBEN, SEL, NBS.   4/9/69.
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      DOUBLE PRECISION
     +   T
      INTEGER
     +   N
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   R(N),Y(N)
      INTEGER
     +   H(N)
C
C  LOCAL SCALARS
      INTEGER
     +   I,IJ,J,K,K2
C
C  EXTERNAL SUBROUTINES
      EXTERNAL SRTIR,SRTRI
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     INTEGER H(N)
C        THE INDICES TO THE HIERARCHY OF R
C     INTEGER I
C        INDEX VARIABLE
C     INTEGER IJ
C        INDEX VARIABLE BASED ON I-J
C     INTEGER J
C        INDEX VARIABLE
C     INTEGER K
C        INDEX VARIABLE
C     INTEGER K2
C          INDEX VARIABLE
C     INTEGER N
C        NUMBER OF OBSERVATIONS
C     DOUBLE PRECISION R(N)
C        FINAL VECTOR CONTAINING RANK
C     DOUBLE PRECISION T
C        12 TIMES THE OMNITAB CORRECTION FOR TIES
C             T = SUM(T-1)*T*(T+1)
C             T = 0 MEANS NO TIES
C     DOUBLE PRECISION Y(N)
C        VECTOR TO BE RANKED
C
C
C     MOVE Y TO R AND PUT I IN H
C
      DO 10 I=1,N
         H(I) = I
         R(I) = Y(I)
   10 CONTINUE
C
C     SORT Y IN R, CARRY ALONG I IN H TO OBTAIN HIERARCHY IN H.
C
      CALL SRTIR(H, N, R)
C
C     REPLACE R(I) BY I*.
C     LET K BE SUCH THAT R(I)=R(I-J+1),J=1,K. THEN I* = I-(K-1)/2.
C
      K = 1
      T = 0
      DO 40 I=2,N
         IF (R(I).EQ.R(I-1)) THEN
            K = K + 1
         ELSE
            DO 30 J=1,K
               IJ = I - J
               R(IJ) = (I-1) - (K-1)/2.0D0
   30       CONTINUE
            T = T + (K-1)*K*(K+1)
            K = 1
         END IF
   40 CONTINUE
      T = T + (K-1)*K*(K+1)
      DO 50 I=1,K
         K2 = N + 1 - I
         R(K2) = N - (K-1)/2.0D0
   50 CONTINUE
C
C     SORT H CARRY ALONG R TO OBTAIN RANKS IN R
C
      CALL SRTRI(R, N, H)
      RETURN
      END
