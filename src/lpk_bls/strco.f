*STRCO
      SUBROUTINE STRCO(T,LDT,N,RCOND,Z,JOB)
C***BEGIN PROLOGUE  STRCO
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D2A3
C***KEYWORDS  CONDITION,FACTOR,LINEAR ALGEBRA,LINPACK,MATRIX,TRIANGULAR
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  ESTIMATES THE CONDITION OF A REAL TRIANGULAR MATRIX.
C***DESCRIPTION
C     STRCO ESTIMATES THE CONDITION OF A REAL TRIANGULAR MATRIX.
C     ON ENTRY
C        T       REAL(LDT,N)
C                T CONTAINS THE TRIANGULAR MATRIX.  THE ZERO
C                ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND
C                THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE
C                USED TO STORE OTHER INFORMATION.
C        LDT     INTEGER
C                LDT IS THE LEADING DIMENSION OF THE ARRAY T.
C        N       INTEGER
C                N IS THE ORDER OF THE SYSTEM.
C        JOB     INTEGER
C                = 0         T  IS LOWER TRIANGULAR.
C                = NONZERO   T  IS UPPER TRIANGULAR.
C     ON RETURN
C        RCOND   REAL
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  T .
C                FOR THE SYSTEM  T*X = B , RELATIVE PERTURBATIONS
C                IN  T  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  T  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C        Z       REAL(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  T  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  SASUM,SAXPY,SSCAL
C***END PROLOGUE  STRCO

C...SCALAR ARGUMENTS
      REAL RCOND
      INTEGER
     +   JOB,LDT,N

C...ARRAY ARGUMENTS
      REAL T(LDT,*),Z(*)

C...LOCAL SCALARS
      REAL EK,S,SM,TNORM,W,WK,WKM,YNORM
      INTEGER
     +   I1,J,J1,J2,K,KK,L
      LOGICAL
     +   LOWER

C...EXTERNAL FUNCTIONS
      REAL SASUM
      EXTERNAL
     +   SASUM

C...EXTERNAL SUBROUTINES
      EXTERNAL
     +   SAXPY,SSCAL

C...INTRINSIC FUNCTIONS
      INTRINSIC
     +   ABS,AMAX1,SIGN


C***FIRST EXECUTABLE STATEMENT  STRCO


      LOWER = JOB .EQ. 0

C     COMPUTE 1-NORM OF T

      TNORM = 0.0E0
      DO 10 J = 1, N
         L = J
         IF (LOWER) L = N + 1 - J
         I1 = 1
         IF (LOWER) I1 = J
         TNORM = AMAX1(TNORM,SASUM(L,T(I1,J),1))
   10 CONTINUE

C     RCOND = 1/(NORM(T)*(ESTIMATE OF NORM(INVERSE(T)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  T*Z = Y  AND  TRANS(T)*Y = E .
C     TRANS(T)  IS THE TRANSPOSE OF T .
C     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C     GROWTH IN THE ELEMENTS OF Y .
C     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.

C     SOLVE TRANS(T)*Y = E

      EK = 1.0E0
      DO 20 J = 1, N
         Z(J) = 0.0E0
   20 CONTINUE
      DO 100 KK = 1, N
         K = KK
         IF (LOWER) K = N + 1 - KK
         IF (Z(K) .NE. 0.0E0) EK = SIGN(EK,-Z(K))
         IF (ABS(EK-Z(K)) .LE. ABS(T(K,K))) GO TO 30
            S = ABS(T(K,K))/ABS(EK-Z(K))
            CALL SSCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = ABS(WK)
         SM = ABS(WKM)
         IF (T(K,K) .EQ. 0.0E0) GO TO 40
            WK = WK/T(K,K)
            WKM = WKM/T(K,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0E0
            WKM = 1.0E0
   50    CONTINUE
         IF (KK .EQ. N) GO TO 90
            J1 = K + 1
            IF (LOWER) J1 = 1
            J2 = N
            IF (LOWER) J2 = K - 1
            DO 60 J = J1, J2
               SM = SM + ABS(Z(J)+WKM*T(K,J))
               Z(J) = Z(J) + WK*T(K,J)
               S = S + ABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               W = WKM - WK
               WK = WKM
               DO 70 J = J1, J2
                  Z(J) = Z(J) + W*T(K,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0E0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)

      YNORM = 1.0E0

C     SOLVE T*Z = Y

      DO 130 KK = 1, N
         K = N + 1 - KK
         IF (LOWER) K = KK
         IF (ABS(Z(K)) .LE. ABS(T(K,K))) GO TO 110
            S = ABS(T(K,K))/ABS(Z(K))
            CALL SSCAL(N,S,Z,1)
            YNORM = S*YNORM
  110    CONTINUE
         IF (T(K,K) .NE. 0.0E0) Z(K) = Z(K)/T(K,K)
         IF (T(K,K) .EQ. 0.0E0) Z(K) = 1.0E0
         I1 = 1
         IF (LOWER) I1 = K + 1
         IF (KK .GE. N) GO TO 120
            W = -Z(K)
            CALL SAXPY(N-KK,W,T(I1,K),1,Z(I1),1)
  120    CONTINUE
  130 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0E0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
      YNORM = S*YNORM

      IF (TNORM .NE. 0.0E0) RCOND = YNORM/TNORM
      IF (TNORM .EQ. 0.0E0) RCOND = 0.0E0
      RETURN
      END
