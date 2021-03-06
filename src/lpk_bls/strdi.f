*STRDI
      SUBROUTINE STRDI(T,LDT,N,DET,JOB,INFO)
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER INFO,JOB,LDT,N
C
C  ARRAY ARGUMENTS
      REAL DET(2),T(LDT,*)
C
C  LOCAL SCALARS
      REAL TEMP,TEN
      INTEGER I,J,K,KB,KM1,KP1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL SAXPY,SSCAL
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,MOD
C
C
C     STRDI COMPUTES THE DETERMINANT AND INVERSE OF A REAL
C     TRIANGULAR MATRIX.
C
C     ON ENTRY
C
C        T       REAL(LDT,N)
C                T CONTAINS THE TRIANGULAR MATRIX. THE ZERO
C                ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND
C                THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE
C                USED TO STORE OTHER INFORMATION.
C
C        LDT     INTEGER
C                LDT IS THE LEADING DIMENSION OF THE ARRAY T.
C
C        N       INTEGER
C                N IS THE ORDER OF THE SYSTEM.
C
C        JOB     INTEGER
C                = 010       NO DET, INVERSE OF LOWER TRIANGULAR.
C                = 011       NO DET, INVERSE OF UPPER TRIANGULAR.
C                = 100       DET, NO INVERSE.
C                = 110       DET, INVERSE OF LOWER TRIANGULAR.
C                = 111       DET, INVERSE OF UPPER TRIANGULAR.
C
C     ON RETURN
C
C        T       INVERSE OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE UNCHANGED.
C
C        DET     REAL(2)
C                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE NOT REFERENCED.
C                DETERMINANT = DET(1) * 10.0**DET(2)
C                WITH  1.0 .LE. ABS(DET(1)) .LT. 10.0
C                OR  DET(1) .EQ. 0.0 .
C
C        INFO    INTEGER
C                INFO CONTAINS ZERO IF THE SYSTEM IS NONSINGULAR
C                AND THE INVERSE IS REQUESTED.
C                OTHERWISE INFO CONTAINS THE INDEX OF
C                A ZERO DIAGONAL ELEMENT OF T.
C
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS SAXPY,SSCAL
C     FORTRAN ABS,MOD
C
C     BEGIN BLOCK PERMITTING ...EXITS TO 180
C
C        COMPUTE DETERMINANT
C
         IF (JOB/100 .EQ. 0) GO TO 70
            DET(1) = 1.0E0
            DET(2) = 0.0E0
            TEN = 10.0E0
            DO 50 I = 1, N
               DET(1) = T(I,I)*DET(1)
C           ...EXIT
               IF (DET(1) .EQ. 0.0E0) GO TO 60
   10          IF (ABS(DET(1)) .GE. 1.0E0) GO TO 20
                  DET(1) = TEN*DET(1)
                  DET(2) = DET(2) - 1.0E0
               GO TO 10
   20          CONTINUE
   30          IF (ABS(DET(1)) .LT. TEN) GO TO 40
                  DET(1) = DET(1)/TEN
                  DET(2) = DET(2) + 1.0E0
               GO TO 30
   40          CONTINUE
   50       CONTINUE
   60       CONTINUE
   70    CONTINUE
C
C        COMPUTE INVERSE OF UPPER TRIANGULAR
C
         IF (MOD(JOB/10,10) .EQ. 0) GO TO 170
            IF (MOD(JOB,10) .EQ. 0) GO TO 120
C              BEGIN BLOCK PERMITTING ...EXITS TO 110
                  DO 100 K = 1, N
                     INFO = K
C              ......EXIT
                     IF (T(K,K) .EQ. 0.0E0) GO TO 110
                     T(K,K) = 1.0E0/T(K,K)
                     TEMP = -T(K,K)
                     CALL SSCAL(K-1,TEMP,T(1,K),1)
                     KP1 = K + 1
                     IF (N .LT. KP1) GO TO 90
                     DO 80 J = KP1, N
                        TEMP = T(K,J)
                        T(K,J) = 0.0E0
                        CALL SAXPY(K,TEMP,T(1,K),1,T(1,J),1)
   80                CONTINUE
   90                CONTINUE
  100             CONTINUE
                  INFO = 0
  110          CONTINUE
            GO TO 160
  120       CONTINUE
C
C              COMPUTE INVERSE OF LOWER TRIANGULAR
C
               DO 150 KB = 1, N
                  K = N + 1 - KB
                  INFO = K
C     ............EXIT
                  IF (T(K,K) .EQ. 0.0E0) GO TO 180
                  T(K,K) = 1.0E0/T(K,K)
                  TEMP = -T(K,K)
                  IF (K .NE. N) CALL SSCAL(N-K,TEMP,T(K+1,K),1)
                  KM1 = K - 1
                  IF (KM1 .LT. 1) GO TO 140
                  DO 130 J = 1, KM1
                     TEMP = T(K,J)
                     T(K,J) = 0.0E0
                     CALL SAXPY(N-K+1,TEMP,T(K,K),1,T(K,J),1)
  130             CONTINUE
  140             CONTINUE
  150          CONTINUE
               INFO = 0
  160       CONTINUE
  170    CONTINUE
  180 CONTINUE
      RETURN
      END
