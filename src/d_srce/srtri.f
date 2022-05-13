*SRTRI
      SUBROUTINE SRTRI(A, LA, IR)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C   FUNCTION     SRTRI  - SORT ARRAY A ON AN INTEGER ARRAY IR.
C                           IF THE INTEGER ARRAY IS A PERMUTATION
C                           VECTOR FOR THE ARRAY, THEN THE
C                           ARRAY IS RESTORED TO ITS ORIGINAL
C                           (UNPERMUTED) ORDER.
C                           PERMUTATIONS RETURNED
C   USAGE               - CALL SRTRI (A,LA,IR)
C   PARAMETERS   A(LA)  - ON INPUT, CONTAINS THE ARRAY TO BE SORTED
C                         ON OUTPUT, A CONTAINS THE SORTED ARRAY
C                LA     - INPUT VARIABLE CONTAINING THE NUMBER OF
C                           ELEMENTS IN THE ARRAY TO BE SORTED
C                IR(LA) - ON INPUT, CONTAINS THE INTEGER KEY ARRAY
C                         ON OUTPUT, CONTAINS THE SORTED KEY ARRAY
C                           1,2,...,LA.
C                       - THEN ON OUTPUT, IR CONTAINS A RECORD OF THE
C                           PERMUTATIONS MADE ON THE VECTOR A.
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   LA
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   A(LA)
      INTEGER
     +   IR(LA)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   R,T,TT
      INTEGER
     +   I,IJ,IT,ITT,J,K,L,M
C
C  LOCAL ARRAYS
      INTEGER
     +   IL(21),IU(21)
C
C  INTRINSIC FUNCTIONS
      INTRINSIC INT
C
C
      M = 1
      I = 1
      J = LA
      R = .375D0
   10 IF (I.EQ.J) GO TO 90
      IF (R.GT.0.5898437D0) GO TO 20
      R = R + 3.90625D-2
      GO TO 30
   20 R = R - .21875D0
   30 K = I
C                                  SELECT A CENTRAL ELEMENT OF THE
C                                  ARRAY AND SAVE IT IN LOCATION IT
      IJ = I + INT((J-I)*R)
      T = A(IJ)
      IT = IR(IJ)
C                                  IF FIRST ELEMENT OF ARRAY IS GREATER
C                                  THAN IT, INTERCHANGE WITH IT
      IF (IR(I).LE.IT) GO TO 40
      A(IJ) = A(I)
      A(I) = T
      T = A(IJ)
      IR(IJ) = IR(I)
      IR(I) = IT
      IT = IR(IJ)
   40 L = J
C                                  IF LAST ELEMENT OF ARRAY IS LESS THAN
C                                  IT, INTERCHANGE WITH IT
      IF (IR(J).GE.IT) GO TO 60
      A(IJ) = A(J)
      A(J) = T
      T = A(IJ)
      IR(IJ) = IR(J)
      IR(J) = IT
      IT = IR(IJ)
C                                  IF FIRST ELEMENT OF ARRAY IS GREATER
C                                  THAN IT, INTERCHANGE WITH IT
      IF (IR(I).LE.IT) GO TO 60
      A(IJ) = A(I)
      A(I) = T
      T = A(IJ)
      IR(IJ) = IR(I)
      IR(I) = IT
      IT = IR(IJ)
      GO TO 60
   50 TT = A(L)
      A(L) = A(K)
      A(K) = TT
      ITT = IR(L)
      IR(L) = IR(K)
      IR(K) = ITT
C                                  FIND AN ELEMENT IN THE SECOND HALF OF
C                                  THE ARRAY WHICH IS SMALLER THAN IT
   60 L = L - 1
      IF (IR(L).GT.IT) GO TO 60
C                                  FIND AN ELEMENT IN THE FIRST HALF OF
C                                  THE ARRAY WHICH IS GREATER THAN IT
   70 K = K + 1
      IF (IR(K).LT.IT) GO TO 70
C                                  INTERCHANGE THESE ELEMENTS
      IF (K.LE.L) GO TO 50
C                                  SAVE UPPER AND LOWER SUBSCRIPTS OF
C                                  THE ARRAY YET TO BE SORTED
      IF (L-I.LE.J-K) GO TO 80
      IL(M) = I
      IU(M) = L
      I = K
      M = M + 1
      GO TO 100
   80 IL(M) = K
      IU(M) = J
      J = L
      M = M + 1
      GO TO 100
C                                  BEGIN AGAIN ON ANOTHER PORTION OF
C                                  THE UNSORTED ARRAY
   90 M = M - 1
      IF (M.EQ.0) RETURN
      I = IL(M)
      J = IU(M)
  100 IF (J-I.GE.1) GO TO 30
      IF (I.EQ.1) GO TO 10
      I = I - 1
  110 I = I + 1
      IF (I.EQ.J) GO TO 90
      T = A(I+1)
      IT = IR(I+1)
      IF (IR(I).LE.IT) GO TO 110
      K = I
  120 A(K+1) = A(K)
      IR(K+1) = IR(K)
      K = K - 1
      IF (IT.LT.IR(K)) GO TO 120
      A(K+1) = T
      IR(K+1) = IT
      GO TO 110
      END
