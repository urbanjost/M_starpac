*SRTIRR
      SUBROUTINE SRTIRR(IR, RR, LA, A)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE SORTS THE LENGTH LA ARRAY A, THE LENGTH LA
C     INTEGER ARRAY IR, AND THE LENGTH LA ARRAY RR INTO
C     ASCENDING ORDER, BASED ON THE VALUES IN A.  THE ARRAY
C     A CONSTITUTES THE SORTING KEY.  THE OTHER TWO ARRAYS ARE
C     CARRIED ALONG.  ORDINARILY THE ARRAY IR CONTAINS THE
C     VALUES 1, ..., LA INITIALLY, SO THAT THE THREE ARRAYS CAN
C     LATER BE SORTED AGAIN WITH IR AS THE KEY, IN ORDER TO
C     RESTORE A AND RR TO THEIR ORIGINAL ORDER.
C
C     WRITTEN BY - JOHN E. KOONTZ
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  MAY 17, 1982
C        (BASED CLOSELY ON THE IMSL CDC LIBRARY 3 ROUTINE VSORTP)
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
     +   A(LA),RR(LA)
      INTEGER
     +   IR(LA)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   R,RT,RTT,T,TT
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
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION A(LA)
C        INPUT/OUTPUT PARAMETER.  THE KEY ARRAY.
C     INTEGER I
C        *
C     INTEGER IJ
C        *
C     INTEGER IL(21)
C        *
C     INTEGER IR(LA)
C        INPUT/OUTPUT PARAMETER.  THE INTEGER ARRAY CARRIED ALONG
C        IN THE SORT.  INITIALLY IT SHOULD CONTAIN 1, ..., LA.
C        ON EXIT IT CONTAINS THE PERMUTATION VECTOR OF THE SORT.
C        SORTING ON THE PERMUTATION VECTOR WILL RESTORE THE KEY
C        ARRAY A AND THE ARRAY RR TO THEIR ORIGINAL ORDERS.
C     INTEGER IT
C        *
C     INTEGER ITT
C        *
C     INTEGER IU(21)
C        *
C     INTEGER J
C        *
C     INTEGER K
C        *
C     INTEGER L
C        *
C     INTEGER LA
C        INPUT PARAMETER.  THE LENGTH OF THE INPUT/OUTPUT PARAMETERS
C        A, IR, AND RR.
C     INTEGER M
C        *
C     DOUBLE PRECISION R
C        *
C     DOUBLE PRECISION RR(LA)
C        INPUT/OUTPUT PARAMETER.  THE ARRAY CARRIED ALONG IN
C        THE SORT.  IT MIGHT BE THE SET OF WEIGHTS FOR A.
C     DOUBLE PRECISION RT
C        *
C     DOUBLE PRECISION RTT
C        *
C     DOUBLE PRECISION T
C        *
C     DOUBLE PRECISION TT
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
C                                  ARRAY AND SAVE IT IN LOCATION T
      IJ = I + INT((J-I)*R)
      T = A(IJ)
      IT = IR(IJ)
      RT = RR(IJ)
C                                  IF FIRST ELEMENT OF ARRAY IS GREATER
C                                  THAN T, INTERCHANGE WITH T
      IF (A(I).LE.T) GO TO 40
      A(IJ) = A(I)
      A(I) = T
      T = A(IJ)
      IR(IJ) = IR(I)
      IR(I) = IT
      IT = IR(IJ)
      RR(IJ) = RR(I)
      RR(I) = RT
      RT = RR(IJ)
   40 L = J
C                                  IF LAST ELEMENT OF ARRAY IS LESS THAN
C                                  T, INTERCHANGE WITH T
      IF (A(J).GE.T) GO TO 60
      A(IJ) = A(J)
      A(J) = T
      T = A(IJ)
      IR(IJ) = IR(J)
      IR(J) = IT
      IT = IR(IJ)
      RR(IJ) = RR(J)
      RR(J) = RT
      RT = RR(IJ)
C                                  IF FIRST ELEMENT OF ARRAY IS GREATER
C                                  THAN T, INTERCHANGE WITH T
      IF (A(I).LE.T) GO TO 60
      A(IJ) = A(I)
      A(I) = T
      T = A(IJ)
      IR(IJ) = IR(I)
      IR(I) = IT
      IT = IR(IJ)
      RR(IJ) = RR(I)
      RR(I) = RT
      RT = RR(IJ)
      GO TO 60
   50 TT = A(L)
      A(L) = A(K)
      A(K) = TT
      ITT = IR(L)
      IR(L) = IR(K)
      IR(K) = ITT
      RTT = RR(L)
      RR(L) = RR(K)
      RR(K) = RTT
C                                  FIND AN ELEMENT IN THE SECOND HALF OF
C                                  THE ARRAY WHICH IS SMALLER THAN T
   60 L = L - 1
      IF (A(L).GT.T) GO TO 60
C                                  FIND AN ELEMENT IN THE FIRST HALF OF
C                                  THE ARRAY WHICH IS GREATER THAN T
   70 K = K + 1
      IF (A(K).LT.T) GO TO 70
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
      RT = RR(I+1)
      IF (A(I).LE.T) GO TO 110
      K = I
  120 A(K+1) = A(K)
      IR(K+1) = IR(K)
      RR(K+1) = RR(K)
      K = K - 1
      IF (T.LT.A(K)) GO TO 120
      A(K+1) = T
      IR(K+1) = IT
      RR(K+1) = RT
      GO TO 110
      END
