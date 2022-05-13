*FFT
      SUBROUTINE FFT(A, B, NTOT, N, NSPAN, ISN)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C  MULTIVARIATE COMPLEX FOURIER TRANSFORM, COMPUTED IN PLACE
C    USING MIXED-RADIX FAST FOURIER TRANSFORM ALGORITHM.
C  BY R. C. SINGLETON, STANFORD RESEARCH INSTITUTE, OCT. 1968
C  ARRAYS A AND B ORIGINALLY HOLD THE DOUBLE PRECISION AND IMAGINARY
C    COMPONENTS OF THE DATA, AND RETURN THE DOUBLE PRECISION AND
C    IMAGINARY COMPONENTS OF THE RESULTING FOURIER COEFFICIENTS.
C  MULTIVARIATE DATA IS INDEXED ACCORDING TO THE FORTRAN
C    ARRAY ELEMENT SUCCESSOR FUNCTION, WITHOUT LIMIT
C    ON THE NUMBER OF IMPLIED MULTIPLE SUBSCRIPTS.
C    THE SUBROUTINE IS CALLED ONCE FOR EACH VARIATE.
C    THE CALLS FOR A MULTIVARIATE TRANSFORM MAY BE IN ANY ORDER.
C  NTOT IS THE TOTAL NUMBER OF COMPLEX DATA VALUES.
C  N IS THE DIMENSION OF THE CURRENT VARIABLE.
C  NSPAN/N IS THE SPACING OF CONSECUTIVE DATA VALUES
C    WHILE INDEXING THE CURRENT VARIABLE.
C  THE SIGN OF ISN DETERMINES THE SIGN OF THE COMPLEX
C    EXPONENTIAL, AND THE MAGNITUDE OF ISN IS NORMALLY ONE.
C  A TRI-VARIATE TRANSFORM WITH A(N1,N2,N3), B(N1,N2,N3)
C    IS COMPUTED BY
C      CALL FFT(A,B,N1*N2*N3,N1,N1,1)
C      CALL FFT(A,B,N1*N2*N3,N2,N1*N2,1)
C      CALL FFT(A,B,N1*N2*N3,N3,N1*N2*N3,1)
C  FOR A SINGLE-VARIATE TRANSFORM,
C    NTOT = N = NSPAN = (NUMBER OF COMPLEX DATA VALUES), F.G.
C      CALL FFT(A,B,N,N,N,1)
C  THE DATA MAY ALTERNATIVELY BE STORED IN A SINGLE COMPLEX
C    ARRAY A, THEN THE MAGNITUDE OF ISN CHANGED TO TWO TO
C    GIVE THE CORRECT INDEXING INCREMENT AND A(2) USED TO
C    PASS THE INITIAL ADDRESS FOR THE SEQUENCE OF IMAGINARY
C    VALUES, E.G.
C      CALL FFT(A,A(2),NTOT,N,NSPAN,2)
C  ARRAYS AT(MAXF), CK(MAXF), BT(MAXF), SK(MAXF), AND NP(MAXP)
C    ARE USED FOR TEMPORARY STORAGE.  IF THE AVAILABEL STORAGE
C    IS INSUFFICIENT, THE PROGRAM IS TERMINATED BY A STOP.
C    MAXF MUST BE .GE. THE MAXIMUM PRIME FACTOR OF N.
C    MAXP MUST BE .GT. THE NUMBER OF PRIME FACTORS OF N.
C
C     NB. THE ABOVE DESCRIPTION OF MAXP APPEARS TO BE INCORRECT.
C         MAXP SEEMS TO BE THE MAXIMUM SIZE OF THE SQUARE FREE
C         PORTION K OF N.
C
C    IN ADDITION, IF THE SQUARE-FREE PORTION K OF N HAS TWO OR
C    MORE PRIME FACTORS, THEN MAXP MUST BE .GE. K-1.
C     DIMENSION A(1), B(1)
C  ARRAY STORAGE IN NFAC FOR A MAXIMUM OF 11 FACTORS OF N.
C  IF N HAS MORE THAN ONE SQUARE-FREE FACTOR, THE PRODUCT OF THE
C    SQUARE-FREE FACTORS MUST BE .LE. 210
C     DIMENSION NFAC(11), NP(209)
C  ARRAY STORAGE FOR MAXIMUM PRIME FACTOR OF 23
C     DIMENSION AT(23), CK(23), BT(23), SK(23)
C
C
C  VARIABLE DECLARATIONS
C
C  PARAMETERS
      INTEGER
     +   MAXF1
      PARAMETER (MAXF1=23)
      INTEGER
     +   MAXP1
      PARAMETER (MAXP1=209)
C
C  SCALAR ARGUMENTS
      INTEGER
     +   ISN,N,NSPAN,NTOT
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   A(*),B(*)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   AA,AJ,AJM,AJP,AK,AKM,AKP,BB,BJ,BJM,BJP,BK,BKM,BKP,C1,C2,C3,
     +   C72,CD,RAD,RADF,S1,S120,S2,S3,S72,SD
      INTEGER
     +   I,II,INC,IPRT,J,JC,JF,JJ,K,K1,K2,K3,K4,KK,KS,KSPAN,KSPNN,
     +   KT,M,MAXF,MAXP,NN,NT
C
C  LOCAL ARRAYS
      DOUBLE PRECISION
     +   AT(MAXF1),BT(MAXF1),CK(MAXF1),SK(MAXF1)
      INTEGER
     +   NFAC(11),NP(MAXP1)
C
C  EXTERNAL SUBROUTINES
      EXTERNAL IPRINT
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ATAN,COS,MOD,SIN,SQRT
C
C  EQUIVALENCES
      EQUIVALENCE (I,II)
C
C  THE FOLLOWING TWO CONSTANTS SHOULD AGREE WITH THE ARRAY DIMENSIONS.
      MAXF = MAXF1
      MAXP = MAXP1
      IF (N.LT.2) RETURN
C
C  INITIALIZE VARIABLES
C
      C1 = 0
      C2 = 0
      C3 = 0
      S1 = 0
      S2 = 0
      S3 = 0
      K1 = 0
      K2 = 0
      K3 = 0
      K4 = 0
C
      INC = ISN
      RAD = 8.0D0*ATAN(1.0D0)
      S72 = RAD/5.0D0
      C72 = COS(S72)
      S72 = SIN(S72)
      S120 = SQRT(0.75D0)
      IF (ISN.GE.0) GO TO 10
      S72 = -S72
      S120 = -S120
      RAD = -RAD
      INC = -INC
   10 NT = INC*NTOT
      KS = INC*NSPAN
      KSPAN = KS
      NN = NT - INC
      JC = KS/N
      RADF = RAD*JC*0.5D0
      I = 0
      JF = 0
C  DETERMINE THE FACTORS OF N
      M = 0
      K = N
      GO TO 30
   20 M = M + 1
      NFAC(M) = 4
      K = K/16
   30 IF (K-(K/16)*16.EQ.0) GO TO 20
      J = 3
      JJ = 9
      GO TO 50
   40 M = M + 1
      NFAC(M) = J
      K = K/JJ
   50 IF (MOD(K,JJ).EQ.0) GO TO 40
      J = J + 2
      JJ = J**2
      IF (JJ.LE.K) GO TO 50
      IF (K.GT.4) GO TO 60
      KT = M
      NFAC(M+1) = K
      IF (K.NE.1) M = M + 1
      GO TO 100
   60 IF (K-(K/4)*4.NE.0) GO TO 70
      M = M + 1
      NFAC(M) = 2
      K = K/4
   70 KT = M
      J = 2
   80 IF (MOD(K,J).NE.0) GO TO 90
      M = M + 1
      NFAC(M) = J
      K = K/J
   90 J = ((J+1)/2)*2 + 1
      IF (J.LE.K) GO TO 80
  100 IF (KT.EQ.0) GO TO 120
      J = KT
  110 M = M + 1
      NFAC(M) = NFAC(J)
      J = J - 1
      IF (J.NE.0) GO TO 110
C  COMPUTE FOURIER TRANSFORM
  120 SD = RADF/KSPAN
      CD = 2.0D0*SIN(SD)**2
      SD = SIN(SD+SD)
      KK = 1
      I = I + 1
      IF (NFAC(I).NE.2) GO TO 170
C  TRANSFORM FOR FACTOR OF 2 (INCLUDING ROTATION FACTOR)
      KSPAN = KSPAN/2
      K1 = KSPAN + 2
  130 K2 = KK + KSPAN
      AK = A(K2)
      BK = B(K2)
      A(K2) = A(KK) - AK
      B(K2) = B(KK) - BK
      A(KK) = A(KK) + AK
      B(KK) = B(KK) + BK
      KK = K2 + KSPAN
      IF (KK.LE.NN) GO TO 130
      KK = KK - NN
      IF (KK.LE.JC) GO TO 130
      IF (KK.GT.KSPAN) GO TO 360
  140 C1 = 1.0D0 - CD
      S1 = SD
  150 K2 = KK + KSPAN
      AK = A(KK) - A(K2)
      BK = B(KK) - B(K2)
      A(KK) = A(KK) + A(K2)
      B(KK) = B(KK) + B(K2)
      A(K2) = C1*AK - S1*BK
      B(K2) = S1*AK + C1*BK
      KK = K2 + KSPAN
      IF (KK.LT.NT) GO TO 150
      K2 = KK - NT
      C1 = -C1
      KK = K1 - K2
      IF (KK.GT.K2) GO TO 150
      AK = C1 - (CD*C1+SD*S1)
      S1 = (SD*C1-CD*S1) + S1
C  THE FOLLOWING THREE STATEMENTS COMPENSATE FOR TRUNCATION ERROR
      C1 = 0.5D0/(AK**2+S1**2) + 0.5D0
      S1 = C1*S1
      C1 = C1*AK
      KK = KK + JC
      IF (KK.LT.K2) GO TO 150
      K1 = K1 + INC + INC
      KK = (K1-KSPAN)/2 + JC
      IF (KK.LE.JC+JC) GO TO 140
      GO TO 120
C  TRANSFORM FOR FACTOR OF 3 (OPTIONAL CODE)
  160 K1 = KK + KSPAN
      K2 = K1 + KSPAN
      AK = A(KK)
      BK = B(KK)
      AJ = A(K1) + A(K2)
      BJ = B(K1) + B(K2)
      A(KK) = AK + AJ
      B(KK) = BK + BJ
      AK = -0.5D0*AJ + AK
      BK = -0.5D0*BJ + BK
      AJ = (A(K1)-A(K2))*S120
      BJ = (B(K1)-B(K2))*S120
      A(K1) = AK - BJ
      B(K1) = BK + AJ
      A(K2) = AK + BJ
      B(K2) = BK - AJ
      KK = K2 + KSPAN
      IF (KK.LT.NN) GO TO 160
      KK = KK - NN
      IF (KK.LE.KSPAN) GO TO 160
      GO TO 320
C  TRANSFORM FOR FACTOR OF 4
  170 IF (NFAC(I).NE.4) GO TO 260
      KSPNN = KSPAN
      KSPAN = KSPAN/4
  180 C1 = 1.0D0
      S1 = 0
  190 K1 = KK + KSPAN
      K2 = K1 + KSPAN
      K3 = K2 + KSPAN
      AKP = A(KK) + A(K2)
      AKM = A(KK) - A(K2)
      AJP = A(K1) + A(K3)
      AJM = A(K1) - A(K3)
      A(KK) = AKP + AJP
      AJP = AKP - AJP
      BKP = B(KK) + B(K2)
      BKM = B(KK) - B(K2)
      BJP = B(K1) + B(K3)
      BJM = B(K1) - B(K3)
      B(KK) = BKP + BJP
      BJP = BKP - BJP
      IF (ISN.LT.0) GO TO 220
      AKP = AKM - BJM
      AKM = AKM + BJM
      BKP = BKM + AJM
      BKM = BKM - AJM
      IF (S1.EQ.0.0D0) GO TO 230
  200 A(K1) = AKP*C1 - BKP*S1
      B(K1) = AKP*S1 + BKP*C1
      A(K2) = AJP*C2 - BJP*S2
      B(K2) = AJP*S2 + BJP*C2
      A(K3) = AKM*C3 - BKM*S3
      B(K3) = AKM*S3 + BKM*C3
      KK = K3 + KSPAN
      IF (KK.LE.NT) GO TO 190
  210 C2 = C1 - (CD*C1+SD*S1)
      S1 = (SD*C1-CD*S1) + S1
      C1 = 0.5D0/(C2**2+S1**2) + 0.5D0
      S1 = C1*S1
      C1 = C1*C2
      C2 = C1**2 - S1**2
      S2 = 2.0D0*C1*S1
      C3 = C2*C1 - S2*S1
      S3 = C2*S1 + S2*C1
      KK = KK - NT + JC
      IF (KK.LE.KSPAN) GO TO 190
      KK = KK - KSPAN + INC
      IF (KK.LE.JC) GO TO 180
      IF (KSPAN.EQ.JC) GO TO 360
      GO TO 120
  220 AKP = AKM + BJM
      AKM = AKM - BJM
      BKP = BKM - AJM
      BKM = BKM + AJM
      IF (S1.NE.0.0D0) GO TO 200
  230 A(K1) = AKP
      B(K1) = BKP
      A(K2) = AJP
      B(K2) = BJP
      A(K3) = AKM
      B(K3) = BKM
      KK = K3 + KSPAN
      IF (KK.LE.NT) GO TO 190
      GO TO 210
C  TRANSFORM FOR FACTOR OF 5 (OPTIONAL CODE)
  240 C2 = C72**2 - S72**2
      S2 = 2.0D0*C72*S72
  250 K1 = KK + KSPAN
      K2 = K1 + KSPAN
      K3 = K2 + KSPAN
      K4 = K3 + KSPAN
      AKP = A(K1) + A(K4)
      AKM = A(K1) - A(K4)
      BKP = B(K1) + B(K4)
      BKM = B(K1) - B(K4)
      AJP = A(K2) + A(K3)
      AJM = A(K2) - A(K3)
      BJP = B(K2) + B(K3)
      BJM = B(K2) - B(K3)
      AA = A(KK)
      BB = B(KK)
      A(KK) = AA + AKP + AJP
      B(KK) = BB + BKP + BJP
      AK = AKP*C72 + AJP*C2 + AA
      BK = BKP*C72 + BJP*C2 + BB
      AJ = AKM*S72 + AJM*S2
      BJ = BKM*S72 + BJM*S2
      A(K1) = AK - BJ
      A(K4) = AK + BJ
      B(K1) = BK + AJ
      B(K4) = BK - AJ
      AK = AKP*C2 + AJP*C72 + AA
      BK = BKP*C2 + BJP*C72 + BB
      AJ = AKM*S2 - AJM*S72
      BJ = BKM*S2 - BJM*S72
      A(K2) = AK - BJ
      A(K3) = AK + BJ
      B(K2) = BK + AJ
      B(K3) = BK - AJ
      KK = K4 + KSPAN
      IF (KK.LT.NN) GO TO 250
      KK = KK - NN
      IF (KK.LE.KSPAN) GO TO 250
      GO TO 320
C  TRANSFORM FOR ODD FACTORS
  260 K = NFAC(I)
      KSPNN = KSPAN
      KSPAN = KSPAN/K
      IF (K.EQ.3) GO TO 160
      IF (K.EQ.5) GO TO 240
      IF (K.EQ.JF) GO TO 280
      JF = K
      S1 = RAD/K
      C1 = COS(S1)
      S1 = SIN(S1)
      IF (JF.GT.MAXF) GO TO 590
      CK(JF) = 1.0D0
      SK(JF) = 0.0D0
      J = 1
  270 CK(J) = CK(K)*C1 + SK(K)*S1
      SK(J) = CK(K)*S1 - SK(K)*C1
      K = K - 1
      CK(K) = CK(J)
      SK(K) = -SK(J)
      J = J + 1
      IF (J.LT.K) GO TO 270
  280 K1 = KK
      K2 = KK + KSPNN
      AA = A(KK)
      BB = B(KK)
      AK = AA
      BK = BB
      J = 1
      K1 = K1 + KSPAN
  290 K2 = K2 - KSPAN
      J = J + 1
      AT(J) = A(K1) + A(K2)
      AK = AT(J) + AK
      BT(J) = B(K1) + B(K2)
      BK = BT(J) + BK
      J = J + 1
      AT(J) = A(K1) - A(K2)
      BT(J) = B(K1) - B(K2)
      K1 = K1 + KSPAN
      IF (K1.LT.K2) GO TO 290
      A(KK) = AK
      B(KK) = BK
      K1 = KK
      K2 = KK + KSPNN
      J = 1
  300 K1 = K1 + KSPAN
      K2 = K2 - KSPAN
      JJ = J
      AK = AA
      BK = BB
      AJ = 0.0D0
      BJ = 0.0D0
      K = 1
  310 K = K + 1
      AK = AT(K)*CK(JJ) + AK
      BK = BT(K)*CK(JJ) + BK
      K = K + 1
      AJ = AT(K)*SK(JJ) + AJ
      BJ = BT(K)*SK(JJ) + BJ
      JJ = JJ + J
      IF (JJ.GT.JF) JJ = JJ - JF
      IF (K.LT.JF) GO TO 310
      K = JF - J
      A(K1) = AK - BJ
      B(K1) = BK + AJ
      A(K2) = AK + BJ
      B(K2) = BK - AJ
      J = J + 1
      IF (J.LT.K) GO TO 300
      KK = KK + KSPNN
      IF (KK.LE.NN) GO TO 280
      KK = KK - NN
      IF (KK.LE.KSPAN) GO TO 280
C  MULTIPLY BY ROTATION FACTOR (EXCEPT FOR FACTORS OF 2 AND 4)
  320 IF (I.EQ.M) GO TO 360
      KK = JC + 1
  330 C2 = 1.0D0 - CD
      S1 = SD
  340 C1 = C2
      S2 = S1
      KK = KK + KSPAN
  350 AK = A(KK)
      A(KK) = C2*AK - S2*B(KK)
      B(KK) = S2*AK + C2*B(KK)
      KK = KK + KSPNN
      IF (KK.LE.NT) GO TO 350
      AK = S1*S2
      S2 = S1*C2 + C1*S2
      C2 = C1*C2 - AK
      KK = KK - NT + KSPAN
      IF (KK.LE.KSPNN) GO TO 350
      C2 = C1 - (CD*C1+SD*S1)
      S1 = S1 + (SD*C1-CD*S1)
      C1 = 0.5D0/(C2**2+S1**2) + 0.5D0
      S1 = C1*S1
      C2 = C1*C2
      KK = KK - KSPNN + JC
      IF (KK.LE.KSPAN) GO TO 340
      KK = KK - KSPAN + JC + INC
      IF (KK.LE.JC+JC) GO TO 330
      GO TO 120
C  PERMUTE THE RESULTS TO NORMAL ORDER--- DONE IN TWO STAGES
C  PERMUTATION FOR SQUARE FACTORS OF N
  360 NP(1) = KS
      IF (KT.EQ.0) GO TO 450
      K = KT + KT + 1
      IF (M.LT.K) K = K - 1
      J = 1
      NP(K+1) = JC
  370 NP(J+1) = NP(J)/NFAC(J)
      NP(K) = NP(K+1)*NFAC(J)
      J = J + 1
      K = K - 1
      IF (J.LT.K) GO TO 370
      K3 = NP(K+1)
      KSPAN = NP(2)
      KK = JC + 1
      K2 = KSPAN + 1
      J = 1
      IF (N.NE.NTOT) GO TO 410
C  PERMUTATION FOR SINGLE-VARIATE TRANSFORM (OPTIONAL CODE)
  380 AK = A(KK)
      A(KK) = A(K2)
      A(K2) = AK
      BK = B(KK)
      B(KK) = B(K2)
      B(K2) = BK
      KK = KK + INC
      K2 = KSPAN + K2
      IF (K2.LT.KS) GO TO 380
  390 K2 = K2 - NP(J)
      J = J + 1
      K2 = NP(J+1) + K2
      IF (K2.GT.NP(J)) GO TO 390
      J = 1
  400 IF (KK.LT.K2) GO TO 380
      KK = KK + INC
      K2 = KSPAN + K2
      IF (K2.LT.KS) GO TO 400
      IF (KK.LT.KS) GO TO 390
      JC = K3
      GO TO 450
C  PERMUTATION FOR MULTIVARIATE TRANSFORM
  410 K = KK + JC
  420 AK = A(KK)
      A(KK) = A(K2)
      A(K2) = AK
      BK = B(KK)
      B(KK) = B(K2)
      B(K2) = BK
      KK = KK + INC
      K2 = K2 + INC
      IF (KK.LT.K) GO TO 420
      KK = KK + KS - JC
      K2 = K2 + KS - JC
      IF (KK.LT.NT) GO TO 410
      K2 = K2 - NT + KSPAN
      KK = KK - NT + JC
      IF (K2.LT.KS) GO TO 410
  430 K2 = K2 - NP(J)
      J = J + 1
      K2 = NP(J+1) + K2
      IF (K2.GT.NP(J)) GO TO 430
      J = 1
  440 IF (KK.LT.K2) GO TO 410
      KK = KK + JC
      K2 = KSPAN + K2
      IF (K2.LT.KS) GO TO 440
      IF (KK.LT.KS) GO TO 430
      JC = K3
  450 IF (2*KT+1.GE.M) RETURN
      KSPNN = NP(KT+1)
C  PERMUTATION FOR SQUARE-FREE FACTORS OF N
      J = M - KT
      NFAC(J+1) = 1
  460 NFAC(J) = NFAC(J)*NFAC(J+1)
      J = J - 1
      IF (J.NE.KT) GO TO 460
      KT = KT + 1
      NN = NFAC(KT) - 1
      IF (NN.GT.MAXP) GO TO 590
      JJ = 0
      J = 0
      GO TO 490
  470 JJ = JJ - K2
      K2 = KK
      K = K + 1
      KK = NFAC(K)
  480 JJ = KK + JJ
      IF (JJ.GE.K2) GO TO 470
      NP(J) = JJ
  490 K2 = NFAC(KT)
      K = KT + 1
      KK = NFAC(K)
      J = J + 1
      IF (J.LE.NN) GO TO 480
C  DETERMINE THE PERMUTATION CYCLES OF LENGTH GREATER THAN 1
      J = 0
      GO TO 510
  500 K = KK
      KK = NP(K)
      NP(K) = -KK
      IF (KK.NE.J) GO TO 500
      K3 = KK
  510 J = J + 1
      KK = NP(J)
      IF (KK.LT.0) GO TO 510
      IF (KK.NE.J) GO TO 500
      NP(J) = -J
      IF (J.NE.NN) GO TO 510
      MAXF = INC*MAXF
C  REORDER A AND B, FOLLOWING THE PERMUTATION CYCLES
      GO TO 580
  520 J = J - 1
      IF (NP(J).LT.0) GO TO 520
      JJ = JC
  530 KSPAN = JJ
      IF (JJ.GT.MAXF) KSPAN = MAXF
      JJ = JJ - KSPAN
      K = NP(J)
      KK = JC*K + II + JJ
      K1 = KK + KSPAN
      K2 = 0
  540 K2 = K2 + 1
      AT(K2) = A(K1)
      BT(K2) = B(K1)
      K1 = K1 - INC
      IF (K1.NE.KK) GO TO 540
  550 K1 = KK + KSPAN
      K2 = K1 - JC*(K+NP(K))
      K = -NP(K)
  560 A(K1) = A(K2)
      B(K1) = B(K2)
      K1 = K1 - INC
      K2 = K2 - INC
      IF (K1.NE.KK) GO TO 560
      KK = K2
      IF (K.NE.J) GO TO 550
      K1 = KK + KSPAN
      K2 = 0
  570 K2 = K2 + 1
      A(K1) = AT(K2)
      B(K1) = BT(K2)
      K1 = K1 - INC
      IF (K1.NE.KK) GO TO 570
      IF (JJ.NE.0) GO TO 530
      IF (J.NE.1) GO TO 520
  580 J = K3 + 1
      NT = NT - KSPNN
      II = NT - INC + 1
      IF (NT.GE.0) GO TO 520
      RETURN
C  ERROR FINISH, INSUFFICIENT ARRAY STORAGE
  590 ISN = 0
      CALL IPRINT(IPRT)
      WRITE(IPRT, 1000)
C
C     NB.  THE FOLLOWING STOP SHOULD BE CHANGED TO A RETURN WHEN
C          THE TIME SERIES ROUTINES ARE MODIFIED FOR STARPAC.
C
      STOP
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (' ', 17('*')/18H * ERROR MESSAGE */1X, 17('*')//
     +   45H ARRAY BOUNDS EXCEEDED WITHIN SUBROUTINE FFT./
     +   44H PLEASE BRING THIS ERROR TO THE ATTENTION OF/
     +   22H    JANET R. DONALDSON/
     +   16H    303-497-5114/
     +   16H    FTS 320-5114)
      END
