*UFSLAG
      SUBROUTINE UFSLAG (ACOV, LAGMAX, LAGS, N, NW, NWUSED, LACOV)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE COMPUTES THE LAG WINDOW TRUNCATION POINTS FOR
C     SPECTRUM ANALYSIS.
C
C     WRITTEN BY - JANET R. DONALDSON
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  DECEMBER 7, 1981
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   LACOV,LAGMAX,N,NW,NWUSED
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   ACOV(LACOV)
      INTEGER
     +   LAGS(NW)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   ACOVMX,P95LIM
      INTEGER
     +   I,J,K,LAG,NWM1
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,DBLE,SQRT
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION ACOV(LACOV)
C        THE ARRAY IN WHICH THE AUTOCOVARIANCES ARE STORED
C     DOUBLE PRECISION ACOVMX
C        THE MAXIMUM AUTOCOVARIANCE VALUE.
C     INTEGER I
C        AN INDEX VARIABLE
C     INTEGER J, K
C        INDEX VARIABLES.
C     INTEGER LACOV
C        THE LENGTH OF VECTOR ACOV.
C     INTEGER LAG, LAGMAX
C        THE INDEXING VARIABLE INDICATING THE LAG VALUE OF THE
C        AUTOCOVARIANCE BEING COMPUTED AND THE MAXIMUM LAG TO BE USED,
C        RESPECTIVELY.
C     INTEGER LAGS(NW)
C        THE ARRAY USED TO STORE THE LAG WINDOW TRUCCATION
C        POINTS USED FOR EACH SET OF SPECTRUM VALUES.
C     INTEGER N
C        THE INTEGER NUMBER OF OBSERVATIONS IN EACH SERIES
C     INTEGER NW
C        THE NUMBER OF DIFFERENT BANDWIDTHS REQUESTED.
C     INTEGER NWM1, NWUSED
C        THE NUMBER OF DIFFERENT BANDWIDTHS MINUS 1, AND THE
C        ACTUAL NUMBER OF BANDWIDTHS ACTUALLY USED.
C     DOUBLE PRECISION P95LIM
C        THE 95 PERCENT CONFIDENT LIMIT FOR WHITE NOISE.
C
      LAGS(NW) = LAGMAX
      IF (LAGS(NW) .LE. 32) GO TO 30
C
C     COMPUTE 95 PERCENT CONFIDENCE LIMITS ON AUTOCOVARIANCES,
C     ASSUMING WHITE NOISE.
C
      P95LIM = 1.96D0 * ACOV(1) / SQRT(DBLE(N))
C
C     CHECK FOR FIRST ACVF EXCEEDING 95 PERCENT LIMIT ON WHITE NOISE
C
      DO 10 I = 1, LAGMAX
         LAG = LAGMAX + 1 - I
         IF (ABS(ACOV(LAG + 1)) .GE. P95LIM) GO TO 30
         LAGS(NW) = LAGS(NW) - 1
   10 CONTINUE
C
C     IF NO ACVF EXCEEDS WHITE NOISE LIMITS, CHECK FOR LARGEST ACVF.
C
      LAGS(NW) = 1
      ACOVMX = ABS(ACOV(2))
      DO 20 LAG = 1, LAGMAX
         IF (ABS(ACOV(LAG + 1)) .LE. ACOVMX) GO TO 20
         LAGS(NW) = LAG
         ACOVMX = ABS(ACOV(LAG + 1))
   20 CONTINUE
C
C     COMPUTE LAG WINDOW TRUNCATION POINTS
C
   30 LAGS(NW) = LAGS(NW) * 3.0D0 / 2.0D0
      IF (LAGS(NW) .LT. 32) LAGS(NW) = 32
      IF (LAGS(NW) .GT. LAGMAX) LAGS(NW) = LAGMAX
      NWUSED = NW
      IF (NW .EQ. 1) RETURN
      NWM1 = NW - 1
      DO 40 I = 1, NWM1
         K = NW - I
         LAGS(K) = LAGS(K + 1) / 2
   40 CONTINUE
C
C     CHECK WHETHER ALL NW LAG WINDOW TRUNCATION POINTS CAN BE USED.
C
      NWUSED = NW
      IF (LAGS(1) .GE. 4) RETURN
C
C     RECONSTURCT -LAGS- VECTOR IF NOT ALL TRUNCATION POINTS ARE
C     TO BE USED
C
      DO 50 I = 2, NW
         NWUSED = NWUSED - 1
         IF (LAGS(I) .GE. 4) GO TO 60
   50 CONTINUE
C
   60 DO 70 I = 1, NWUSED
         J = NW - NWUSED + I
         LAGS(I) = LAGS(J)
   70 CONTINUE
C
      RETURN
      END
