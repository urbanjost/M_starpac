*CCVFM
      SUBROUTINE CCVFM(Y1, Y1MISS, Y2, Y2MISS, N, NC, Y1MEAN, Y2MEAN,
     +   CCOV12, CCOV21, ICCOV, NLPP12, NLPP21)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE COMPUTES THE CROSS COVARIANCE FUNCTION
C     BETWEEN TWO SERIES WHEN MISSING DATA ARE INVOLVED.
C
C     WRITTEN BY - JANET R. DONALDSON
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  NOVEMBER 21, 1980
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      DOUBLE PRECISION
     +   Y1MEAN,Y1MISS,Y2MEAN,Y2MISS
      INTEGER
     +   ICCOV,N,NC
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   CCOV12(ICCOV),CCOV21(ICCOV),Y1(N),Y2(N)
      INTEGER
     +   NLPP12(ICCOV),NLPP21(ICCOV)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   DOTXY,FPLM
      INTEGER
     +   LAG,NDOTXY
C
C  EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   D1MACH
      EXTERNAL D1MACH
C
C  EXTERNAL SUBROUTINES
      EXTERNAL DOTCM
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION CCOV12(ICCOV), CCOV21(ICCOV)
C        THE ARRAYS IN WHICH THE CCVF FOR SERIES 1 LAGGED
C        BEHIND SERIES 2 AND VISA VERSA, RESPECTIVELY, ARE
C        STORED.
C     DOUBLE PRECISION DOTXY
C        VARIOUS CROSS PRUDUCTS BETWEEN SERIES Y1 AND Y2.
C     DOUBLE PRECISION FPLM
C        THE FLOATING POINT LARGEST MAGNITUDE.
C     INTEGER ICCOV
C        THE ROW DIMENSION OF THE ARRAYS CCOV12 AND CCOV21.
C     INTEGER LAG
C        THE INDEXING VARIABLE INDICATING THE LAG VALUE OF THE
C        AUTOCORRELATION BEING COMPUTED.
C     INTEGER N
C        THE INTEGER NUMBER OF OBSERVATIONS IN THE SERIES
C     INTEGER NC
C        THE NUMBER OF CROSS CORRELATIONS DESIRED.
C     INTEGER NLPP12(ICCOV), NLPP21(ICCOV)
C        THE NUMBER OF LAGGED PRODUCT PAIRS USED TO COMPUTE THE CCVF
C        FOR EACH PAIR OF SERIES AT EACH LAG.
C     INTEGER NDOTXY
C        THE NUMBER OF OBSERVATIONS USED TO COMPUTE DOTXY.
C     DOUBLE PRECISION Y1(N), Y1MEAN, Y1MISS
C        THE FIRST SERIES, AND ITS MEAN, AND MISSING VALUE CODE.
C     DOUBLE PRECISION Y2(N), Y2MEAN, Y2MISS
C        THE SECOND SERIES, AND ITS MEAN, AND MISSING VALUE CODE.
C
      FPLM = D1MACH(2)
C
C     COMPUTE THE CROSS COVARIANCES
C
      CALL DOTCM (Y1, Y1MEAN, Y1MISS, N, Y2, Y2MEAN, Y2MISS, N, DOTXY,
     +   NDOTXY)
C
      NLPP12(1) = NDOTXY
      CCOV12(1) = FPLM
      IF (NDOTXY .GE. 1) CCOV12(1) = DOTXY / NDOTXY
C
      CCOV21(1) = CCOV12(1)
      NLPP21(1) = NDOTXY
C
      DO 10 LAG = 1, NC
C
         CALL DOTCM (Y1, Y1MEAN, Y1MISS, N, Y2(LAG+1), Y2MEAN, Y2MISS,
     +       N-LAG, DOTXY, NDOTXY)
C
         NLPP12(LAG+1) = NDOTXY
         CCOV12(LAG+1) = FPLM
         IF (NDOTXY .GE. 1)
     +      CCOV12(LAG+1) = DOTXY * (N-LAG) / (N*NDOTXY)
C
         CALL DOTCM (Y2, Y2MEAN, Y2MISS, N, Y1(LAG+1), Y1MEAN, Y1MISS,
     +      N-LAG, DOTXY, NDOTXY)
C
         NLPP21(LAG+1) = NDOTXY
         CCOV21(LAG+1) = FPLM
         IF (NDOTXY .GE. 1)
     +      CCOV21(LAG+1) = DOTXY * (N-LAG) / (N*NDOTXY)
C
   10 CONTINUE
C
      RETURN
      END
