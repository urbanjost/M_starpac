*CCFMNM
      SUBROUTINE CCFMNM (Y1, Y1MISS, Y2, Y2MISS, N, LAGMAX, NCC,
     +   CCOV11, CCOV22, CCOV12, CCOV21, ICCOV, NLPP11, NLPP22,
     +   NLPP12, NLPP21, INLPPC, Y1MEAN, Y2MEAN, RHOC, SDRHOC, NPRT,
     +   LAGLST)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE MAIN SUBROUTINE FOR COMPUTING CROSS CORRELATIONS AND
C     THEIR STANDARD ERRORS WHEN MISSING DATA ARE INVOLVED.
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
      REAL
     +   Y1MEAN,Y1MISS,Y2MEAN,Y2MISS
      INTEGER
     +   ICCOV,INLPPC,LAGLST,LAGMAX,N,NCC,NPRT
C
C  ARRAY ARGUMENTS
      REAL
     +   CCOV11(ICCOV),CCOV12(ICCOV),CCOV21(ICCOV),CCOV22(ICCOV),
     +   RHOC(NCC),SDRHOC(NCC),Y1(N),Y2(N)
      INTEGER
     +   NLPP11(INLPPC),NLPP12(INLPPC),NLPP21(INLPPC),
     +   NLPP22(INLPPC)
C
C  LOCAL SCALARS
      REAL
     +   FAC,FPLM
      INTEGER
     +   I,I0,IM,IP
C
C  EXTERNAL FUNCTIONS
      REAL
     +   R1MACH
      EXTERNAL R1MACH
C
C  EXTERNAL SUBROUTINES
      EXTERNAL CCFSDM,CCVFM
C
C  INTRINSIC FUNCTIONS
      INTRINSIC SQRT
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     REAL CCOV11(ICCOV), CCOV12(ICCOV)
C     REAL CCOV21(ICCOV), CCOV22(ICCOV)
C        THE ARRAY CONTAINING THE AUTOCOVARIANCE AND CROSS COVARIANCE
C        ESTIMATES FOR SERIES 1 AND 2.
C     REAL FAC
C        THE INVERSE OF THE SQUARE ROOT OF THE PRODUCT OF THE
C        AUTOCOVARIANCES AT LAG ZERO.
C     REAL FPLM
C        THE FLOATING POINT LARGEST MAGNITUDE.
C     INTEGER I
C        THE INDEXING VARIABLE FOR THE LAG.
C     INTEGER ICCOV
C        THE DIMENSION OF THE COVARIANCE VECTORS.
C     INTEGER IM
C        THE LOCATIONS IN THE VARIOUS CCF RELATED ARRAYS OF LAG -I.
C     INTEGER INLPPC
C        THE DIMENSION OF THE LAGGED PRODUCT PAIR COUNT VECTORS.
C     INTEGER IP
C        THE LOCATION IF THE VARIOUS CCF RELATED ARRAYS OF LAG I.
C     INTEGER I0
C        THE LOCATION IF THE VARIOUS CCF RELATED ARRAYS OF LAG 0.
C     INTEGER LAGLST
C        THE LAST LAG BEFORE MISSING DATA CAUSED THE ACVF OF EITHER
C        SERIES 1 OR 2 NOT TO BE COMPUTED.
C     INTEGER LAGMAX
C        THE MAXIMUM LAG VALUE TO BE USED.
C     INTEGER N
C        THE INTEGER NUMBER OF OBSERVATIONS IN THE SERIES
C     INTEGER NCC
C        THE NUMBER OF CROSS CORRELATIONS TO BE COMPUTED.
C     INTEGER NLPP11(INLPPC), NLPP12(INLPPC), NLPP21(INLPPC),
C    +        NLPP22(INLPPC)
C        THE NUMBERS OF LAGGED PRODUCT PAIRS USED TO COMPUTE
C        THE AUTOCOVARIANCE AND CROSS COVARIANCE ESTIMATES.
C     INTEGER NPRT
C        THE VARIABLE USED TO CONTROL PRINTED OUTPUT.
C     REAL RHOC(NCC)
C        THE ARRAY IN WHICH THE AUTO AND CROSS CORRELATIONS ARE STORED
C     REAL SDRHOC(NCC)
C        THE ARRAY CONTAINING THE STD. ERRORS OF THE CROSS CORRELATIONS.
C        ARE STORED
C     REAL Y1(N), Y1MEAN, Y1MISS
C        THE FIRST SERIES, AND ITS MEAN, AND MISSING VALUE CODE.
C     REAL Y2(N), Y2MEAN, Y2MISS
C        THE SECOND SERIES, AND ITS MEAN, AND MISSING VALUE CODE.
C
      FPLM = R1MACH(2)
C
C     COMPUTE AUTOCORRELATIONS AND STANDARD DEVIATION OF THE SERIES.
C
      CALL CCVFM(Y1, Y1MISS, Y2, Y2MISS, N, LAGMAX, Y1MEAN, Y2MEAN,
     +   CCOV12, CCOV21, ICCOV, NLPP12, NLPP21)
C
      IF (NPRT .EQ. 0 .OR. NLPP11(1) .EQ. 0) RETURN
      IF (CCOV11(1) *CCOV22(1) .EQ. 0.0E0) RETURN
C
      FAC = 1.0E0 / SQRT(CCOV11(1) * CCOV22(1))
C
      I0 = LAGMAX + 1
      RHOC(I0) = FPLM
      IF (NLPP12(1).GE.1) RHOC(I0) = CCOV12(1) * FAC
C
      DO 10 I = 1, LAGMAX
         IP = I0 + I
         RHOC(IP) = FPLM
         IF (NLPP12(I+1).GE.1) RHOC(IP) = CCOV12(I+1) * FAC
C
         IM = I0 - I
         RHOC(IM) = FPLM
         IF (NLPP21(I+1).GE.1) RHOC(IM) = CCOV21(I+1) * FAC
   10 CONTINUE
C
C     COMPUTE STANDARD ERROR OF AUTOCORRELATIONS.
C
      CALL CCFSDM (CCOV11, CCOV22, SDRHOC, LAGMAX, NCC, LAGLST, N,
     +   NLPP12, NLPP21, ICCOV, INLPPC)
C
      RETURN
      END
