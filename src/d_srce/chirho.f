*CHIRHO
      SUBROUTINE CHIRHO (RHO, N, NC, CHI, CHIP)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE COMPUTES THE CHI SQUARED STATISTIC AND ITS
C     PROBABILITY BASED IN A VECTOR OF AUTOCORRELATIONS.
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
     +   CHI,CHIP
      INTEGER
     +   N,NC
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   RHO(*)
C
C  LOCAL SCALARS
      INTEGER
     +   LAG
C
C  EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   CDFCHI
      EXTERNAL CDFCHI
C
C  INTRINSIC FUNCTIONS
      INTRINSIC DBLE
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION CHI, CHIP
C        THE VARIABLES IN WHICH THE CHI SQUARE STATISTIC AND
C        CHI SQUARE STATISTIC PROBABILITY ARE STORED.
C     INTEGER LAG
C        THE INDEXING VARIABLE INDICATING THE LAG VALUE OF THE
C        AUTOCORRELATION BEING EXAMINED.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS IN THE SERIES.
C     INTEGER NC
C        THE NUMBER OF AUTOCORRELATIONS COMPUTED.
C     DOUBLE PRECISION RHO(NC)
C        THE ARRAY IN WHICH THE AUTOCORRELATIONS ARE STORED
C
      CHI = 0.0D0
      DO 10 LAG = 1, NC
         CHI = CHI + RHO(LAG) * RHO(LAG)
   10 CONTINUE
      CHI = CHI * N
      CHIP = 1.0D0 - CDFCHI(CHI, DBLE(NC))
      RETURN
      END
