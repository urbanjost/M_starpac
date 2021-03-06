*DFBW
      SUBROUTINE DFBW (N, LAG, W, LW, DF, BW)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS SUBROUTINE COMPUTES AND STORES THE ASSOCIATED DEGREES OF
C     FREEDOM AND BANDWIDTH FOR A GIVEN LAG WINDOW.
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
      REAL
     +   BW,DF
      INTEGER
     +   LAG,LW,N
C
C  ARRAY ARGUMENTS
      REAL
     +   W(LW)
C
C  LOCAL SCALARS
      INTEGER
     +   K
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     REAL BW
C        THE BANDWIDTH.
C     REAL DF
C        THE EFFECTIVE DEGREES OF FREEDOM.
C     INTEGER K
C        AN INDEX VARIABLE
C     INTEGER LAG
C        THE LAG WINDOW TRUNCATION POINT USED FOR A SPECIFIC WINDOW.
C     INTEGER LW
C        THE LENGTH OF THE VECTOR W.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS IN THE SERIES.
C     REAL W(LW)
C        THE VECTOR OF LAG WINDOWS.
C
      BW = 0.0E0
      DO 10 K = 1, LAG
         BW = BW + W(K+1) * W(K+1) * (N-K)
   10 CONTINUE
C
      BW = 1.0E0 / (W(1)*W(1) + 2.0E0*BW/N)
      DF = 2.0E0 * BW * N
      RETURN
      END
