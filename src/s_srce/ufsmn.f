*UFSMN
      SUBROUTINE UFSMN(ACOV, NLPPA, LAG, DF, NF, FREQ, ALPHA, BW, SPCF,
     +   ALOW, AUP, LACOV, ISPCF, WINDOW, W, LW, N, DELTA, MISS, LNLPPA)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE MAIN SUBROUTINE FOR COMPUTING AUTOCORRELATIONS AND
C     PARTIAL AUTOCORRELATIONS OF A TIME SERIES
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
     +   ALOW,ALPHA,AUP,BW,DELTA,DF
      INTEGER
     +   ISPCF,LACOV,LAG,LNLPPA,LW,N,NF
      LOGICAL
     +   MISS
C
C  ARRAY ARGUMENTS
      REAL
     +   ACOV(LACOV),FREQ(NF),SPCF(ISPCF),W(LW)
      INTEGER
     +   NLPPA(LNLPPA)
C
C  SUBROUTINE ARGUMENTS
      EXTERNAL WINDOW
C
C  EXTERNAL FUNCTIONS
      REAL
     +   PPFCHS
      EXTERNAL PPFCHS
C
C  EXTERNAL SUBROUTINES
      EXTERNAL DFBW,DFBWM,UFSEST
C
C  INTRINSIC FUNCTIONS
      INTRINSIC NINT
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     REAL ACOV(LACOV)
C        THE AUTOCOVARIANCES OF THE SERIES.
C     REAL ALOW
C        A FACTOR USED TO COMPUTE THE LOWER CONFIDENCE LIMITS.
C     REAL ALPHA
C        THE DESIRED CONFIDENCE LEVEL.
C     REAL AUP
C        A FACTOR USED TO COMPUTE THE UPPER CONFIDENCE LIMITS.
C     REAL BW
C        THE BANDWIDTH.
C     REAL DELTA
C        THE SAMPLING INTERVAL.
C     REAL DF
C        THE EFFECTIVE DEGREES OF FREEDOM.
C     REAL FREQ(NF)
C        THE VECTOR OF FREQUENCIES AT WHICH THE SPECTRUM IS TO BE
C        COMPUTED.
C     INTEGER ISPCF
C         THE ACTUAL FIRST DIMENSION OF THE SPECTRUM ARRAYS.
C     INTEGER LACOV
C        THE LENGTH OF VECTOR ACOV.
C     INTEGER LAG
C        THE VARIABLE INDICATING THE LAG VALUE BEING EXAMINED.
C     INTEGER LNLPPA
C        THE LENGTH OF THE VECTOR NLPPA.
C     INTEGER LW
C        THE LENGTH OF THE VECTOR W.
C     LOGICAL MISS
C        AN INDICATOR VARIABLE WHICH DESIGNATES WHETHER THERE ARE
C        MISSING VALUES (TRUE) OR NOT (FALSE)
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS IN THE TIME SERIES.
C     INTEGER NF
C        THE NUMBER OF FREQUENCIES AT WHICH THE SPECTRUM IS
C        TO BE COMPUTED.
C     INTEGER NLPPA(LNLPPA)
C        THE NUMBERS OF LAGGED PRODUCT PAIRS IN EACH ACVF VALUE.
C     REAL SPCF(ISPCF)
C        THE ARRAY IN WHICH THE SPECTRUM IS STORED.
C     REAL W(LW)
C        THE VECTOR OF LAG WINDOWS.
C     EXTERNAL WINDOW
C        THE NAME OF THE WINDOW COMPUTING SUBROUTINE.
C
C     COMPUTE THE WINDOW, EFFECTIVE DEGREES OF FREEDOM AND
C     BANDWIDTH BASED ON THE WINDOW.
C
      CALL WINDOW(LAG, W, LW)
      IF (.NOT.MISS) CALL DFBW(N, LAG, W, LW, DF, BW)
      IF (MISS) CALL DFBWM(N, LAG, W, LW, NLPPA, NLPPA, LNLPPA, DF, BW)
C
C     COMPUTE THE SPECTRUM
C
      CALL UFSEST(ACOV, W, LAG, SPCF, ISPCF, LACOV, LW, NF, FREQ, DELTA)
C
C     COMPUTE -ALPHA- PERCENT POINT FUNCTION VALUE FOR
C     SPECTRUM WINDOW BEING USED.
C
      ALOW = DF/PPFCHS(0.5E0+ALPHA/2.0E0,NINT(DF))
      AUP = DF/PPFCHS(0.5E0-ALPHA/2.0E0,NINT(DF))
C
      RETURN
      END
