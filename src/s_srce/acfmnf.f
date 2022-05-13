*ACFMNF
      SUBROUTINE ACFMNF (YFFT, N, NFFT, LAGMAX, RHO, SDRHO, YMEAN,
     +   PRHO, AIC, FTEST, PHI, IAR, OSPVAR, ACOV, LACOV, LAIC,
     +   CHIA, CHIAP, LYFFT, WORK, LWORK, NPRT)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE MAIN SUBROUTINE FOR COMPUTING AUTOCORRELATIONS AND
C     PARTIAL AUTOCORRELATIONS OF A TIME SERIES .
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
     +   CHIA,CHIAP,OSPVAR,YMEAN
      INTEGER
     +   IAR,LACOV,LAGMAX,LAIC,LWORK,LYFFT,N,NFFT,NPRT
C
C  ARRAY ARGUMENTS
      REAL
     +   ACOV(*),AIC(*),FTEST(2,*),PHI(*),PRHO(*),RHO(*),SDRHO(*),
     +   WORK(*),YFFT(*)
C
C  LOCAL SCALARS
      INTEGER
     +   I
C
C  EXTERNAL SUBROUTINES
      EXTERNAL ACFSD,ACVFF,AOS,CHIRHO
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     REAL ACOV(LACOV)
C        THE AUTOCOVARIANCE FUNCTION ESTIMATE VECTOR.
C     REAL AIC(LAIC)
C        THE AKAIKES INFORMATION CRITERION VECTOR.
C     REAL CHIA, CHIAP
C        THE VARIABLES IN WHICH THE CHI SQUARE STATISTIC AND
C        CHI SQUARED STATISTIC PROBABILITY FOR THE AUTOCORRELATIONS
C        ARE STORED.
C     REAL FTEST(2, LAGMAX)
C        THE ARRAY IN WHICH THE PARTIAL F RATIOS AND PROBABILITIES
C        ARE STORED.
C     INTEGER I
C        AN INDEX VARIABLE.
C     INTEGER IAR
C        THE ORDER OF THE AUTOREGRESSIVE PROCESS CHOSEN.
C     INTEGER LACOV
C        THE LENGTH OF THE VECTOR ACOV.
C     INTEGER LAGMAX
C        THE MAXIMUM LAG VALUE TO BE USED.
C     INTEGER LAIC
C        THE LENGTH OF THE VECTOR AIC.
C     INTEGER LWORK
C        THE LENGTH OF THE VECTOR WORK.
C     INTEGER LYFFT
C        THE LENGTH OF THE VECTOR YFFT.
C     INTEGER N
C        THE INTEGER NUMBER OF OBSERVATIONS IN THE SERIES
C     INTEGER NFFT
C        THE NUMBER OF OBSERVATIONS IN THE EXTENDED SERIES.
C     INTEGER NPRT
C        THE INDICATOR VARIABLE USED TO SPECIFY WHETHER OR NOT
C        PRINTED OUTPUT IS TO BE GIVEN, WHERE IF THE VALUE OF
C        NPRT IS ZERO, NO OUTPUT IS MADE.
C     REAL OSPVAR
C        THE ONE STEP PREDICTION VARIANCE.
C     REAL PHI(LAGMAX)
C        THE ARRAY OF AUTOREGRESSIVE COEFFICIENTS FOR THE SELECTED
C        ORDER (IAR).
C     REAL PRHO(LAGMAX)
C        THE ARRAY IN WHICH THE PARTIAL AUTOCORRELATIONS ARE STORED
C     REAL RHO(LAGMAX)
C        THE ARRAY IN WHICH THE AUTOCORRELATIONS ARE STORED
C     REAL SDRHO(LAGMAX)
C        THE ARRAY IN WHICH THE STANDARD ERRORS OF THE AUTOCORRELATIONS
C        ARE STORED
C     REAL WORK(LWORK)
C        A WORK ARRAY.
C     REAL YFFT(LYFFT)
C        THE VECTOR CONTAINING THE OBSERVED SERIES
C     REAL YMEAN
C        THE MEAN OF THE OBSERVED TIME SERIES
C
C     COMPUTE AUTOCOVARIANCESS AND STANDARD DEVIATION OF THE SERIES.
C
      CALL ACVFF(YFFT, N, NFFT, YMEAN, ACOV, LAGMAX, LACOV, LYFFT, WORK,
     +   LWORK)
C
      IF (ACOV(1) .EQ. 0.0E0) RETURN
C
C     COMPUTE PARTIAL AUTOCORRELATIONS AND AUTOREGRESSIVE ORDER
C     SELECTION STATISTICS.
C
      CALL AOS (N, LAGMAX, ACOV, PRHO, IAR, OSPVAR, PHI, WORK,
     +   AIC, FTEST, LACOV, LAIC)
C
      IF (NPRT .EQ. 0) RETURN
C
C     COMPUTE AUTOCORRELATIONS
C
      DO 10 I = 1, LAGMAX
         RHO(I) = ACOV(I+1) / ACOV(1)
   10 CONTINUE
C
C     COMPUTE STANDARD ERROR OF AUTOCORRELATIONS.
C
      CALL ACFSD (RHO, SDRHO, LAGMAX, N)
C
C     COMPUTE CHI STATISTIC BASED ON AUTOCORRELATION VALUES
C
      CALL CHIRHO (RHO, N, LAGMAX, CHIA, CHIAP)
C
      RETURN
      END
