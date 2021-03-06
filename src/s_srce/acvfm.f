*ACVFM
      SUBROUTINE ACVFM (Y, YMISS, N, YMEAN, ACOV, LAGMAX,
     +   LAGLST, NLPPA, LACOV)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE COMPUTES THE AUTOCOVARIANCES WHEN MISSING DATA ARE
C     INVOLVED.
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
     +   YMEAN,YMISS
      INTEGER
     +   LACOV,LAGLST,LAGMAX,N
C
C  ARRAY ARGUMENTS
      REAL
     +   ACOV(*),Y(*)
      INTEGER
     +   NLPPA(*)
C
C  LOCAL SCALARS
      REAL
     +   DOTXY,DOTYY,FPLM
      INTEGER
     +   LAG,NDOTXY,NDOTYY,NUSED
C
C  EXTERNAL FUNCTIONS
      REAL
     +   R1MACH
      INTEGER
     +   LSTLAG
      EXTERNAL R1MACH,LSTLAG
C
C  EXTERNAL SUBROUTINES
      EXTERNAL AMEANM,DOTCM
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     REAL ACOV(LACOV)
C        THE ARRAY IN WHICH THE AUTOCOVARIANCES ARE STORED
C     REAL DOTXY, DOTYY
C        THE DOT PRODUCT BETWEEN VECTORS (Y(I) - YMEAN)) AND
C        (Y(LAG) - YMEAN)), AND (Y(I) - YMEAN)) AND (Y(I) - YMEAN)),
C        RESPECTIVELY.
C     REAL FPLM
C        THE FLOATING POINT LARGEST MAGNITUDE.
C     INTEGER LACOV
C        THE LENGTH OF THE VECTOR ACOV.
C     INTEGER LAG, LAGLST, LAGMAX
C        THE INDEXING VARIABLE INDICATING THE LAG VALUE OF THE
C        AUTOCORRELATION BEING COMPUTED, THE NUMBER OF AUTOCORRELATIONS
C        COMPUTED BEFORE A MISSING AUTOCORRELATION, AND THE NUMBER OF
C        AUTOCORRELATIONS DESIRED, RESPECTIVELY.
C     INTEGER N
C        THE INTEGER NUMBER OF OBSERVATIONS IN THE SERIES
C     INTEGER NDOTXY, NDOTYY
C        THE NUMBER OF OBSERVATIONS USED TO COMPUTE DOTXY AND
C        DOTYY, RESPECTIVELY.
C     INTEGER NLPPA(LACOV)
C        THE ARRAY CONTAINING THE NUMBERS OF LAGGED PRODUCT PAIRS
C        USED TO COMPUTE THE ACVF AT EACH LAG.
C     INTEGER NUSED
C        THE NUMBER OF ACTIVE OBSERVATIONS IN THE SERIES.
C     REAL Y(N)
C        THE VECTOR CONTAINING THE OBSERVED SERIES
C     REAL YMEAN
C        THE MEAN OF THE OBSERVED TIME SERIES
C     REAL YMISS
C        THE USER SUPPLIED CODE WHICH IS USED TO DETERMINE WHETHER OR
C        NOT AN OBSERVATION IN THE SERIES IS MISSING.  IF Y(I) = YMISS,
C        THE VALUE IS ASSUMED MISSING, OTHERWISE IT IS NOT.
C
C
      FPLM = R1MACH(2)
C
C     COMPUTE ARITHMETIC MEAN, WITH MISSING VALUES TAKEN INTO ACCOUNT
C
      CALL AMEANM (Y, YMISS, N, NUSED, YMEAN)
C
C     COMPUTE THE VARIANCE OF THE SERIES Y
C
      CALL DOTCM (Y, YMEAN, YMISS, N, Y, YMEAN, YMISS, N,
     +   DOTYY, NDOTYY)
      NLPPA(1) = NDOTYY
      IF (NLPPA(1).EQ.0) THEN
         LAGLST = 0
      ELSE
         ACOV(1) = DOTYY / NDOTYY
C
C     COMPUTE AUTOCORRELATIONS, WITH MISSING VALUES TAKEN INTO ACCOUNT
C
         DO 10 LAG = 1, LAGMAX
            CALL DOTCM (Y, YMEAN, YMISS, N, Y(LAG+1), YMEAN,
     +         YMISS, N - LAG, DOTXY, NDOTXY)
            NLPPA(LAG + 1) = NDOTXY
            ACOV(LAG + 1) = FPLM
            IF (NLPPA(LAG + 1) .LE. 0) GO TO 10
            ACOV(LAG + 1) = DOTXY * (N-LAG) / (NLPPA(LAG + 1) * N)
   10    CONTINUE
C
C     FIND THE LAST AUTOCORRELATION TO BE COMPUTED BEFORE
C     ONE COULD NOT BE COMPUTED DUE TO MISSING DATA
C
         LAGLST = LSTLAG(NLPPA, LAGMAX, LACOV)
      END IF
      RETURN
      END
