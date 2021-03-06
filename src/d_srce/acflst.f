*ACFLST
      SUBROUTINE ACFLST (RHO, SDRHO, NLPPA, LAGMAX, IFMISS, CHIA,
     +   NDFCHI, CHIAP)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE ACTUALLY LISTS THE AUTOCORRELATIONS OR
C     PARTIAL AUTOCORRELATIONS AND OTHER PERTINENT INFORMATION.
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
     +   CHIA,CHIAP
      INTEGER
     +   LAGMAX,NDFCHI
      LOGICAL
     +   IFMISS
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   RHO(*),SDRHO(*)
      INTEGER
     +   NLPPA(*)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   FPLM
      INTEGER
     +   I,IMAX,IMIN,IPRT,LAG,NPERL
C
C  EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   D1MACH
      EXTERNAL D1MACH
C
C  EXTERNAL SUBROUTINES
      EXTERNAL IPRINT
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MIN
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION CHIA, CHIAP
C        THE VARIABLES IN CHICH THE CHI SQUARE STATISTIC AND
C        PROBABILITY FOR THE AUTOCORRELATIONS ARE STORED.
C     DOUBLE PRECISION FPLM
C        THE FLOATING POINT LARGEST MAGNITUDE.
C     INTEGER I
C        AN INDEX VARIABLE.
C     LOGICAL IFMISS
C        THE INDICATOR VARIABLE USED TO DETERMINE
C        WHETHER THE INPUT SERIES HAS MISSING DATA OR NOT.
C     INTEGER IMAX, IMIN
C        THE INDEX VALUES OF THE FIRST AND LAST OBSERVATION
C        TO BE PRINTED PER LINE
C     INTEGER IPRT
C        THE UNIT NUMBER FOR PRINTED OUTPUT
C     INTEGER LAG
C        THE LAG VALUE OF THE AUTOCORRELATION OR PARTIAL
C        AUTOCORRELATION BEING PRINTED.
C     INTEGER LAGMAX
C        THE NUMBER OF AUTOCORRELATIONS OR PARTIAL AUTOCORRELATIONS
C        TO BE PRINTED.
C     INTEGER NDFCHI
C        THE DEGREES OF FREEDOM FOR THE CHI SQUARED STATISTIC.
C     INTEGER NLPPA(LAGMAX)
C        THE ARRAY IN WHICH THE NUMBER OF LAGGED PRODUCT PAIRS USED TO
C        COMPUTE EACH AUTOCORRELATION IS STORED
C     INTEGER NPERL
C        THE NUMBER OF VALUES TO BE PRINTED PER LINE.
C     DOUBLE PRECISION RHO(LAGMAX)
C        THE ARRAY IN WHICH THE AUTOCORRELATIONS ARE STORED.
C     DOUBLE PRECISION SDRHO(LAGMAX)
C        THE ARRAY IN WHICH THE STANDARD ERRORS OF THE AUTOCORRELATIONS
C        ARE STORED
C
C
      FPLM = D1MACH(2)
C
      CALL IPRINT(IPRT)
C
      NPERL = 12
      DO 30 I = 1, LAGMAX, NPERL
         IMIN = I
         IMAX = MIN(I + NPERL - 1, LAGMAX)
         WRITE(IPRT, 1000) (LAG, LAG = IMIN, IMAX)
         WRITE(IPRT, 1001) (RHO(LAG), LAG = IMIN, IMAX)
         WRITE(IPRT, 1002) (SDRHO(LAG), LAG = IMIN, IMAX)
         IF (IFMISS) WRITE(IPRT, 1003) (NLPPA(LAG), LAG = IMIN, IMAX)
   30 CONTINUE
C
      IF (SDRHO(LAGMAX) .EQ. FPLM) WRITE(IPRT, 1004) FPLM
C
      WRITE (IPRT, 1005) CHIA, NDFCHI, CHIAP
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT(/19H LAG               , 12(1X, I6))
 1001 FORMAT( 19H ACF               , 12(2X, F5.2))
 1002 FORMAT( 19H STANDARD ERROR    , 12(2X, F5.2))
 1003 FORMAT( 19H NO. OF OBS. USED  , 12(1X, I6))
 1004 FORMAT(///5X, F5.2, 38H INDICATES VALUE COULD NOT BE COMPUTED,
     +   ' DUE TO MISSING DATA.')
 1005 FORMAT(///33H THE CHI SQUARE TEST STATISTIC OF/
     +   40H THE NULL HYPOTHESIS OF WHITE NOISE    =, G21.4/
     +   40H DEGREES OF FREEDOM                    =, I17/
     +   40H OBSERVED SIGNIFICANCE LEVEL           =, F17.4)
      END
