*CCFLST
      SUBROUTINE CCFLST (RHOC, SDRHOC, NLPP12, NLPP21, LAGMAX, LCCOV,
     +   NCC, IFMISS)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE ACTUALLY LISTS THE CROSS CORRELATIONS AND THEIR
C     STANDARD ERRORS, AND MISCELLANEOUS SUMMARY INFORMATION.
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
      INTEGER
     +   LAGMAX,LCCOV,NCC
      LOGICAL
     +   IFMISS
C
C  ARRAY ARGUMENTS
      REAL
     +   RHOC(*),SDRHOC(*)
      INTEGER
     +   NLPP12(*),NLPP21(*)
C
C  LOCAL SCALARS
      REAL
     +   FPLM
      INTEGER
     +   I,I1,IMAX,IMIN,IPRT,K,K0,K1,LAGN,NPERL
C
C  LOCAL ARRAYS
      REAL
     +   RLST(12),SDRLST(12)
      INTEGER
     +   LAG(12),NLPLST(12)
C
C  EXTERNAL FUNCTIONS
      REAL
     +   R1MACH
      EXTERNAL R1MACH
C
C  EXTERNAL SUBROUTINES
      EXTERNAL IPRINT
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MIN,MOD
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     REAL FPLM
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
C     INTEGER I1
C        AN INDEX VARIABLE.
C     INTEGER K, K0, K1
C        INDEX VARIABLES.
C     INTEGER LAG(12)
C        THE LAG VALUE OF THE CORRELATION BEING PRINTED.
C     INTEGER LAGMAX
C        THE LARGEST LAG VALUE TO BE USED.
C     INTEGER LAGN
C        THE NUMBER OF LAG VALUES TO BE PRINTED PER LINE.
C     INTEGER LCCOV
C        THE NUMBER OF LOCATIONS ALLOWED FOR STORING THE NLPPC.
C     INTEGER NCC
C        THE NUMBER OF CROSS CORRELATIONS COMPUTED (FROM -LAGMAX
C        TO +LAGMAX).
C     INTEGER NLPLST(12)
C        THE ARRAY WHICH CONTAINS THE VALUES OF NLPPC TO BE PRINTED
C        ON EACH LINE, ORDERED PROPERLY.
C     INTEGER NLPP12(LCCOV), NLPP21(LCCOV)
C        THE NUMBER OF LAGGED PRODUCT PAIRS USED TO COMPUTE EACH
C        CCVF AT EACH LAG.
C     INTEGER NPERL
C        THE NUMBER OF VALUES TO BE PRINTED PER LINE.
C     REAL RHOC(NCC)
C        THE ARRAY IN WHICH THE AUTOCORRELATIONS OR PARTIAL
C        AUTOCORRELATIONS WILL BE PASSED TO THIS ROUTINE.
C     REAL RLST(12)
C        THE ARRAY WHICH CONTAINS THE VALUES OF RHO TO BE PRINTED
C        ON EACH LINE, ORDERED PROPERLY.
C     REAL SDRHOC(NCC)
C        THE ARRAY IN WHICH THE STANDARD ERRORS OF THE AUTOCORRELATIONS
C        ARE STORED
C     REAL SDRLST(12)
C        THE ARRAY WHICH CONTAINS THE VALUES OF SDRHO TO BE PRINTED
C        ON EACH LINE, ORDERED PROPERLY.
C
C
      CALL IPRINT(IPRT)
      NPERL = 12
C
      K0 = LAGMAX + 1
C
      LAGN = MOD(LAGMAX, NPERL)
      IF (LAGN .EQ. 0) LAGN = NPERL
      I1 = LAGN + 1
C
      DO 20 I = I1, K0, NPERL
         DO 10 K = 1, LAGN
            LAG(K) = I - K0 - K
            K1 = I - K
            RLST(K) = RHOC(K1)
            SDRLST(K) = SDRHOC(K1)
            IF (.NOT. IFMISS) GO TO 10
            K1 = K0 - K1
            NLPLST(K) = NLPP21(K1+1)
   10    CONTINUE
         WRITE(IPRT, 1000) (LAG(K), K = 1, LAGN)
         WRITE(IPRT, 1001) (RLST(K), K = 1, LAGN)
         WRITE(IPRT, 1002) (SDRLST(K), K = 1, LAGN)
         IF (IFMISS) WRITE(IPRT, 1003) (NLPLST(K), K = 1, LAGN)
         LAGN = NPERL
   20 CONTINUE
C
      LAG(1) = 0
      WRITE(IPRT, 1000) LAG(1)
      WRITE(IPRT, 1001) RHOC(K0)
      WRITE(IPRT, 1002) SDRHOC(K0)
      IF (IFMISS) WRITE(IPRT, 1003) NLPP12(1)
C
      DO 40 I = 1, LAGMAX, NPERL
         IMIN = I + K0
         IMAX = MIN(IMIN + NPERL - 1, 2*LAGMAX+1)
         LAGN = IMAX - IMIN + 1
         DO 30 K = 1, LAGN
            LAG(K) = I - 1 + K
   30    CONTINUE
         WRITE(IPRT, 1000) (LAG(K), K = 1, LAGN)
         WRITE(IPRT, 1001) (RHOC(K), K = IMIN, IMAX)
         WRITE(IPRT, 1002) (SDRHOC(K), K = IMIN, IMAX)
         IF (.NOT. IFMISS) GO TO 40
         IMIN = I
         IMAX = MIN(I + NPERL - 1, LAGMAX)
         WRITE (IPRT,1003) (NLPP12(K+1), K=IMIN,IMAX)
   40 CONTINUE
C
      FPLM = R1MACH(2)
C
      IF (SDRHOC(1).EQ.FPLM .OR. SDRHOC(2*LAGMAX+1).EQ.FPLM)
     +   WRITE(IPRT, 1004) FPLM
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT(/19H LAG               , 12(1X, I6))
 1001 FORMAT( 19H CCF               , 12(2X, F5.2))
 1002 FORMAT( 19H STANDARD ERROR    , 12(2X, F5.2))
 1003 FORMAT( 19H NO. OF OBS. USED  , 12(1X, I6))
 1004 FORMAT(///5X, F5.2, 38H INDICATES VALUE COULD NOT BE COMPUTED,
     +   21H DUE TO MISSING DATA.)
      END
