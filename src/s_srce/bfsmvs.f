*BFSMVS
      SUBROUTINE BFSMVS(CCOV, NLPPC, INDEX1, INDEX2, N, ICCOV, JCCOV,
     +   INLPPC, JNLPPC, NW, LAGS, NF, FMIN,
     +   FMAX, NPRT, CSPC2, ICSPC2, PHAS, IPHAS, FREQ, LDSTAK)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE USER CALLABLE ROUTINE FOR TIME SERIES BIVARIATE
C     FOURIER SPECTRUM ANALYSIS OF SERIES WITH MISSING OBSERVATIONS
C     WITH USER INPUT OF THE COVARIANCES RATHER THAN THE SERIES
C     (LONG CALL)
C
C     WRITTEN BY - JANET R. DONALDSON
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  DECEMBER 2, 1985
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      REAL
     +   FMAX,FMIN
      INTEGER
     +   ICCOV,ICSPC2,INDEX1,INDEX2,INLPPC,IPHAS,JCCOV,JNLPPC,
     +   LDSTAK,N,NF,NPRT,NW
C
C  ARRAY ARGUMENTS
      REAL
     +   CCOV(*),CSPC2(*),FREQ(*),PHAS(*)
      INTEGER
     +   LAGS(*),NLPPC(*)
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  ARRAYS IN COMMON
      DOUBLE PRECISION DSTAK(12)
C
C  LOCAL SCALARS
      REAL
     +   ALPHA,DELTA,YMISS1,YMISS2
      INTEGER
     +   CEVEN,CODD,I,IFP,IO,IPRT,ISYM,LAGMAX,LAGMX1,LDSMIN,LPCV,
     +   LW,LWORK,LY,M,NALL0,SPCF1,SPCF2,W,WORK,XAXIS,YAXIS
C
C  LOCAL ARRAYS
      REAL
     +   RSTAK(12),Y1(1),Y2(1)
      INTEGER
     +   ISTAK(12)
      LOGICAL
     +   OPTION(4)
      CHARACTER
     +   NMSUB(6)*1
C
C  EXTERNAL FUNCTIONS
      INTEGER
     +   STKGET,STKST
      EXTERNAL STKGET,STKST
C
C  EXTERNAL SUBROUTINES
      EXTERNAL BFSDRV,ECVF,IPRINT,LDSCMP,PARZEN,STKCLR,STKSET
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MAX
C
C  COMMON BLOCKS
      COMMON /CSTAK/DSTAK
      COMMON /ERRCHK/IERR
C
C  EQUIVALENCES
      EQUIVALENCE (DSTAK(1),RSTAK(1))
      EQUIVALENCE (DSTAK(1),ISTAK(1))
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     REAL ALPHA
C        THE DESIRED CONFIDENCE LEVEL.
C     REAL CCOV(ICCOV,JCCOV,*)
C        THE COVARIANCES.
C     INTEGER CEVEN
C        THE STARTING LOCATION IN THE WORK AREA FOR
C        THE SUMS OF THE AUTOCOVARIANCES FOR EACH LAG.
C     INTEGER CODD
C        THE STARTING LOCATION IN THE WORK AREA FOR
C        THE DIFFERENCES OF THE AUTOCOVARIANCES FOR EACH LAG.
C     REAL CSPC2(ICSPC2,NW)
C        THE SQUARED COHERENCY COMPONENT OF THE BIVARIATE SPECTRA.
C     REAL DELTA
C        THE SAMPLING INTERVAL.
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     REAL FMAX, FMIN
C        THE MAXIMUM AND MINIMUM FREQUENCY AT WHICH THE
C        SPECTRUM IS TO BE COMPUTED.
C     REAL FREQ(NF)
C        THE VECTOR OF FREQUENCIES AT WHICH THE SPECTRUM IS TO BE
C        COMPUTED.
C     INTEGER I
C        AN INDEX VALUE.
C     INTEGER ICCOV
C        THE FIRST DIMENSION OF THE ARRAY CCOV.
C     INTEGER ICSPC2
C        THE FIRST DIMENSION OF THE ARRAY CSPC2.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED
C        IF IERR .EQ. 1, ERRORS HAVE BEEN DETECTED
C     INTEGER IFP
C        AN INDICATOR FOR STACK ALLOCATION TYPE, WHERE IFP=3 INDICATES
C        REAL AND IFP=4 INDICATES DOUBLE PRECISION.
C     INTEGER INDEX1, INDEX2
C        THE INDICES OF THE COVARIANCES OF THE TWO SERIES.
C     INTEGER IO
C        A VARIABLE USED TO DETERMINE THE AMOUNT OF STORAGE REQUIRED,
C        BASED ON PRINTED OUTPUT REQUESTED.
C     INTEGER IPHAS
C        THE FIRST DIMENSION OF THE ARRAY PHAS.
C     INTEGER IPRT
C        THE LOGICAL UNIT USED FOR PRINTED OUTPUT.
C     INTEGER ISTAK(12)
C        THE INTEGER VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER ISYM
C        THE STARTING LOCATION IN THE WORK AREA FOR
C        THE ARRAY CONTAINING THE CODE FOR THE PLOT SYMBOLS.
C     INTEGER JCCOV
C        THE SECOND DIMENSION OF THE ARRAY CCOV.
C     INTEGER LAGMAX
C        THE MAXIMUM LAG VALUE TO BE USED.
C     INTEGER LAGMX1
C        LAGMAX+1.
C     INTEGER LAGS(NW)
C        THE ARRAY USED TO STORE THE LAG WINDOW TRUCCATION
C        POINTS USED FOR EACH SET OF SPECTRUM VALUES.
C     INTEGER LDSMIN
C        THE MINIMUM LENGTH ALLOWED FOR DSTAK.
C     INTEGER LDSTAK
C        THE LENGTH OF THE VECTOR DSTAK IN COMMON CSTAK.
C     INTEGER LPCV
C        THE LENGTH OF THE PLOT CO-ORDINATE VECTORS.
C     INTEGER LW
C        THE LENGTH OF THE VECTOR W.
C     INTEGER LWORK
C        THE LENGTH OF THE VECTOR WORK.
C     INTEGER LY
C        THE LENGTH OF THE VECTORS Y1 AND Y2.
C     INTEGER M
C        THE NUMBER OF SERIES FOR WHICH THE COVARIANCES WERE
C        COMPUTED
C     INTEGER N
C        THE INTEGER NUMBER OF OBSERVATIONS IN EACH SERIES
C     INTEGER NALL0
C        THE NUMBER OF STACK ALLOCATIONS ON ENTRY.
C     INTEGER NF
C        THE NUMBER OF FREQUENCIES AT WHICH THE SPECTRUM IS
C        TO BE COMPUTED.
C     INTEGER NLPPC(INLPPC,JNLPPC,*)
C        THE NUMBER OF OBSERVATIONS IN EACH COVARIANCE ESTIMATE
C     CHARACTER*1 NMSUB(6)
C        THE ARRAY CONTAINING THE NAME OF THIS SUBROUTINE.
C     INTEGER NPRT
C        A CODE USED TO SPECIFY THE TYPE OF PLOT.
C        IF NPRT < 0 THE PLOT IS DECIBLES/LINEAR
C        IF NPRT = 0 THE PLOT IS SUPPRESSED.
C        IF NPRT > 0 THE PLOT IS LOG/LINEAR
C     INTEGER NW
C        THE ARGUMENT USED TO DETERMINE THE NUMBER OF DIFFERENT
C        BANDWIDTHS TO BE USED.
C     LOGICAL OPTION(4)
C        AN INDICATOR ARRAY USED TO DESIGNATE WHETHER ANY OF THE
C        FOUR POSSIBLE OPTIONS (F, M, V, OR S) HAVE BEEN USED (TRUE)
C        OR NOT (FALSE).
C     EXTERNAL PARZEN
C        THE SUBROUTINE USED TO COMPUTE THE WINDOW.
C     REAL PHAS(IPHAS,NW)
C        THE PHASE COMPONENT OF THE BIVARIATE SPECTRA.
C     REAL RSTAK(12)
C        THE REAL VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER SPCF1, SPCF2
C        THE STARTING LOCATION IN THE WORK AREA FOR
C        THE ARRAYS IN WHICH THE SPECTRUM IS STORED.
C     INTEGER W
C        THE STARTING LOCATION IN THE WORK AREA FOR
C        THE WINDOWS.
C     INTEGER WORK
C        THE STARTING LOCATION IN THE WORK AREA FOR THE VECTOR WORK.
C     INTEGER XAXIS
C        THE STARTING LOCATION IN THE WORK AREA FOR
C        THE X AXIS VALUES FOR THE SPECTRUM PLOTS.
C     INTEGER YAXIS
C        THE STARTING LOCATION IN THE WORK AREA FOR
C        THE Y AXIS VALUES FOR THE SPECTRUM PLOTS.
C     REAL YMISS1, YMISS2
C        DUMMY VARIABLES.
C     REAL Y1(1)
C        THE FIRST TIME SERIES.
C     REAL Y2(1)
C        THE SECOND TIME SERIES.
C
C     SET UP NAME ARRAYS
C
      DATA NMSUB(1), NMSUB(2), NMSUB(3), NMSUB(4), NMSUB(5), NMSUB(6) /
     +   'B','F','S','M','V','S'/
C
C     SET UP FOR ERROR CHECKING
C
      OPTION(1) = .FALSE.
      OPTION(2) = .TRUE.
      OPTION(3) = .TRUE.
      OPTION(4) = .TRUE.
C
C     SET MAXIMUM LAG VALUE TO BE USED (LAGMAX).
C
      LAGMAX = N - 1
      IF (NW.LE.0) GO TO 20
      LAGMAX = LAGS(1)
      DO 10 I=1,NW
         LAGMAX = MAX(LAGMAX,LAGS(I))
   10 CONTINUE
   20 CONTINUE
      LAGMX1 = LAGMAX + 1
C
      M = 2
C
C     COMPUTE THE MINIMUM ALLOWABLE STACK AREA (LDSMIN)
C
      IO = 1
      IF (NPRT.EQ.0) IO = 0
C
      CALL LDSCMP(8, 0, IO*4*NF, 0, 0, 0, 'S',
     +  7*LAGMAX+7+2*NF+IO*8*NF, LDSMIN)
C
      LY = N
      LPCV = 4*NF
      LW = LAGMAX + 1
C
C     SET SIZE OF WORK AREA.
C     SET THE NUMBER OF OUTSTANDING STACK ALLOCATIONS (NALL0).
C     SET THE STACK ALLOCATION TYPE (IFP)
C
      CALL STKSET(LDSTAK, 4)
      NALL0 = STKST(1)
C
      IFP = 3
C
C     SET STARTING LOCATIONS IN THE WORK AREA FOR VARIOUS ARRAYS
C
      IF ((LDSMIN.GT.LDSTAK) .OR. (LDSMIN.LE.6)) THEN
         CEVEN = 1
         CODD = 1
         SPCF1 = 1
         SPCF2 = 1
         W = 1
C
         ISYM = 1
         XAXIS = 1
         YAXIS = 1
      ELSE
         CEVEN = STKGET(LAGMX1,IFP)
         CODD = STKGET(LAGMX1,IFP)
         SPCF1 = STKGET(NF,IFP)
         SPCF2 = STKGET(NF,IFP)
         W = STKGET(LW,IFP)
         IF (NPRT.EQ.0) THEN
            ISYM = W
            XAXIS = W
            YAXIS = W
         ELSE
            ISYM = STKGET(LPCV,2)
            XAXIS = STKGET(LPCV,IFP)
            YAXIS = STKGET(LPCV,IFP)
         END IF
      END IF
C
      WORK = W
      LWORK = LW
C
C     CALL THE CONTROLING ROUTINE FOR THE BIVARIATE SPECTRUM ROUTINES
C
      CALL BFSDRV(Y1, Y2, YMISS1, YMISS2, CCOV, NLPPC,
     +   RSTAK(SPCF1), RSTAK(SPCF2), NF, FMIN, FMAX, FREQ, N, NW,
     +   LAGMAX, LAGS, LAGMX1, RSTAK(WORK), LWORK, DELTA, ISTAK(ISYM),
     +   RSTAK(XAXIS), RSTAK(YAXIS), LPCV, ALPHA, NPRT, PARZEN, ICCOV,
     +   JCCOV, M, INDEX1, INDEX2, CSPC2, PHAS, ICSPC2, IPHAS,
     +   RSTAK(CODD), RSTAK(CEVEN), RSTAK(W), LW, NMSUB, LDSMIN,
     +   LDSTAK, OPTION, N, INLPPC, JNLPPC, LY)
C
      CALL STKCLR(NALL0)
C
C     CHECK FOR ERRORS
C
      IF (IERR.NE.0) THEN
        IF (IERR.EQ.2) CALL ECVF(NMSUB)
        IERR = 1
        CALL IPRINT(IPRT)
        WRITE (IPRT,1000)
      END IF
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (/42H THE CORRECT FORM OF THE CALL STATEMENT IS//
     +   '       CALL BFSMVS(CCOV, NLPPC, INDEX1, INDEX2, N,'/
     +   '      +            ICCOV, JCCOV, INLPPC, JNLPPC,'/
     +   '      +            NW, LAGS, NF, FMIN, FMAX, NPRT,'/
     +   '      +            CSPC2, ICSPC2, PHAS, IPHAS, FREQ, LDSTAK)')
      END
