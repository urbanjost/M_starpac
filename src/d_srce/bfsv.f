*BFSV
      SUBROUTINE BFSV(CCOV, INDEX1, INDEX2, N, LAGMAX, ICCOV, JCCOV)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE USER CALLABLE ROUTINE FOR TIME SERIES BIVARIATE
C     FOURIER SPECTRUM ANALYSIS OF SERIES WITH
C     COVARIANCES INPUT RATHER THAN ORIGINAL SERIES
C     (SHORT CALL)
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
      INTEGER
     +   ICCOV,INDEX1,INDEX2,JCCOV,LAGMAX,N
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   CCOV(*)
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   ALPHA,DELTA,FMAX,FMIN,YMISS1,YMISS2
      INTEGER
     +   ICSPC2,INLPPC,IPHAS,IPRT,JNLPPC,LAGMX1,LAGMXU,LDSMIN,
     +   LDSTAK,LPCV,LW,LY,M,NF,NPRT,NW
C
C  LOCAL ARRAYS
      DOUBLE PRECISION
     +   CEVEN(101),CODD(101),CSPC2(101,4),FREQ(101),PHAS(101,4),
     +   SPCF1(101),SPCF2(101),W(101),XAXIS(404),Y1(1),Y2(1),
     +   YAXIS(404)
      INTEGER
     +   ISYM(404),LAGS(4),NLPPC(1,1,1)
      LOGICAL
     +   OPTION(4)
      CHARACTER
     +   NMSUB(6)*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL BFSDRV,IPRINT,PARZEN,SETLAG
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MIN
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION ALPHA
C        THE DESIRED CONFIDENCE LEVEL.
C     DOUBLE PRECISION CCOV(ICCOV,JCCOV,*)
C        THE COVARIANCES.
C     DOUBLE PRECISION CEVEN(101)
C        THE SUMS OF THE AUTOCOVARIANCES FOR EACH LAG.
C     DOUBLE PRECISION CODD(101)
C        THE DIFFERENCES OF THE AUTOCOVARIANCES FOR EACH LAG.
C     DOUBLE PRECISION CSPC2(101,4)
C        THE SQUARED COHERENCY COMPONENT OF THE BIVARIATE SPECTRA.
C     DOUBLE PRECISION DELTA
C        THE SAMPLING INTERVAL.
C     DOUBLE PRECISION FMAX, FMIN
C        THE MAXIMUM AND MINIMUM FREQUENCY AT WHICH THE
C        SPECTRUM IS TO BE COMPUTED.
C     DOUBLE PRECISION FREQ(101)
C        THE VECTOR OF FREQUENCIES AT WHICH THE SPECTRUM IS TO BE
C        COMPUTED.
C     INTEGER ICCOV
C        THE FIRST DIMENSION OF THE ARRAY CCOV.
C     INTEGER ICSPC2
C        THE FIRST DIMENSION OF THE ARRAY CSPC2.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED
C        IF IERR .EQ. 1, ERRORS HAVE BEEN DETECTED
C     INTEGER INDEX1, INDEX2
C        THE INDICES OF THE COVARIANCES OF THE TWO SERIES.
C     INTEGER INLPPC
C        THE FIRST DIMENSION OF THE ARRAY NLPPC.
C     INTEGER IPHAS
C        THE FIRST DIMENSION OF THE ARRAY PHAS.
C     INTEGER IPRT
C        THE LOGICAL UNIT USED FOR PRINTED OUTPUT.
C     INTEGER ISYM(404)
C        THE ARRAY CONTAINING THE CODE FOR THE PLOT SYMBOLS.
C     INTEGER JCCOV
C        THE SECOND DIMENSION OF THE ARRAY CCOV.
C     INTEGER JNLPPC
C        THE SECOND DIMENSION OF THE ARRAY NLPPC.
C     INTEGER LAGMAX, LAGMXU
C        THE MAXIMUM LAG VALUE TO BE USED.
C     INTEGER LAGMX1
C        LAGMAX+1.
C     INTEGER LAGS(4)
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
C     INTEGER LY
C        THE LENGTH OF THE VECTORS Y1 AND Y2.
C     INTEGER M
C        THE NUMBER OF SERIES FOR WHICH THE COVARIANCES WERE
C        COMPUTED
C     INTEGER N
C        THE INTEGER NUMBER OF OBSERVATIONS IN EACH SERIES
C     INTEGER NF
C        THE NUMBER OF FREQUENCIES AT WHICH THE SPECTRUM IS
C        TO BE COMPUTED.
C     INTEGER NLPPC(1,1,1)
C        A DUMMY ARRAY.
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
C     DOUBLE PRECISION PHAS(101,4)
C        THE PHASE COMPONENT OF THE BIVARIATE SPECTRA.
C     DOUBLE PRECISION SPCF1(101), SPCF2(101)
C        THE ARRAYS IN WHICH THE SPECTRUM IS STORED.
C     DOUBLE PRECISION W(101)
C        THE WINDOWS.
C     DOUBLE PRECISION XAXIS(404)
C        THE X AXIS VALUES FOR THE SPECTRUM PLOTS.
C     DOUBLE PRECISION YAXIS(404)
C        THE Y AXIS VALUES FOR THE SPECTRUM PLOTS.
C     DOUBLE PRECISION YMISS1, YMISS2
C        THE MISSING VALUE CODES
C     DOUBLE PRECISION Y1(1)
C        THE FIRST TIME SERIES.
C     DOUBLE PRECISION Y2(1)
C        THE SECOND TIME SERIES.
C
C     SET UP NAME ARRAYS
C
      DATA NMSUB(1), NMSUB(2), NMSUB(3), NMSUB(4), NMSUB(5), NMSUB(6) /
     +   'B','F','S','V',' ',' '/
C
C     SET UP FOR ERROR CHECKING
C
      OPTION(1) = .FALSE.
      OPTION(2) = .FALSE.
      OPTION(3) = .TRUE.
      OPTION(4) = .FALSE.
C
      M = 2
C
      INLPPC = 1
      JNLPPC = 1
      ICSPC2 = 101
      IPHAS = 101
C
      LDSTAK = 0
      LDSMIN = 0
C
      NF = 101
      LW = 101
      LY = N
      LPCV = 404
C
C     SET MAXIMUM LAG VALUE USED (LAGMXU)
C     SET NUMBER OF LAG WINDOW TRUCCATION POINTS (NW)
C
      CALL SETLAG(N, LAGMXU)
      LAGMXU = MIN(LAGMXU,LAGMAX)
      NW = 4
C
C     CALL THE CONTROLING ROUTINE FOR THE BIVARIATE SPECTRUM ROUTINES
C
      CALL BFSDRV(Y1, Y2, YMISS1, YMISS2, CCOV, NLPPC, SPCF1, SPCF2,
     +   NF, FMIN, FMAX, FREQ, N, NW, LAGMXU, LAGS, LAGMX1, W, LW,
     +   DELTA, ISYM, XAXIS, YAXIS, LPCV, ALPHA, NPRT, PARZEN, ICCOV,
     +   JCCOV, M, INDEX1, INDEX2, CSPC2, PHAS, ICSPC2, IPHAS, CODD,
     +   CEVEN, W, LW, NMSUB, LDSMIN, LDSTAK, OPTION, N, INLPPC,
     +   JNLPPC, LY)
C
      IF (IERR.EQ.0) RETURN
C
      CALL IPRINT(IPRT)
      WRITE (IPRT,1000)
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (/42H THE CORRECT FORM OF THE CALL STATEMENT IS//
     +  '       CALL BFSV(CCOV, INDEX1, INDEX2, N, LAGMAX, ICCOV,',
     +  ' JCCOV)')
      END
