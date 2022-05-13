*UFSM
      SUBROUTINE UFSM(Y, YMISS, N)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE USER CALLABLE ROUTINE FOR TIME SERIES FOURIER
C     SPECTRUM ANALYSIS WITH MISSING DATA (SHORT CALL).
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
      DOUBLE PRECISION
     +   YMISS
      INTEGER
     +   N
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   Y(*)
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   ALPHA,DELTA,FMAX,FMIN
      INTEGER
     +   IPRT,ISPCF,LACOV,LAGMAX,LDSMIN,LDSTAK,LNLPPA,LPCV,LWORK,
     +   LY,NF,NPRT,NW
C
C  LOCAL ARRAYS
      DOUBLE PRECISION
     +   ACOV(101),FREQ(101),SPCF(101,4),WORK(101),XAXIS(106),
     +   YAXIS(106)
      INTEGER
     +   ISORT(101),ISYM(106),LAGS(4),NLPPA(101)
      LOGICAL
     +   OPTION(4)
      CHARACTER
     +   NMSUB(6)*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL ECVF,IPRINT,PARZEN,SETLAG,UFSDRV
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION ACOV(101)
C        THE AUTOCOVARIANCE AT LAG ZERO (BIASED VARIANCE).
C     DOUBLE PRECISION ALPHA
C        THE DESIRED CONFIDENCE LEVEL.
C     DOUBLE PRECISION DELTA
C        THE SAMPLING INTERVAL.
C     DOUBLE PRECISION FMAX, FMIN
C        THE MAXIMUM AND MINIMUM FREQUENCES AT WHICH THE
C        SPECTRUM IS TO BE COMPUTED.
C     DOUBLE PRECISION FREQ(101)
C        THE VECTOR OF FREQUENCIES AT WHICH THE SPECTRUM IS TO BE
C        COMPUTED.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED
C        IF IERR .EQ. 1, ERRORS HAVE BEEN DETECTED
C     INTEGER IPRT
C        THE LOGICAL UNIT USED FOR PRINTED OUTPUT.
C     INTEGER ISORT(101)
C        AN ARRAY USED FOR SORTING.
C     INTEGER ISPCF
C         THE ACTUAL FIRST DIMENSION OF THE SPECTRUM ARRAYS.
C     INTEGER ISYM(106)
C        THE ARRAY CONTAINING THE CODE FOR THE PLOT SYMBOLS.
C     INTEGER LACOV
C        THE LENGTH OF THE VECTOR ACOV.
C     INTEGER LAGMAX
C        THE MAXIMUM LAG VALUE TO BE USED.
C     INTEGER LAGS(4)
C        THE ARRAY USED TO STORE THE LAG WINDOW TRUCCATION
C        POINTS USED FOR EACH SET OF SPECTRUM VALUES.
C     INTEGER LDSMIN
C        THE MINIMUM LENGTH ALLOWED FOR DSTAK.
C     INTEGER LDSTAK
C        THE LENGTH OF THE VECTOR DSTAK IN COMMON CSTAK.
C     INTEGER LNLPPA
C        THE LENGTH OF THE VECTOR NLPPA.
C     INTEGER LPCV
C        THE LENGTH OF THE VECTORS USED FOR PLOTTING.
C     INTEGER LWORK
C        THE LENGTH OF THE VECTOR W.
C     INTEGER LY
C        THE LENGTH OF THE VECTOR Y.
C     INTEGER N
C        THE INTEGER NUMBER OF OBSERVATIONS IN EACH SERIES
C     INTEGER NF
C        THE NUMBER OF FREQUENCIES AT WHICH THE SPECTRUM IS
C        TO BE COMPUTED.
C     INTEGER NLPPA(101)
C        THE ARRAY CONTAINING THE NUMBER OF LAG PRODUCT PAIRS.
C     CHARACTER*1 NMSUB(6)
C        THE ARRAY CONTAINING THE NAME OF THIS SUBROUTINE.
C     INTEGER NPRT
C        A CODE USED TO SPECIFY THE TYPE OF PLOT, WHERE IF
C        NPRT < 0 THE PLOT IS DECIBLES/LINEAR
C        NPRT = 0 THE PLOT IS SUPPRESSED
C        NPRT > 0 THE PLOT IS LOG/LINEAR
C     INTEGER NW
C        THE VARIABLE USED TO DETERMINE THE NUMBER OF DIFFERENT
C        BANDWIDTHS TO BE USED.
C     LOGICAL OPTION(4)
C        AN INDICATOR ARRAY USED TO DESIGNATE WHETHER ANY OF THE
C        FOUR POSSIBLE OPTIONS (F, M, V, OR S) HAVE BEEN USED (TRUE)
C        OR NOT (FALSE).
C     EXTERNAL PARZEN
C        THE SUBROUTINE USED TO COMPUTE THE WINDOW.
C     DOUBLE PRECISION SPCF(101,4)
C        THE ARRAYS IN WHICH THE SPECTRUM IS STORED.
C     DOUBLE PRECISION WORK(101)
C        THE VECTOR OF LAG WINDOWS.
C     DOUBLE PRECISION XAXIS(106)
C        THE X AXIS VALUES FOR THE SPECTRUM PLOT.
C     DOUBLE PRECISION Y(N)
C        THE ARRAY CONTAINING THE OBSERVED TIME SERIES.
C     DOUBLE PRECISION YAXIS(106)
C        THE Y AXIS VALUES FOR THE SPECTRUM PLOT.
C     DOUBLE PRECISION YMISS
C        THE USER SUPPLIED CODE WHICH IS USED TO DETERMINE WHETHER OR
C        NOT AN OBSERVATION IN THE SERIES IS MISSING.  IF Y(I) = YMISS,
C        THE VALUE IS ASSUMED MISSING, OTHERWISE IT IS NOT.
C
C     SET UP NAME ARRAYS
C
      DATA NMSUB(1), NMSUB(2), NMSUB(3), NMSUB(4), NMSUB(5), NMSUB(6) /
     +   'U','F','S','M',' ',' '/
C
C     SET UP
C
      OPTION(4) = .FALSE.
      OPTION(3) = .FALSE.
      OPTION(2) = .TRUE.
      OPTION(1) = .FALSE.
C
      LDSTAK = 0
      LDSMIN = 0
C
      ISPCF = 101
      LACOV = 101
      LNLPPA = 101
      LY = N
      LPCV = 106
      LWORK = 101
      NF = 101
C
C     SET MAXIMUM LAG VALUE TO BE USED.
C     AND NUMBER OF LAG WINDOW TRUNCATION POINTS TO USE.
C
      CALL SETLAG(N, LAGMAX)
      NW = 4
C
C     CALL THE CONTROLLING ROUTINE FOR FOURIER SPECTRUM ROUTINES
C     FOR SERIES WITH MISSING DATA.
C
      CALL UFSDRV(Y, LY, YMISS, ACOV, NLPPA, SPCF, ISPCF, NF, FMIN,
     +   FMAX, FREQ, N, NW, LAGMAX, LAGS, WORK, LACOV, LWORK, DELTA,
     +   ISORT, ISYM, XAXIS, YAXIS, LPCV, ALPHA, NPRT, PARZEN, NMSUB,
     +   LDSMIN, LDSTAK, OPTION, LNLPPA, LY)
C
C     CHECK FOR ERRORS
C
      IF (IERR.EQ.0) RETURN
C
      IF (IERR.EQ.2) CALL ECVF(NMSUB)
      IERR = 1
      CALL IPRINT(IPRT)
      WRITE (IPRT,1000)
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (/42H THE CORRECT FORM OF THE CALL STATEMENT IS//
     +   '       CALL UFSM (Y, YMISS, N)')
      END
