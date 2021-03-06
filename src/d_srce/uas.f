*UAS
      SUBROUTINE UAS (Y, N)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE USER CALLABLE ROUTINE FOR AUTOREGRESSIVE
C     SPECTRUM ESTIMATION (SHORT CALL).
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  DECEMBER 7, 1981
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
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
     +   ALPHA,DELTA,FMAX,FMIN,VAR,YMEAN
      INTEGER
     +   IAR,IPRT,LACOV,LAG,LAGMAX,LAIC,LDSMIN,LDSTAK,LPCV,LPHI,
     +   LSPC,LWORK,NF,NPRT
C
C  LOCAL ARRAYS
      DOUBLE PRECISION
     +   ACOV(101),AIC(101),FREQ(101),FTEST(2,100),PHI(100),SPCA(101),
     +   SPCF(101),WORK(101),XAXIS(207),YAXIS(207)
      INTEGER
     +   ISORT(101),ISYM(207)
      LOGICAL
     +   OPTION(4)
      CHARACTER
     +   NMSUB(6)*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL ACVF,IPRINT,PARZEN,SETLAG,UASDV,UASER
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION ACOV(101)
C        THE AUTOCOVARIANCE COMPUTED FROM THE LAG PRODUCT PAIRS.
C     DOUBLE PRECISION AIC(101)
C        THE ARRAY CONTANING AKIAKES CRITERIA FOR EACH ORDER.
C     DOUBLE PRECISION ALPHA
C        THE DESIRED CONFIDENCE LEVEL.
C     DOUBLE PRECISION DELTA
C        THE SAMPLING INTERVAL.
C     DOUBLE PRECISION FMAX, FMIN
C        THE MAXIMUM AND MINIMUM FREQUENCY FOR WHICH THE
C        SPECTRUM ESTIMATES ARE TO BE COMPUTED.
C     DOUBLE PRECISION FREQ(101)
C        THE VECTOR OF FREQUENCIES AT WHICH THE SPECTRUM IS TO BE
C        ESTIMATED.
C     DOUBLE PRECISION FTEST(2, 100)
C        THE ARRAY IN WHICH THE F RATIO AND PROBABILITY ARE STORED.
C     INTEGER IAR
C        THE ORDER OF THE AUTOREGRESSIVE PROCESS CHOSEN.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED
C        IF IERR .NE. 0, ERRORS HAVE BEEN DETECTED
C     INTEGER IPRT
C        THE LOGICAL UNIT USED FOR PRINTED OUTPUT.
C     INTEGER ISORT(101)
C        AN ARRAY USED FOR SORTING.
C     INTEGER ISYM(207)
C        THE ARRAY CONTAINING THE CODE FOR THE PLOT SYMBOLS.
C     INTEGER LACOV
C        THE LENGTH OF THE COVARIANCE ARRAYS.
C     INTEGER LAG
C        THE LAG WINDOW TRUNCATION POINT USED FOR A SPECIFIC WINDOW.
C     INTEGER LAGMAX
C        THE MAXIMUM LAG VALUE TO BE USED.
C     INTEGER LAIC
C        THE LENGTH OF THE ARRAY AIC.
C     INTEGER LDSMIN
C        THE MINIMUM LENGTH ALLOWED FOR DSTAK.
C     INTEGER LDSTAK
C        THE LENGTH OF THE VECTOR DSTAK IN COMMON CSTAK.
C     INTEGER LPCV
C        THE LENGTH OF THE PLOT CO-ORDINATE VECTORS.
C     INTEGER LPHI
C        THE LENGTH OF THE VECTOR PHI.
C     INTEGER LSPC
C         THE LENGTH OF THE SPECTRUM ARRAYS.
C     INTEGER LWORK
C        THE LENGTH OF THE VECTOR WORK.
C     INTEGER N
C        THE INTEGER NUMBER OF OBSERVATIONS IN THE SERIES.
C     INTEGER NF
C        THE NUMBER OF FREQUENCIES FOR WHICH THE SPECTRUM ESTIMATES
C        ARE TO BE ESTIMATED.
C     CHARACTER*1 NMSUB(6)
C        THE ARRAY CONTAINING THE NAME OF THIS SUBROUTINE.
C     INTEGER NPRT
C        A CODE USED TO SPECIFY THE TYPE OF PLOT, WHERE IF
C        NPRT < 0 THE PLOT IS DECIBLES/LINEAR
C        NPRT = 0 THE PLOT IS SUPPRESSED
C        NPRT > 0 THE PLOT IS LOG/LINEAR
C     LOGICAL OPTION(4)
C        AN INDICATOR ARRAY USED TO DESIGNATE WHETHER ANY OF THE
C        FOUR POSSIBLE OPTIONS (F, M, V, OR S) HAVE BEEN USED (TRUE)
C        OR NOT (FALSE).
C     EXTERNAL PARZEN
C        THE TYPE OF WINDOW TO BE USED.
C     DOUBLE PRECISION PHI(100)
C        THE ARRAY OF AUTOREGRESSIVE COEFFICIENTS FOR THE
C        SELECTED ORDER.
C     DOUBLE PRECISION SPCA(101)
C        THE ARAY CONTAINING THE AUTOREGRESSIVE SPECTRUM ESTIMATES.
C     DOUBLE PRECISION SPCF(101)
C        THE ARRAY CONTAINING THE FOURIER SPECTRUM ESTIMATES.
C     DOUBLE PRECISION VAR
C        THE ONE STEP PREDICTION VARIANCE FOR THE SELECTED MODEL.
C     DOUBLE PRECISION WORK(101)
C        A DOUBLE PRECISION WORK AREA USED FOR THE LAG WINDOWS AND FOR
C        COMPUTING THE AUTOREGRESSIVE COEFFICIENTS.
C     DOUBLE PRECISION XAXIS(207)
C        THE X AXIS VALUES FOR THE SPECTRUM PLOT.
C     DOUBLE PRECISION Y(N)
C        THE ARRAY CONTAINING THE OBSERVED TIME SERIES.
C     DOUBLE PRECISION YAXIS(207)
C        THE Y AXIS VALUES FOR THE SPECTRUM PLOT.
C     DOUBLE PRECISION YMEAN
C        THE MEAN OF THE OBSERVED TIME SERIES
C
C     SET UP NAME ARRAYS
C
      DATA
     +  NMSUB(1),  NMSUB(2),  NMSUB(3),  NMSUB(4),  NMSUB(5),  NMSUB(6)
     + /     'U',       'A',       'S',       ' ',       ' ',       ' '/
C
C     SET UP FOR ERROR CHECKING
C
      OPTION(1) = .FALSE.
      OPTION(2) = .FALSE.
      OPTION(3) = .FALSE.
      OPTION(4) = .FALSE.
C
      LAG = 0
      IAR = 0
      NF = 101
      FMIN = 0.0D0
      FMAX = 0.5D0
      NPRT = -1
      LDSTAK = 0
      LDSMIN = 0
C
C     SET THE MAXIMUM NUMBER OF LAGS TO BE USED.
C
      CALL SETLAG(N, LAGMAX)
C
C     CALL ERROR CHECKING ROUTINE
C
      CALL UASER(NMSUB, N, ACOV, IAR, PHI, LAGMAX, LAG, LACOV,
     +   NF, LDSTAK, LDSMIN, N, N, OPTION)
C
      IF (IERR.NE.0) THEN
         IERR = 1
         CALL IPRINT (IPRT)
         WRITE (IPRT, 1000)
         RETURN
      END IF
C
C     SET VARIOUS PROGRAM PARAMETERS.
C
      LPCV = 207
      LSPC = 101
      LPHI = 100
      LAIC = 101
      LACOV = 101
      LWORK = 101
C
      ALPHA = .95D0
      DELTA = 1.0D0
C
C     COMPUTE AUTOCOVARIANCES
C
      CALL ACVF (Y, N, YMEAN, ACOV, LAGMAX, LACOV)
C
C     CALL THE MAIN DRIVER FOR AUTOREGRESSIVE SPECTRUM ROUTINES.
C
      CALL UASDV(ACOV, SPCA, SPCF, LSPC, IAR, PHI, NF, FMIN, FMAX, FREQ,
     +   N, LAGMAX, FTEST, AIC, WORK, LACOV, LWORK, DELTA, ISORT,
     +   ISYM, XAXIS, YAXIS, LPCV, ALPHA, LAG, LAIC, LPHI, NPRT, VAR,
     +   PARZEN, NMSUB)
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (/42H THE CORRECT FORM OF THE CALL STATEMENT IS//
     +   22H       CALL UAS (Y, N))
      END
