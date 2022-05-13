*UASV
      SUBROUTINE UASV (ACOV, LAGMAX, N)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE USER CALLABLE ROUTINE FOR AUTOREGRESSIVE
C     SPECTRUM ESTIMATION WHEN THE ACVF HAVE PREVIOUSLY BEEN
C     COMPUTED AND STORED (SHORT CALL).
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
     +   LAGMAX,N
C
C  ARRAY ARGUMENTS
      REAL
     +   ACOV(*)
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      REAL
     +   ALPHA,DELTA,FMAX,FMIN,VAR
      INTEGER
     +   IAR,IPRT,LACOV,LAG,LAIC,LDSMIN,LDSTAK,LPCV,LPHI,LSPC,
     +   LWORK,NF,NPRT
C
C  LOCAL ARRAYS
      REAL
     +   AIC(101),FREQ(101),FTEST(2,100),PHI(100),SPCA(101),SPCF(101),
     +   WORK(101),XAXIS(207),YAXIS(207)
      INTEGER
     +   ISORT(101),ISYM(207)
      LOGICAL
     +   OPTION(4)
      CHARACTER
     +   NMSUB(6)*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL IPRINT,PARZEN,UASDV,UASER
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     REAL ACOV(LAGMAX+1)
C        THE AUTOCOVARIANCE COMPUTED FROM THE LAG PRODUCT PAIRS.
C     REAL AIC(101)
C        THE ARRAY CONTANING AKIAKES CRITERIA FOR EACH ORDER.
C     REAL ALPHA
C        THE DESIRED CONFIDENCE LEVEL.
C     REAL DELTA
C        THE SAMPLING INTERVAL.
C     REAL FMAX, FMIN
C        THE MAXIMUM AND MINIMUM FREQUENCY FOR WHICH THE
C        SPECTRUM ESTIMATES ARE TO BE COMPUTED.
C     REAL FREQ(101)
C        THE VECTOR OF FREQUENCIES AT WHICH THE SPECTRUM IS TO BE
C        ESTIMATED.
C     REAL FTEST(2, 100)
C        THE ARRAY CONTAINING THE F RATIO AND F TEST.
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
C        THE LENGTH OF THE WORK ARRAY.
C     INTEGER N
C        THE INTEGER NUMBER OF OBSERVATIONS IN EACH SERIES
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
C     REAL PHI(100)
C        THE ARRAY OF AUTOREGRESSIVE COEFFICIENTS FOR THE
C        SELECTED ORDER.
C     REAL SPCA(101)
C        THE ARAY CONTAINING THE AUTOREGRESSIVE SPECTRUM ESTIMATES.
C     REAL SPCF(101)
C        THE ARRAY CONTAINING THE FOURIER SPECTRUM ESTIMATES.
C     REAL VAR
C        THE ONE STEP PREDICTION VARIANCE.
C     REAL WORK(101)
C        A WORK AREA USED FOR THE LAG WINDOWS AND FOR
C        COMPUTING THE AUTOREGRESSIVE COEFFICIENTS.
C     REAL XAXIS(207)
C        THE X AXIS VALUES FOR THE SPECTRUM PLOT.
C     REAL YAXIS(207)
C        THE Y AXIS VALUES FOR THE SPECTRUM PLOT.
C
C     SET UP NAME ARRAYS
C
      DATA
     +  NMSUB(1),  NMSUB(2),  NMSUB(3),  NMSUB(4),  NMSUB(5),  NMSUB(6)
     + /     'U',       'A',       'S',       'V',       ' ',       ' '/
C
C     SET UP FOR ERROR CHECKING
C
      OPTION(1) = .FALSE.
      OPTION(2) = .FALSE.
      OPTION(3) = .TRUE.
      OPTION(4) = .FALSE.
C
      LAG = 0
      IAR = 0
      NF = 101
      FMIN = 0.0E0
      FMAX = 0.5E0
      NPRT = -1
      LACOV = LAGMAX+1
      LDSTAK = 0
      LDSMIN = 0
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
C     SET VARIOUS PROGRAM PARAMETERS
C
      LPCV = 207
      LAIC = 101
      LSPC = 101
      LPHI = 100
      LWORK = 101
C
      ALPHA = 0.95E0
      DELTA = 1.0E0
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
     +  '       CALL UASV (ACOV, LAGMAX, N)')
      END
