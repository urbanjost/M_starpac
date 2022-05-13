*UASVS
      SUBROUTINE UASVS (ACOV, LAGMAX, Y, N, IAR, PHI, LAG, NF,
     +   FMIN, FMAX, NPRT, SPCA, SPCF, FREQ, LDSTAK)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE USER CALLABLE ROUTINE FOR AUTOREGRESSIVE
C     SPECTRUM ESTIMATION WHEN THE ACVF HAVE PREVIOUSLY BEEN
C     COMPUTED AND STORED (LONG CALL).
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
      DOUBLE PRECISION
     +   FMAX,FMIN
      INTEGER
     +   IAR,LAG,LAGMAX,LDSTAK,N,NF,NPRT
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   ACOV(*),FREQ(*),PHI(*),SPCA(*),SPCF(*),Y(*)
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  ARRAYS IN COMMON
      DOUBLE PRECISION DSTAK(12)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   ALPHA,DELTA,FMN,FMX,VAR,YMEAN
      INTEGER
     +   AIC,FTEST,IA,IFP,IO,IPRT,ISORT,ISYM,LACOV,LAIC,LDSMIN,
     +   LPCV,LPHI,LSPC,LWORK,NALL0,WORK,XAXIS,YAXIS
C
C  LOCAL ARRAYS
      DOUBLE PRECISION
     +   RSTAK(12)
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
      EXTERNAL AMEAN,IPRINT,LDSCMP,PARZEN,STKCLR,STKSET,UASDV,UASER,
     +   UASVAR
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MAX,MIN
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
C     DOUBLE PRECISION ACOV(LAGMAX+1)
C        THE ARRAY OF AUTOCOVARIANCE ESTIMATES.
C     INTEGER AIC
C        THE STARTING LOCATION IN THE STACK FOR
C        THE ARRAY CONTAINING THE AIC.
C     DOUBLE PRECISION ALPHA
C        THE DESIRED CONFIDENCE LEVEL.
C     DOUBLE PRECISION DELTA
C        THE SAMPLING INTERVAL.
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     DOUBLE PRECISION FMAX, FMIN
C        THE MAXIMUM AND MINIMUM FREQUENCY FOR WHICH THE
C        SPECTRUM ESTIMATES ARE TO BE COMPUTED.
C     DOUBLE PRECISION FMN, FMX
C        THE MAXIMUM AND MINIMUM FREQUENCY ACTUALLY USED.
C     DOUBLE PRECISION FREQ(NF)
C        THE VECTOR OF FREQUENCIES AT WHICH THE SPECTRUM IS TO BE
C        ESTIMATED.
C     INTEGER FTEST
C        THE STARTING LOCATION IN THE STACK FOR
C        THE ARRAY IN WHICH THE F RATIO AND PROBABILITY ARE STORED.
C     INTEGER IA
C        A VARIABLE USED TO DETERMINE THE AMOUNT OF STORAGE REQUIRED,
C        BASED ON WHETHER OR NOT THE MODEL ORDER IS TO BE SELECTED OR
C        HAS BEEN PROVIDED.
C     INTEGER IAR
C        THE ORDER OF THE AUTOREGRESSIVE PROCESS CHOSEN.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED
C        IF IERR .NE. 0, ERRORS HAVE BEEN DETECTED
C     INTEGER IFP
C        AN INDICATOR FOR STACK ALLOCATION TYPE, WHERE IFP=3 INDICATES
C        REAL AND IFP=4 INDICATES DOUBLE PRECISION.
C     INTEGER IO
C        A VARIABLE USED TO DETERMINE THE AMOUNT OF STORAGE REQUIRED,
C        BASED ON PRINTED OUTPUT REQUESTED.
C     INTEGER IPRT
C        THE LOGICAL UNIT USED FOR PRINTED OUTPUT.
C     INTEGER ISORT
C        THE STARTING LOCATION IN ISTAK FOR
C        AN ARRAY USED FOR SORTING.
C     INTEGER ISTAK(12)
C        THE INTEGER VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER ISYM
C        THE STARTING LOCATION IN ISTAK FOR
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
C        THE ACTUAL LENGTH OF THE WORK ARRAY.
C     INTEGER N
C        THE INTEGER NUMBER OF OBSERVATIONS IN THE SERIES.
C     INTEGER NALL0
C        THE NUMBER OF STACK ALLOCATIONS OUTSTANDING WHEN THIS ROUTINE
C        WAS CALLED.
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
C     DOUBLE PRECISION PHI(LAGMAX)
C        THE ARRAY OF AUTOREGRESSIVE COEFFICIENTS FOR THE
C        SELECTED ORDER.
C     DOUBLE PRECISION RSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     DOUBLE PRECISION SPCA(NF)
C        THE ARAY CONTAINING THE AUTOREGRESSIVE SPECTRUM ESTIMATES.
C     DOUBLE PRECISION SPCF(NF)
C        THE ARRAY CONTAINING THE FOURIER SPECTRUM ESTIMATES.
C     DOUBLE PRECISION VAR
C        THE ONE STEP PREDICTION VARIANCE.
C     INTEGER WORK
C        THE STARTING LOCATION IN THE STACK FOR
C        THE WORK ARRAY.
C     INTEGER XAXIS
C        THE STARTING LOCATION IN RSTAK FOR
C        THE X AXIS VALUES FOR THE SPECTRUM PLOT.
C     DOUBLE PRECISION Y(N)
C        THE ARRAY CONTAINING THE OBSERVED TIME SERIES.
C     INTEGER YAXIS
C        THE STARTING LOCATION IN RSTAK FOR
C        THE Y AXIS VALUES FOR THE SPECTRUM PLOT.
C     DOUBLE PRECISION YMEAN
C        THE MEAN OF THE OBSERVED TIME SERIES
C
C     SET UP NAME ARRAYS
C
      DATA
     +  NMSUB(1),  NMSUB(2),  NMSUB(3),  NMSUB(4),  NMSUB(5),  NMSUB(6)
     + /     'U',       'A',       'S',       'V',       'S',       ' '/
C
      IFP = 4
C
C     SET UP FOR ERROR CHECKING
C
      OPTION(1) = .FALSE.
      OPTION(2) = .FALSE.
      OPTION(3) = .TRUE.
      OPTION(4) = .TRUE.
C
      LACOV = LAGMAX+1
C
      IO = 1
      IF (NPRT .EQ. 0) IO = 0
      IA = 1
      IF (IAR .NE. 0) IA = 0
C
      CALL LDSCMP(6, 0, IO*(2*NF+5), 0, 0, 0, 'D',
     +   LAGMAX + 1 + IA*(3*LAGMAX+1) + IO*(4*NF+10), LDSMIN)
C
C     CALL ERROR CHECKING ROUTINE
C
      CALL UASER(NMSUB, N, ACOV, IAR, PHI, LAGMAX, LAG, LACOV,
     +   NF, LDSTAK, LDSMIN, N, N, OPTION)
C
    5 IF (IERR.NE.0) THEN
         IERR = 1
         CALL IPRINT (IPRT)
         WRITE (IPRT, 1000)
         RETURN
      END IF
C
C     SET SIZE OF WORK AREA.
C
      CALL STKSET (LDSTAK, 4)
C
C     SAVE NUMBER OF OUTSTANDING STACK ALLOCATIONS.
C
      NALL0 = STKST(1)
C
C     SET VARIOUS PROGRAM PARAMETERS.
C
      LSPC = NF
      LPCV = 2*NF + 5
      LPHI = LAGMAX
      LWORK = LAGMAX+1
C
      FMN = MAX(FMIN, 0.0D0)
      FMX = MIN(FMAX, 0.5D0)
      IF (FMN.GE.FMX) THEN
        FMN = 0.0D0
        FMX = 0.5D0
      END IF
C
      ALPHA = 0.95D0
      DELTA = 1.0D0
C
      IF (IAR.GE.1) THEN
C
C     USER HAS CHOSEN ORDER AND SUPPLIED COEFFICIENTS.
C     COMPUTE RESIDUAL VARIANCE.
C
         CALL AMEAN (Y, N, YMEAN)
         CALL UASVAR (Y, YMEAN, N, IAR, PHI, VAR)
      END IF
C
C     SET UP ADDITIONAL STACK WORK AREA, IF NEEDED.
C
      WORK = STKGET(LWORK,IFP)
      IF (IAR.EQ.0) THEN
         LAIC = LAGMAX+1
         AIC = STKGET(LAIC, IFP)
         FTEST = STKGET(2*LAGMAX, IFP)
      ELSE
         LAIC = LWORK
         AIC = WORK
         FTEST = WORK
      END IF
      IF (NPRT.EQ.0) THEN
         XAXIS = WORK
         YAXIS = WORK
         ISYM = WORK
         ISORT = WORK
      ELSE
         XAXIS = STKGET(LPCV, IFP)
         YAXIS = STKGET(LPCV, IFP)
         ISYM = STKGET(LPCV, 2)
         ISORT = ISYM
      END IF
C
      IF (IERR.EQ.1) GO TO 5
C
C     CALL THE MAIN DRIVER FOR AUTOREGRESSIVE SPECTRUM ROUTINES.
C
      CALL UASDV(ACOV, SPCA, SPCF, LSPC, IAR, PHI, NF, FMN,
     +   FMX, FREQ, N, LAGMAX, RSTAK(FTEST), RSTAK(AIC), RSTAK(WORK),
     +   LACOV, LWORK, DELTA, ISTAK(ISORT), ISTAK(ISYM), RSTAK(XAXIS),
     +   RSTAK(YAXIS), LPCV, ALPHA, LAG, LAIC, LPHI, NPRT, VAR, PARZEN,
     +   NMSUB)
C
      CALL STKCLR(NALL0)
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (/42H THE CORRECT FORM OF THE CALL STATEMENT IS//
     +  '       CALL UASVS (ACOV, LAGMAX, Y, N,'/
     +  '      +            IAR, PHI, LAG, NF, FMIN, FMAX, NPRT,'/
     +  '      +            SPCA, SPCF, FREQ, LDSTAK)')
      END
