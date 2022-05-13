*UASDV
      SUBROUTINE UASDV(ACOV, SPCA, SPCF, LSPC, IAR, PHI, NF, FMIN, FMAX,
     +   FREQ, N, LAGMAX, FTEST, AIC, WORK, LACOV, LWORK, DELTA, ISORT,
     +   ISYM, XAXIS, YAXIS, LPCV, ALPHA, LAG, LAIC, LPHI, NPRT, VAR,
     +   WINDOW, NMSUB)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE MAIN DRIVER FOR COMPUTING THE AUTOREGRESSIVE
C     (AND FOURIER) SPECTRUMS.
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
     +   ALPHA,DELTA,FMAX,FMIN,VAR
      INTEGER
     +   IAR,LACOV,LAG,LAGMAX,LAIC,LPCV,LPHI,LSPC,LWORK,N,NF,NPRT
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   ACOV(LACOV),AIC(LAIC),FREQ(NF),FTEST(2,LAGMAX),PHI(LPHI),
     +   SPCA(LSPC),SPCF(LSPC),WORK(LWORK),XAXIS(LPCV),YAXIS(LPCV)
      INTEGER
     +   ISORT(NF),ISYM(LPCV)
      CHARACTER
     +   NMSUB(6)*1
C
C  SUBROUTINE ARGUMENTS
      EXTERNAL WINDOW
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   ALOW,AUP,BW,DF,SPCAMN,SPCAMX,SPCFMN,SPCFMX,XPLTMN,XPLTMX,
     +   YPLTMN,YPLTMX
      INTEGER
     +   ISPCER,NPTS,NSPCA,NSPCF,NW
      LOGICAL
     +   AICPRT
C
C  LOCAL ARRAYS
      INTEGER
     +   LAGS(1),NLPPA(1)
C
C  EXTERNAL SUBROUTINES
      EXTERNAL AOS,SETFRQ,SPCCK,UASCFT,UASEST,UASORD,UASOUT,UFSLAG,UFSMN
C
C  INTRINSIC FUNCTIONS
      INTRINSIC IABS,INT
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION ACOV(LACOV)
C        THE AUTOCOVARIANCE COMPUTED FROM THE LAG PRODUCT PAIRS.
C     DOUBLE PRECISION AIC(LAIC)
C        THE ARRAY CONTAINING THE AKAIKES CRITERIA FOR EACH ORDER(?).
C     LOGICAL AICPRT
C        AN INDICATOR VARIABLE USED TO DETERMINE IF THE AKIAKE
C        INFORMATION CRITERIA AND CHI SQUARED STATISTICS SHOULD
C        BE PRINTED.
C     DOUBLE PRECISION ALOW
C        A FACTOR USED TO COMPUTE THE LOWER CONFIDENCE LIMITS.
C     DOUBLE PRECISION ALPHA
C        THE DESIRED CONFIDENCE LEVEL.
C     DOUBLE PRECISION AUP
C        A FACTOR USED TO COMPUTE THE UPPER CONFIDENCE LIMITS.
C     DOUBLE PRECISION BW
C        THE BANDWIDTH.
C     DOUBLE PRECISION DELTA
C        THE SAMPLING INTERVAL.
C     DOUBLE PRECISION DF
C        THE EFFECTIVE DEGREES OF FREEDOM.
C     DOUBLE PRECISION FMAX, FMIN
C        THE MAXIMUM AND MINIMUM FREQUENCY FOR WHICH THE
C        SPECTRUM ESTIMATES ARE TO BE COMPUTED.
C     DOUBLE PRECISION FREQ(NF)
C        THE VECTOR OF FREQUENCIES AT WHICH THE SPECTRUM IS TO BE
C        ESTIMATED.
C     DOUBLE PRECISION FTEST(2,LAGMAX)
C        THE ARRAY IN WHICH THE F RATIO AND PROBABILITY ARE STORED.
C        THE ORDER OF THE AUTOREGRESSIVE PROCESS CHOSEN.
C     INTEGER ISORT(NF)
C        AN ARRAY USED FOR SORTING.
C     INTEGER LSPC
C         THE ACTUAL FIRST DIMENSION FOR THE SPECTRUM ARRAYS.
C     INTEGER ISPCER
C        AN ERROR FLAG USED FOR THE SPECTRUM PLOTS.
C     INTEGER ISYM(LPCV)
C        THE ARRAY CONTAINING THE CODE FOR THE PLOT SYMBOLS.
C     INTEGER LACOV
C        THE LENGTH OF THE COVARIANCE ARRAYS.
C     INTEGER LAG
C        THE LAG WINDOW TRUNCATION POINT USED FOR A SPECIFIC WINDOW.
C     INTEGER LAGMAX
C        THE MAXIMUM LAG VALUE TO BE USED.
C     INTEGER LAGS(1)
C        THE LAG WINDOW TRUNCATION POINT RETURNED FROM UFSLAG.
C     INTEGER LAIC
C        THE LENGTH OF THE ARRAY AIC.
C     INTEGER LPCV
C        THE LENGTH OF THE PLOT CO-ORDINATE VECTORS.
C     INTEGER LPHI
C        THE LENGTH OF THE VECTOR PHI.
C     INTEGER LWORK
C        THE ACTUAL LENGTH OF THE WORK ARRAY.
C     INTEGER N
C        THE INTEGER NUMBER OF OBSERVATIONS IN THE SERIES
C     INTEGER NF
C        THE NUMBER OF FREQUENCIES FOR WHICH THE SPECTRUM ESTIMATES
C        ARE TO BE ESTIMATED.
C     INTEGER NLPPA(1)
C        A DUMMY ARRAY
C     CHARACTER*1 NMSUB(6)
C        THE NAME OF THE CALLING SUBROUTINE.
C     INTEGER NPRT
C        A CODE USED TO SPECIFY THE TYPE OF PLOT, WHERE IF
C        NPRT < 0 THE PLOT IS DECIBLES/LINEAR
C        NPRT = 0 THE PLOT IS SUPPRESSED
C        NPRT > 0 THE PLOT IS LOG/LINEAR
C     INTEGER NPTS
C        THE NUMBER OF X, Y CO-ORDINATES TO BE PLOTTED.
C     INTEGER NSPCA, NSPCF
C        THE NUMBER OF VALID SPECTRUM ESTIMATES FOR THE AUTOREGRESSIVE
C        AND FOURIER SPECTRUMS, RESPECTIVELY.
C     INTEGER NW
C        THE NUMBER OF LAG WINDOW TRUNCATION POINTS SELCTED.
C     DOUBLE PRECISION PHI(LPHI)
C        THE ARRAY OF AUTOREGRESSIVE COEFFICIENTS FOR THE
C        SELECTED ORDER.
C     DOUBLE PRECISION SPCA(LSPC)
C        THE ARAY CONTAINING THE AUTOREGRESSIVE SPECTRUM ESTIMATES.
C     DOUBLE PRECISION SPCAMN, SPCAMX
C        THE MINIMUM AND MAXIMUM AUTOREGRESSIVE SPECTRUM VALUE TO BE
C        PLOTTED.
C     DOUBLE PRECISION SPCF(LSPC)
C        THE ARRAY CONTAINING THE FOURIER SPECTRUM ESTIMATES.
C     DOUBLE PRECISION SPCFMN, SPCFMX
C        THE MINIMUM AND MAXIMUM FOURIER SPECTRUM VALUE TO BE PLOTTED.
C     DOUBLE PRECISION VAR
C        THE ONE STEP PREDICTION VARIANCE.
C     EXTERNAL WINDOW
C        THE TYPE OF WINDOW TO BE USED.
C     DOUBLE PRECISION WORK(LWORK)
C        THE WORK ARRAY.
C     DOUBLE PRECISION XAXIS(LPCV)
C        THE X AXIS VALUES FOR THE SPECTRUM PLOT.
C     DOUBLE PRECISION XPLTMN, XPLTMX
C        THE MINIMUM AND MAXIMUM VALUES TO BE PLOTTED FOR THE X AXIS.
C     DOUBLE PRECISION YAXIS(LPCV)
C        THE Y AXIS VALUES FOR THE SPECTRUM PLOT.
C     DOUBLE PRECISION YPLTMN, YPLTMX
C        THE MINIMUM AND MAXIMUM VALUES TO BE PLOTTED FOR THE Y AXIS.
C
      NW = 1
C
      IF (LAG.LE.0) THEN
C
C     SET THE LAG WINDOW TRUNCATION POINT TO BE USED FOR THE
C     FOURIER SPECTRUM ESTIMATES.
C
         CALL UFSLAG(ACOV, LAGMAX, LAGS, N, NW, NW, LACOV)
         LAG = LAGS(1)/2
      END IF
C
C     SET FREQUENCIES FOR THE SPECTRUM.
C
      CALL SETFRQ(FREQ, NF, 1, FMIN, FMAX, DELTA)
C
C     COMPUTE THE FOURIER SPECTRUM ESTIMATES
C
      CALL UFSMN(ACOV, NLPPA, LAG, DF, NF, FREQ, ALPHA, BW, SPCF,
     +  ALOW, AUP, LACOV, LSPC, WINDOW, WORK, LWORK, N, DELTA,
     +  .FALSE., 1)
C
      AICPRT = .FALSE.
C
      IF (IAR.LT.0) THEN
C
C     USER HAS CHOSEN ORDER.
C     COMPUTE COEFFICIENTS AND VARIANCE USING DURBINS RECURSIVE METHOD.
C
         CALL UASCFT(ACOV, LAGMAX, LACOV, IABS(IAR), PHI, N, VAR)
C
      ELSE IF (IAR.EQ.0) THEN
C
C     SELECT MODEL ORDER AND COMPUTE COEFFICIENTS AND VARIANCE.
C
         AICPRT = .TRUE.
         CALL AOS(N, LAGMAX, ACOV, WORK, IAR, VAR, PHI,
     +            WORK, AIC, FTEST, LACOV, LAIC)
      END IF
C
C     COMPUTE THE AUTOREGRESSIVE SPECTRUM ESTIMATES.
C
      CALL UASEST(IABS(IAR), VAR, PHI, NF, FREQ, DELTA, SPCA, LPHI,
     +   LSPC)
C
      IF (NPRT.EQ.0) RETURN
C
C     SET PLOTTING VECTORS.
C
      XPLTMN = FMIN
      XPLTMX = FMAX
C
      YPLTMN = 0.0D0
      YPLTMX = 0.0D0
C
      CALL SPCCK(SPCF, ISORT, NF, SPCFMN, SPCFMX, NSPCF, ISPCER)
      IF (ISPCER.NE.0) GO TO 40
      CALL SPCCK(SPCA, ISORT, NF, SPCAMN, SPCAMX, NSPCA, ISPCER)
      IF (ISPCER.NE.0) GO TO 40
C
      CALL UASORD(SPCF, SPCA, SPCFMN, SPCFMX, SPCAMN, SPCAMX, FREQ, NF,
     +   XAXIS, YAXIS, ISYM, NPTS, LSPC, LPCV, NSPCF, NSPCA, BW, ALOW,
     +   AUP, XPLTMN, XPLTMX, YPLTMN, YPLTMX, NPRT)
C
C     PRINT RESULTS
C
   40 CALL UASOUT(XAXIS, YAXIS, ISYM, NPTS, BW, INT(DF+0.5D0), LAG,
     +   IABS(IAR), PHI, ISPCER, LPCV, XPLTMN, XPLTMX, YPLTMN, YPLTMX,
     +   FTEST, AIC, LAIC, VAR, NPRT, LAGMAX, AICPRT, N, NMSUB)
C
      RETURN
      END
