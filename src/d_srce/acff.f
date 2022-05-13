*ACFF
      SUBROUTINE ACFF (YFFT, N, LYFFT, LDSTAK)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE USER CALLABLE ROUTINE FOR COMPUTING THE AUTO-
C     CORRELATIONS AND PARTIAL AUTOCORRELATIONS OF A TIME SERIES
C     USING AN FFT (SHORT CALL).
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
     +   LDSTAK,LYFFT,N
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   YFFT(*)
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
     +   CHIA,CHIAP,OSPVAR,YMEAN,YSD
      INTEGER
     +   IAR,IFP,IPRT,LACOV,LAGMAX,LAIC,LDSMIN,NALL0,NFAC,NFFT,
     +   SDRHO,WORK
      LOGICAL
     +   DIFFER,ISFFT
C
C  LOCAL ARRAYS
      DOUBLE PRECISION
     +   ACOV(101),AIC(101),FTEST(2,100),PHI(100),PRHO(100),RHO(100),
     +   RSTAK(12)
      INTEGER
     +   IOD(1),ND(1),NDUM(1)
      CHARACTER
     +   NMSUB(6)*1
C
C  EXTERNAL FUNCTIONS
      INTEGER
     +   STKGET,STKST
      EXTERNAL STKGET,STKST
C
C  EXTERNAL SUBROUTINES
      EXTERNAL ACFER,ACFMNF,ACFOUT,FFTLEN,IPRINT,LDSCMP,SETLAG,STKCLR,
     +   STKSET
C
C  INTRINSIC FUNCTIONS
      INTRINSIC SQRT
C
C  COMMON BLOCKS
      COMMON /CSTAK/DSTAK
      COMMON /ERRCHK/IERR
C
C  EQUIVALENCES
      EQUIVALENCE (DSTAK(1),RSTAK(1))
      EQUIVALENCE (ACOV(2),RHO(1))
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION ACOV(101)
C        THE AUTOCOVARIANCE VECTOR.
C     DOUBLE PRECISION AIC(101)
C       THE ARRAY CONTAINING AKAIKES CRITERIA FOR EACH ORDER.
C     DOUBLE PRECISION CHIA, CHIAP
C        THE VARIABLES IN WHICH THE CHI SQUARE STATISTIC AND
C        CHI SQUARED STATISTIC PROBABILITY FOR THE AUTOCORRELATIONS
C        ARE STORED.
C     LOGICAL DIFFER
C        THE INDICATOR VARIABLE USED TO DESIGNATE WHETHER THE CALLING
C        ROUTINE IS ACFD (DIFFER = TRUE) OR NOT (DIFFER = FALSE)
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     DOUBLE PRECISION FTEST(2, 100)
C        THE ARRAY IN WHICH THE F RATIO AND PROBABILITY ARE STORED.
C     INTEGER IAR
C        THE ORDER OF THE AUTOREGRESSIVE PROCESS CHOSEN.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED
C        IF IERR .NE. 0 ERRORS WERE DETECTED.
C     INTEGER IFP
C        THE INDICATOR VARIABLE USED TO DESIGNATE WHETHER THE FLOATING
C        POINT VARIABLES ARE SINGLE (IFP=3) OR DOUBLE (IFP=4) PRECISION.
C     INTEGER IOD(1)
C        THE ORDER OF EACH OF THE DIFFERENCE FACTORS.
C     INTEGER IPRT
C        THE UNIT NUMBER USED FOR OUTPUT.
C     LOGICAL ISFFT
C        THE INDICATOR VARIABLE USED TO DESIGNATE WHETHER THE CALLING
C        ROUTINE HAS SUFFIX F (ISFFT = TRUE) OR NOT (ISFFT = FALSE)
C     INTEGER LACOV
C        THE LENGTH OF THE VECTOR ACOV.
C     INTEGER LAGMAX
C        THE MAXIMUM LAG VALUE REQUESTED.
C     INTEGER LAIC
C        THE LENGTH OF THE VECTOR AIC.
C     INTEGER LDSMIN
C        THE MINIMUM LENGTH ALLOWED FOR THE ARRAY DSTAK.
C     INTEGER LDSTAK
C        THE LENGTH OF THE ARRAY DSTAK.
C     INTEGER LYFFT
C        THE LENGTH OF THE VECTOR YFFT.
C     INTEGER N
C        THE INTEGER NUMBER OF OBSERVATIONS IN EACH SERIES
C     INTEGER NALL0
C        THE NUMBER OF OUTSTANDING STACK ALLOCATIONS
C     INTEGER ND(1)
C        THE NUMBER OF TIMES EACH DIFFERENCE FACTOR IS TO BE APPLIED
C     INTEGER NDUM(1)
C        A DUMMY ARRAY.
C     INTEGER NFAC
C        THE NUMBER OF DIFFERENCE FACTORS
C     INTEGER NFFT
C        THE NUMBER OF OBSERVATIONS IN THE EXTENDED SERIES.
C     CHARACTER*1 NMSUB(6)
C        THE NAME OF THE SUBROUTINE CALLING THE ERROR CHECKING
C        SUBROUTINE.
C     DOUBLE PRECISION OSPVAR
C        THE ONE STEP PREDICTION VARIANCE FOR THE SELECTED ORDER (IAR).
C     DOUBLE PRECISION PHI(100)
C        THE ARRAY OF AUTOREGRESSIVE COEFFICIENTS FOR THE SELECTED
C        ORDER.
C     DOUBLE PRECISION PRHO(100)
C        THE ARRAY CONTAINING THE PARITAL ACF ESTIMATES.
C     DOUBLE PRECISION RHO(100)
C        THE ARRAY CONTAINING THE ACF ESTIMATES.
C     DOUBLE PRECISION RSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER SDRHO
C        THE STARTING LOCATION IN DSTAK FOR
C        THE ARRAY CONTAINING THE STANDARD ERRORS OF THE ACF ESTIMATES.
C     INTEGER WORK
C        THE STARTING LOCATION IN THE WORK AREA FOR WORK.
C     DOUBLE PRECISION YFFT(LYFFT)
C        THE VECTOR CONTAINING THE OBSERVED TIME SERIES
C     DOUBLE PRECISION YMEAN
C        THE MEAN OF THE OBSERVED TIME SERIES
C     DOUBLE PRECISION YSD
C        THE STANDARD DEVIATION OF THE OBSERVED TIME SERIES
C
C
C     SET UP NAME ARRAYS
C
      DATA
     +  NMSUB(1),  NMSUB(2),  NMSUB(3),  NMSUB(4),  NMSUB(5),  NMSUB(6)
     + /     'A',       'C',       'F',       'F',       ' ',       ' '/
C
C     SET UP FOR ERROR CHECKING
C
      IERR = 0
      LACOV = 101
      LAIC = 101
      DIFFER = .FALSE.
      NFAC = 1
      ND(1) = 0
      IOD(1) = 0
      ISFFT = .TRUE.
C
      IF (N.GE.3) THEN
C
C     SET NUMBER OF ACF TO BE COMPUTED
C     AND LENGTH OF EXTENDED SERIES
C
         CALL SETLAG(N, LAGMAX)
         CALL FFTLEN(N+LAGMAX, 4, NFFT)
      END IF
C
      CALL LDSCMP(1, 0, 0, 0, 0, 0, 'D', NFFT, LDSMIN)
C
C     CALL ERROR CHECKING ROUTINES
C
      CALL ACFER(NMSUB, N, LAGMAX, LACOV, LDSTAK, LDSMIN,
     +  DIFFER, NFAC, ND, IOD, ISFFT, LYFFT, NFFT)
C
C     CHECK WHETHER AN ERROR HAS BEEN DETECTED
C
      IF (IERR.EQ.0) THEN
C
C       SET UP THE WORK AREA.
C
        CALL STKSET (LDSTAK, 4)
        NALL0 = STKST(1)
C
        IFP = 4
C
        WORK = STKGET(NFFT, IFP)
        SDRHO = WORK
C
        IF (IERR.EQ.0) THEN
C
C         CALL ROUTINE FOR MAIN AUTOCORRELATION COMPUTATIONS.
C
          CALL ACFMNF (YFFT, N, NFFT, LAGMAX, RHO, RSTAK(SDRHO), YMEAN,
     +       PRHO, AIC, FTEST, PHI, IAR, OSPVAR, ACOV, LACOV, LAIC,
     +       CHIA, CHIAP, LYFFT, RSTAK(WORK), NFFT, 1)
C
          YSD = SQRT (ACOV(1) * N / (N - 1))
C
C         CALL ROUTINE TO PRINT OUT AUTOCORRELATIONS
C
          CALL ACFOUT(YMEAN, YSD, N, N, LAGMAX, RHO, RSTAK(SDRHO), PRHO,
     +       NDUM, AIC, LAIC, FTEST, IAR, PHI, OSPVAR, CHIA, CHIAP,
     +       LAGMAX, .FALSE., 0.0D0, .FALSE., .FALSE., 0, NDUM, NDUM,
     +       0)
        END IF
C
        CALL STKCLR(NALL0)
      END IF
C
      IF (IERR.NE.0) THEN
C
C     PRINT PROPER CALL SEQUENCE AND RETURN
C
        IERR = 1
        CALL IPRINT (IPRT)
        WRITE (IPRT, 1000)
      END IF
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT(/42H THE CORRECT FORM OF THE CALL STATEMENT IS//
     +  '       CALL ACFF (YFFT, N, LYFFT, LDSTAK)')
      END
