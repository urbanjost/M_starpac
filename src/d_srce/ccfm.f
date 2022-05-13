*CCFM
      SUBROUTINE CCFM (Y1, Y1MISS, Y2, Y2MISS, N)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE USER CALLABLE ROUTINE FOR COMPUTING THE CROSS
C     CORRELATIONS OF TWO TIME SERIES WITH MISSING OBSERVATIONS
C     (SHORT CALL).
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
     +   Y1MISS,Y2MISS
      INTEGER
     +   N
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   Y1(*),Y2(*)
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   Y1MEAN,Y1SD,Y2MEAN,Y2SD
      INTEGER
     +   ICCOV,INLPPC,IPRT,IYM,IYMFFT,JCCOV,JNLPPC,LAGMAX,LDSMIN,
     +   LDSTAK,LGLST1,LGLST2,LYFFT,M,NFFT,NUSED1,NUSED2
      LOGICAL
     +   ISFFT,ISLONG
C
C  LOCAL ARRAYS
      DOUBLE PRECISION
     +   CCOV(101,2,2),RHOC(201),SDRHOC(201)
      INTEGER
     +   NLPPC(101,2,2)
      CHARACTER
     +   NMSUB(6)*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL ACVFM,CCFER,CCFMNM,CCFOUT,IPRINT,SETLAG
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MIN,SQRT
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION CCOV(101, 2, 2)
C        THE ARRAY USED FOR THE CCVF ESTIMATES.
C     INTEGER ICCOV
C        THE ACTUAL FIRST DIMENSION OF THE ARRAY CCOV, AS
C        SPECIFIED IN THE USERS PROGRAM.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED
C        IF IERR .EQ. 1, ERRORS HAVE BEEN DETECTED
C     INTEGER INLPPC
C        THE ACTUAL FIRST DIMENSION OF THE ARRAY NLPPC AS SPECIFIEC
C        IN THE USERS PROGRAM.
C     INTEGER IPRT
C        THE UNIT NUMBER USED FOR OUTPUT.
C     LOGICAL ISFFT
C        THE INDICATOR VARIABLE USED TO DESIGNATE WHETHER THE CALLING
C        ROUTINE HAS SUFFIX F (ISFFT = TRUE) OR NOT (ISFFT = FALSE)
C     LOGICAL ISLONG
C        THE INDICATOR VARIABLE USED TO DESIGNATE WHETHER THE CALLING
C        ROUTINE HAS SUFFIX S (ISLONG = TRUE) OR NOT (ISLONG = FALSE)
C     INTEGER IYM
C        THE ACTUAL FIRST DIMENSION OF THE MATRIX YM AS
C        SPECIFIED IN THE USERS PROGRAM.
C     INTEGER IYMFFT
C        THE ACTUAL FIRST DIMENSION OF THE MATRIX YMFFT AS
C        SPECIFIED IN THE USERS PROGRAM.
C     INTEGER JCCOV
C        THE ACTUAL SECOND DIMENSION OF THE ARRAY CCOV, AS
C        SPECIFIED IN THE USERS PROGRAM.
C     INTEGER JNLPPC
C        THE SECOND DIMENSION OF THE ARRAY NLPPC AS SPECIFIED
C        IN THE USERS PROGRAM.
C     INTEGER LAGMAX
C        THE NUMBER OF AUTOCORRELATIONS DESIRED.
C     INTEGER LDSMIN
C        THE MINIMUM LENGTH ALLOWED FOR THE ARRAY DSTAK.
C     INTEGER LDSTAK
C        THE LENGTH OF THE ARRAY DSTAK.
C     INTEGER LGLST1, LGLST2
C        THE LAG VALUE OF THE LAST ACVF WHICH COULD BE COMPUTED
C        FOR SERIES 1 AND 2, RESPECTIVELY, BEFORE MISSING DATA
C        CAUSED A MISSING ACVF.
C     INTEGER LYFFT
C        THE LENGTH OF THE VECTOR YFFT.
C     INTEGER M
C        THE NUMBER OF SERIES BEING COMPARED, IE THE
C        NUMBER OF COLUMNS OF DATA IN YM.
C     INTEGER N
C        THE INTEGER NUMBER OF OBSERVATIONS IN EACH SERIES
C     INTEGER NFFT
C        THE NUMBER OF OBSERVATIONS IN THE EXTENDED SERIES.
C     INTEGER NLPPC(101, 2, 2)
C        THE NUMBER OF LAGGED PRODUCT PAIRS USED TO COMPUTE THE CCVF.
C     CHARACTER*1 NMSUB(6)
C        THE NAME OF THE SUBROUTINE CALLING THE ERROR CHECKING
C        SUBROUTINE.
C     INTEGER NUSED1, NUSED2
C        THE NUMBER OF ACTIVE (NON MISSING) OBSERVATIONS IN SERIES
C        1 AND 2, RESPECTIVELY.
C     DOUBLE PRECISION RHOC(201)
C        THE CCF ESTIMATES.
C     DOUBLE PRECISION SDRHOC(201)
C        THE ARRAY CONTAINING THE SD OF THE CCFM.
C     DOUBLE PRECISION Y1(N), Y1MEAN, Y1MISS, Y1SD
C        THE FIRST SERIES, AND ITS MEAN, MISSING VALUE CODE AND
C        STANDARD DEVIATION.
C     DOUBLE PRECISION Y2(N), Y2MEAN, Y2MISS, Y2SD
C        THE SECOND SERIES, AND ITS MEAN, MISSING VALUE CODE AND
C        STANDARD DEVIATION.
C
C     SET UP NAME ARRAYS
C
      DATA
     +  NMSUB(1),  NMSUB(2),  NMSUB(3),  NMSUB(4),  NMSUB(5),  NMSUB(6)
     + /     'C',       'C',       'F',       'M',       ' ',       ' '/
C
C     SET UP FOR ERROR CHECKING
C
      IERR = 0
      ICCOV = 101
      INLPPC = 101
      IYM = N
      JCCOV = 2
      JNLPPC = 2
      LDSMIN = 0
      LDSTAK = 0
      LAGMAX = 1
      LYFFT = N + LAGMAX
      IYMFFT = LYFFT
      M = 2
      NFFT = N
      ISFFT = .FALSE.
      ISLONG = .FALSE.
C
C     CALL ERROR CHECKING ROUTINES
C
      CALL CCFER(NMSUB, N, LAGMAX, LDSTAK, LDSMIN, ICCOV, JCCOV,
     +  INLPPC, JNLPPC, M, LYFFT, NFFT, IYM, IYMFFT, ISFFT, ISLONG)
C
C     CHECK WHETHER AN ERROR HAS BEEN DETECTED
C
      IF (IERR.EQ.0) THEN
C
C       SET THE MAXIMUM LAG VALUE TO BE USED.
C
        CALL SETLAG (N, LAGMAX)
C
C       COMPUTE THE SERIES ACVF AND SD
C
        CALL ACVFM (Y1, Y1MISS, N, Y1MEAN, CCOV(1,1,1), LAGMAX, LGLST1,
     +     NLPPC(1,1,1), 101)
C
        NUSED1 = NLPPC(1,1,1)
        Y1SD = SQRT(CCOV(1,1,1) * N / (NUSED1-1))
C
        CALL ACVFM (Y2, Y2MISS, N, Y2MEAN, CCOV(1,2,2), LAGMAX, LGLST2,
     +     NLPPC(1,2,2), 101)
C
        NUSED2 = NLPPC(1,2,2)
        Y2SD = SQRT(CCOV(1,2,2) * N / (NUSED2-1))
C
C       CALL ROUTINE FOR MAIN AUTOCORRELATION COMPUTATIONS.
C
        IF (CCOV(1,1,1)*CCOV(1,2,2) .NE. 0.0D0)
     +    CALL CCFMNM (Y1, Y1MISS, Y2, Y2MISS, N, LAGMAX, 201,
     +    CCOV(1,1,1),
     +    CCOV(1,2,2), CCOV(1,1,2), CCOV(1,2,1), 101, NLPPC(1,1,1),
     +    NLPPC(1,2,2), NLPPC(1,1,2), NLPPC(1,2,1), 101, Y1MEAN, Y2MEAN,
     +    RHOC, SDRHOC, 1, MIN(LGLST1, LGLST2))
C
C       CALL ROUTINE TO PRINT OUT AUTOCORRELATIONS
C
        CALL CCFOUT (1, Y1MEAN, Y1SD, N, NUSED1, 2, Y2MEAN, Y2SD, N,
     +     NUSED2, LAGMAX, 201, RHOC, SDRHOC, .TRUE., NLPPC(1,1,2),
     +     NLPPC(1,2,1), 101, Y1MISS, Y2MISS, .FALSE.)
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
     +   44H       CALL CCFM (Y1, YMISS1, Y2, YMISS2, N))
      END
