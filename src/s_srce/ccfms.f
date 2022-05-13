*CCFMS
      SUBROUTINE CCFMS (YM, YMMISS, N, M, IYM, LAGMAX, CCOV,
     +   CMISS, ICCOV, JCCOV, NLPPC, INLPPC, JNLPPC, NPRT, LDSTAK)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE USER CALLABLE ROUTINE FOR COMPUTING THE CROSS
C     CORRELATIONS AND COVARIANCES OF A MULTIVARIATE SERIES (LONG CALL).
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
      REAL
     +   CMISS
      INTEGER
     +   ICCOV,INLPPC,IYM,JCCOV,JNLPPC,LAGMAX,LDSTAK,M,N,NPRT
C
C  ARRAY ARGUMENTS
      REAL
     +   CCOV(*),YM(*),YMMISS(*)
      INTEGER
     +   NLPPC(*)
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
     +   FPLM
      INTEGER
     +   I,IFP,IO,IPRT,IYMFFT,J,K,LDSMIN,LGLST,LGLSTJ,LGLSTK,LYFFT,
     +   NALL0,NFFT,RHOC,SDRHOC,YMEAN,YMEANJ,YMEANK,YSD,YSDJ,YSDK
      LOGICAL
     +   ISFFT,ISLONG,NEWPG
C
C  LOCAL ARRAYS
      REAL
     +   RSTAK(12)
      INTEGER
     +   ISTAK(12)
      CHARACTER
     +   NMSUB(6)*1
C
C  EXTERNAL FUNCTIONS
      REAL
     +   R1MACH
      INTEGER
     +   STKGET,STKST
      EXTERNAL R1MACH,STKGET,STKST
C
C  EXTERNAL SUBROUTINES
      EXTERNAL ACVFM,CCFER,CCFMNM,CCFOUT,IPRINT,LDSCMP,STKCLR,STKSET
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MIN,SQRT
C
C  COMMON BLOCKS
      COMMON /CSTAK/DSTAK
      COMMON /ERRCHK/IERR
C
C  STATEMENT FUNCTIONS
      INTEGER
     +   I2,I3C,I3N
C
C  EQUIVALENCES
      EQUIVALENCE (DSTAK(1),ISTAK(1))
      EQUIVALENCE (DSTAK(1),RSTAK(1))
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     REAL CCOV(ICCOV, JCCOV, M)
C        THE CROSS COVARIANCE MATRIX.
C     REAL CMISS
C        THE MISSING VALUE CODE FOR THE CCVF ESTIMATES.
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     REAL FPLM
C        THE FLOATING POINT LARGEST MAGNITUDE.
C     INTEGER ICCOV
C        THE ACTUAL FIRST DIMENSION OF THE ARRAY CCOV, AS
C        SPECIFIED IN THE USERS PROGRAM.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED
C        IF IERR .EQ. 1, ERRORS HAVE BEEN DETECTED
C     INTEGER IFP
C        AN INDICATOR FOR STACK ALLOCATION TYPE, WHERE IFP=3 INDICATES
C        SINGLE PRECISION AND IFP=4 INDICATES DOUBLE PRECISION.
C     INTEGER INLPPC
C        THE ACTUAL FIRST DIMENSION OF THE ARRAY NLPPC AS SPECIFIEC
C        IN THE USERS PROGRAM.
C     INTEGER IO
C        A VARIABLE USED TO DETERMINE THE AMOUNT OF STORAGE REQUIRED,
C        BASED ON WHETHER PRINTED OUTPUT IS DESIRED.
C     INTEGER IPRT
C        THE UNIT NUMBER USED FOR OUTPUT.
C     LOGICAL ISFFT
C        THE INDICATOR VARIABLE USED TO DESIGNATE WHETHER THE CALLING
C        ROUTINE HAS SUFFIX F (ISFFT = TRUE) OR NOT (ISFFT = FALSE)
C     LOGICAL ISLONG
C        THE INDICATOR VARIABLE USED TO DESIGNATE WHETHER THE CALLING
C        ROUTINE HAS SUFFIX S (ISLONG = TRUE) OR NOT (ISLONG = FALSE)
C     INTEGER ISTAK(12)
C        THE INTEGER VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER IYM
C        THE ACTUAL FIRST DIMENSION OF THE MATRIX YM AS
C        SPECIFIED IN THE USERS PROGRAM.
C     INTEGER IYMFFT
C        THE ACTUAL FIRST DIMENSION OF THE MATRIX YMFFT AS
C        SPECIFIED IN THE USERS PROGRAM.
C     INTEGER I2
C        STATEMENT FUNCTION SPECIFYING THE DESIRED LOCATION WITHIN YM
C     INTEGER I3C
C        STATEMENT FUNCTION SPECIFYING THE DESIRED LOCATION WITHIN CCOV
C     INTEGER I3N
C        STATEMENT FUNCTION SPECIFYING THE DESIRED LOCATION WITHIN NLPPC
C     INTEGER J
C        THE INDEX OF -SERIES 1- IN THE ARRAY YM.
C     INTEGER JCCOV
C        THE ACTUAL SECOND DIMENSION OF THE ARRAY CCOV, AS
C        SPECIFIED IN THE USERS PROGRAM.
C     INTEGER JNLPPC
C        THE SECOND DIMENSION OF THE ARRAY NLPPC AS SPECIFIED
C        IN THE USERS PROGRAM.
C     INTEGER K
C        THE INDEX OF -SERIES 2- IN THE ARRAY YM.
C     INTEGER LAGMAX
C        THE MAXIMUM LAG VALUE REQUESTED.
C     INTEGER LDSMIN
C        THE MINIMUM LENGTH ALLOWED FOR THE ARRAY DSTAK.
C     INTEGER LDSTAK
C        THE LENGTH OF THE ARRAY DSTAK.
C     INTEGER LGLST, LGLSTJ, LGLSTK
C        THE STARTING LOCATION IN ISTAK FOR THE ARRAY LGLST,
C        AND THE LOCATIONS IN LGLST IN EHICH THE NUMBER OF THE
C        LAG OF THE LAST ACVF WHICH COULD BE COMPUTED FOR SERIES
C        J AND K, RESPECTIVELY, BEFORE A MISSNG ACVF (DUE TO MISSING
C        DATA).
C     INTEGER LYFFT
C        THE LENGTH OF THE VECTOR YFFT.
C     INTEGER M
C        THE NUMBER OF SERIES BEING COMPARED, IE THE
C        NUMBER OF COLUMNS OF DATA IN YM.
C     INTEGER N
C        THE INTEGER NUMBER OF OBSERVATIONS IN EACH SERIES
C     INTEGER NALL0
C        THE NUMBER OF OUTSTANDING STACK ALLOCATIONS
C     LOGICAL NEWPG
C        AN INDICATOR VARIABLE USED TO DETERMINE WHEN A NEW PAGE
C        IS APPROPRIATE FOR THE OUTPUT.
C     INTEGER NFFT
C        THE NUMBER OF OBSERVATIONS IN THE EXTENDED SERIES.
C     INTEGER NLPPC(INLPPC, JNLPPC, M)
C        THE ARRAY CONTAINING THE NUMBER OF LAGGED PRODUCT PAIRS
C        USED TO COMPUTE EACH ACVF ESTIMATE.
C     CHARACTER*1 NMSUB(6)
C        THE NAME OF THE SUBROUTINE CALLING THE ERROR CHECKING
C        SUBROUTINE.
C     INTEGER NPRT
C        THE INDICATOR VARIABLE USED TO SPEDIFY WHETHER OR NOT
C        PRINTED OUTPUT IS TO BE GIVEN, WHERE IF THE VALUE OF
C        NPRT IS ZERO, NO OUTPUT IS GIVEN.
C     INTEGER RHOC
C        THE STARTING LOCATION IN DSTAK OF THE ARRAY RHOC.
C     REAL RSTAK(12)
C        THE REAL VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER SDRHOC
C        THE STARTING LOCATION IN THE WORK AREA FOR SDRHOC.
C     REAL YM(IYM, M)
C        THE MATRIX CONTAINING THE OBSERVED TIME SERIES
C     INTEGER YMEAN, YMEANJ, YMEANK
C        THE STARTING LOCATION FOR THE ARRAY CONTAINING THE
C        MEANS OF EACH OF THE SERIES, AND THE LOCATION IN
C        THE ARRAY FOR THE MEAN OF THE JTH AND KTH SERIES,
C        RESPECTIVELY.
C     REAL YMMISS(M)
C        THE MISSING VALUE CODES FOR EACH OF THE SERIES IN YM.
C     INTEGER YSD, YSDJ, YSDK
C        THE STARTING LOCATION FOR THE ARRAY CONTAINING THE
C        STANDARD DEVIATIONS OF EACH OF THE SERIES, AND THE
C        LOCATION IN THE ARRAY FOR THE STANDARD DEVIATION OF
C        THE JTH AND KTH SERIES, RESPECTIVELY.
C
C     SET UP NAME ARRAYS
C
      DATA
     +  NMSUB(1),  NMSUB(2),  NMSUB(3),  NMSUB(4),  NMSUB(5),  NMSUB(6)
     + /     'C',       'C',       'F',       'M',       'S',       ' '/
C
      I2(I,J) = I + (J-1)*IYM
      I3C(I,J,K) = I + (J-1)*ICCOV + (K-1)*JCCOV*ICCOV
      I3N(I,J,K) = I + (J-1)*INLPPC + (K-1)*JNLPPC*INLPPC
C
C     SET UP FOR ERROR CHECKING
C
      IERR = 0
      LYFFT = N + LAGMAX
      IYMFFT = LYFFT
      NFFT = N
      ISFFT = .FALSE.
      ISLONG = .TRUE.
C
      IO = 1
      IF (NPRT .EQ. 0) IO = 0
C
      CALL LDSCMP(3+2*IO, 0, M, 0, 0, 0, 'S',
     +            2*M+IO*(4*LAGMAX+2), LDSMIN)
C
      CALL CCFER(NMSUB, N, LAGMAX, LDSTAK, LDSMIN, ICCOV, JCCOV,
     +  INLPPC, JNLPPC, M, LYFFT, NFFT, IYM, IYMFFT, ISFFT, ISLONG)
C
C     CHECK WHETHER AN ERROR HAS BEEN DETECTED
C
      IF (IERR.EQ.0) THEN
C
        FPLM = R1MACH(2)
C
C       SET UP THE WORK AREA.
C
        CALL STKSET (LDSTAK, 4)
        NALL0 = STKST(1)

        IFP = 3
C
        YMEAN = STKGET(M, IFP)
        YSD = STKGET(M, IFP)
        LGLST = STKGET(M, 2)
        IF (NPRT.EQ.0) THEN
          RHOC = YSD
          SDRHOC = YSD
        ELSE
          RHOC = STKGET(2*LAGMAX+1, IFP)
          SDRHOC = STKGET(2*LAGMAX+1, IFP)
        END IF
C
        IF (IERR.EQ.0) THEN
C
C         BEGIN LOOP FOR COMPUTATIONS
C
          NEWPG = .FALSE.
C
          DO 40 K = 1, M
C
             YMEANK = YMEAN + K - 1
             YSDK = YSD + K - 1
             LGLSTK = LGLST + K - 1
C
             CALL ACVFM (YM(I2(1,K)), YMMISS(K), N, RSTAK(YMEANK),
     +                   CCOV(I3C(1,K,K)), LAGMAX, ISTAK(LGLSTK),
     +                   NLPPC(I3N(1,K,K)), LAGMAX+1)
             RSTAK(YSDK) = SQRT(CCOV(I3C(1,K,K)) * N / (N - 1))
C
             IF ((K-1).LE.0) GO TO 40
C
             DO 30 J = 1, (K-1)
C
                YMEANJ = YMEAN + J - 1
                YSDJ = YSD + J - 1
                LGLSTJ = LGLST + J - 1
C
                CALL CCFMNM (YM(I2(1,J)), YMMISS(J),
     +                       YM(I2(1,K)), YMMISS(K),
     +                       N,  LAGMAX, 2*LAGMAX+1,
     +                       CCOV(I3C(1,J,J)), CCOV(I3C(1,K,K)),
     +                       CCOV(I3C(1,J,K)), CCOV(I3C(1,K,J)), ICCOV,
     +                       NLPPC(I3N(1,J,J)), NLPPC(I3N(1,K,K)),
     +                       NLPPC(I3N(1,J,K)), NLPPC(I3N(1,K,J)),
     +                       INLPPC,
     +                       RSTAK(YMEANJ), RSTAK(YMEANK), RSTAK(RHOC),
     +                       RSTAK(SDRHOC),
     +                       NPRT, MIN(ISTAK(LGLSTJ), ISTAK(LGLSTK)))
C
                IF (NPRT .EQ. 0) GO TO 30
C
C         CALL ROUTINE TO PRINT OUT CORRELATIONS
C
                CALL CCFOUT (J, RSTAK(YMEANJ), RSTAK(YSDJ), N,
     +                       NLPPC(I3N(1,J,J)), K, RSTAK(YMEANK),
     +                       RSTAK(YSDK), N, NLPPC(I3N(1,K,K)), LAGMAX,
     +                       2*LAGMAX+1, RSTAK(RHOC), RSTAK(SDRHOC),
     +                       .TRUE., NLPPC(I3N(1,J,K)),
     +                       NLPPC(I3N(1,K,J)),
     +                       INLPPC, YMMISS(J), YMMISS(K), NEWPG)
C
                NEWPG = .TRUE.
C
   30        CONTINUE
   40     CONTINUE
C
          CMISS = FPLM
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
     +  '       CALL CCFMS (YM, YMMISS, N, M, IYM,'/
     +  '      +            LAGMAX, CCOV, CMISS, ICCOV, JCCOV,'/
     +  '      +            NLPPC, INLPPC, JNLPPC, NPRT, LDSTAK)')
      END
