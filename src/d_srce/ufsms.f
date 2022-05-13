*UFSMS
      SUBROUTINE UFSMS(Y, YMISS, N, NW, LAGS, NF, FMIN, FMAX, NPRT,
     +   SPCF, ISPCF, FREQ, LDSTAK)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE USER CALLABLE ROUTINE FOR TIME SERIES FOURIER
C     SPECTRUM ANALYSIS WITH MISSING DATA (LONG CALL).
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
     +   FMAX,FMIN,YMISS
      INTEGER
     +   ISPCF,LDSTAK,N,NF,NPRT,NW
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   FREQ(*),SPCF(*),Y(*)
      INTEGER
     +   LAGS(*)
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
     +   ALPHA,DELTA
      INTEGER
     +   ACOV,I,IFP,IO,IPRT,ISORT,ISYM,LACOV,LAGMAX,LDSMIN,LNLPPA,
     +   LPCV,LWORK,LY,NALL0,NLPPA,WORK,XAXIS,YAXIS
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
      EXTERNAL ECVF,IPRINT,LDSCMP,PARZEN,STKCLR,STKSET,UFSDRV
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
C     INTEGER ACOV
C        THE STARTING LOCATION IN RSTAK FOR THE ACVF VECTOR.
C     DOUBLE PRECISION ALPHA
C        THE DESIRED CONFIDENCE LEVEL.
C     DOUBLE PRECISION DELTA
C        THE SAMPLING INTERVAL.
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     DOUBLE PRECISION FMAX, FMIN
C        THE MAXIMUM AND MINIMUM FREQUENCY AT WHICH THE SPECTRUM
C        IS TO BE COMPUTED.
C     DOUBLE PRECISION FREQ(NF)
C        THE VECTOR OF FREQUENCIES AT WHICH THE SPECTRUM IS TO BE
C        COMPUTED.
C     INTEGER I
C        AN INDEX VARIABLE
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST
C        IF ERR01, NO ERRORS WERE DETECTED
C        IF IERR .EQ. 1, ERRORS HAVE BEEN DETECTED
C     INTEGER IFP
C        AN INDICATOR FOR STACK ALLOCATION TYPE, WHERE IFP=3 INDICATES
C        REAL AND IFP=4 INDICATES DOUBLE PRECISION.
C     INTEGER IO
C        A VARIABLE USED TO DETERMINE THE AMOUNT OF STORAGE REQUIRED
C        BASED ON PRINTED OUTPUT REQUESTED.
C     INTEGER IPRT
C        THE LOGICAL UNIT USED FOR PRINTED OUTPUT.
C     INTEGER ISORT
C        THE STARTING LOCATION FOR THE ARRAY USED FOR SORTING.
C     INTEGER ISPCF
C         THE ACTUAL FIRST DIMENSION OF THE SPECTRUM ARRAYS.
C     INTEGER ISTAK(12)
C        THE INTEGER VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER ISYM
C        THE STARTING LOCATION IN THE WORK AREA FOR ARRAY ISYM.
C     INTEGER LACOV
C        THE LENGTH OF VECTOR ACOV.
C     INTEGER LAGMAX
C        THE MAXIMUM LAG VALUE TO BE USED.
C     INTEGER LAGS(NW)
C        THE ARRAY USED TO SPECIFY THE LAG WINDOW TRUNCATION
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
C        THE INTEGER NUMBER OF OBSERVATIONS IN THE SERIES.
C     INTEGER NALL0
C        THE NUMBER OF ALLOCATIONS OUTSTANDING AT THE TIME THAT
C        THIS ROUTINE WAS CALLED.
C     INTEGER NF
C        THE NUMBER OF FREQUENCIES AT WHICH THE SPECTRUM IS
C        TO BE COMPUTED.
C     INTEGER NLPPA
C        THE STARTING LOCATION IN ISTAK FOR THE ARRAY CONTAINING
C        THE NUMBERS OF LAGGED PRODUCT PAIRS USED FOR EACH ACVF.
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
C     DOUBLE PRECISION RSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     DOUBLE PRECISION SPCF(ISPCF,NW)
C        THE ARRAYS IN WHICH THE SPECTRUM IS STORED
C        FOR EACH LAG WINDOW.
C     INTEGER WORK
C        THE STARTING LOCATION IN THE WORK AREA FOR ARRAY WINDOW.
C     INTEGER XAXIS
C        THE STARTING LOCATION IN THE WORK AREA FOR ARRAY XAXIS.
C     DOUBLE PRECISION Y(N)
C        THE ARRAY CONTAINING THE OBSERVED TIME SERIES.
C     INTEGER YAXIS
C        THE STARTING LOCATION IN THE WORK AREA FOR ARRAY YAXIS.
C     DOUBLE PRECISION YMISS
C        THE MISSING VALUE CODE FOR THE SERIES.
C
C     SET UP NAME ARRAYS
C
      DATA NMSUB(1), NMSUB(2), NMSUB(3), NMSUB(4), NMSUB(5), NMSUB(6) /
     +   'U','F','S','M','S',' '/
C
C     SET UP
C
      OPTION(4) = .TRUE.
      OPTION(3) = .FALSE.
      OPTION(2) = .TRUE.
      OPTION(1) = .FALSE.
C
C     SET MAXIMUM LAG VALUE TO BE USED.
C
      LAGMAX = N - 1
      IF (NW.LE.0) GO TO 20
      LAGMAX = LAGS(1)
      DO 10 I=1,NW
         LAGMAX = MAX(LAGMAX,LAGS(I))
   10 CONTINUE
   20 CONTINUE
      LACOV = LAGMAX + 1
      LNLPPA = LAGMAX + 1
C
C     COMPUTE MINIMUM ALLOWABLE STACK LENGTH
C
      IO = 1
      IF (NPRT.EQ.0) IO = 0
C
      CALL LDSCMP(6, 0, LAGMAX+1+IO*(NF+5), 0, 0, 0, 'D',
     +   2*LAGMAX+2+IO*(2*NF+10), LDSMIN)
C
      LY = N
      LNLPPA = LACOV
      LPCV = NF + 5
      LWORK = LAGMAX+1
C
C     SET SIZE OF WORK AREA.
C     SET THE NUMBER OF OUTSTANDING ALLOCATIONS.
C     SET THE STACK ALLOCATION TYPE.
C
      CALL STKSET(LDSTAK, 4)
      NALL0 = STKST(1)
C
      IFP = 4
C
C     SET STARTING LOCATIONS IN THE WORK AREA FOR VARIOUS ARRAYS.
C
      IF ((LDSMIN.GT.LDSTAK) .OR. (LDSMIN.LE.6)) THEN
         NLPPA = 1
         ACOV = 1
         WORK = 1
         XAXIS = 1
         YAXIS = 1
         ISYM = 1
         ISORT = 1
      ELSE
         NLPPA = STKGET(LACOV,2)
         ACOV = STKGET(LACOV,IFP)
         WORK = STKGET(LWORK,IFP)
         IF (NPRT.NE.0) THEN
            XAXIS = STKGET(LPCV,IFP)
            YAXIS = STKGET(LPCV,IFP)
            ISYM = STKGET(LPCV,2)
            ISORT = ISYM
         ELSE
            XAXIS = WORK
            YAXIS = WORK
            ISYM = WORK
            ISORT = ISYM
         END IF
      END IF
C
C     CALL THE CONTROLLING ROUTINE FOR FOURIER SPECTRUM ROUTINES
C     FOR SERIES WITH MISSING DATA.
C
      CALL UFSDRV(Y, LY, YMISS, RSTAK(ACOV), ISTAK(NLPPA), SPCF, ISPCF,
     +   NF, FMIN, FMAX, FREQ, N, NW, LAGMAX, LAGS, RSTAK(WORK), LACOV,
     +   LWORK, DELTA, ISTAK(ISORT), ISTAK(ISYM), RSTAK(XAXIS),
     +   RSTAK(YAXIS), LPCV, ALPHA, NPRT, PARZEN, NMSUB, LDSMIN, LDSTAK,
     +   OPTION, LNLPPA, LY)
C
      CALL STKCLR(NALL0)
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
     +   '       CALL UFSMS (Y, YMISS, N,'/
     +   '      +            NW, LAGS, NF, FMIN, FMAX, NPRT,'/
     +   '      +            SPCF, ISPCF, FREQ, LDSTAK)')
      END
