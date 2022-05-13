*LLCNT
      SUBROUTINE LLCNT(Y, WT, LWT, XM, N, M, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, LPAR, RSD, PV, LPV, SDPV, LSDPV, SDRES, LSDRES, VCV,
     +   IVCV, LLHDR, IFIT, NMSUB, WEIGHT, SAVE)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE CONTROLLING SUBROUTINE FOR LINEAR LEAST
C     SQUARES.
C
C     WRITTEN BY - JANET R. DONALDSON
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  DECEMBER 29, 1982
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      REAL
     +   RSD
      INTEGER
     +   IFIT,IVCV,IXM,LDSTAK,LPAR,LPV,LSDPV,LSDRES,LWT,M,N,NPAR,
     +   NPRT
      LOGICAL
     +   SAVE,WEIGHT
C
C  ARRAY ARGUMENTS
      REAL
     +   PAR(*),PV(*),RES(*),SDPV(*),SDRES(*),VCV(*),WT(*),XM(*),Y(*)
      CHARACTER
     +   NMSUB(6)*1
C
C  SUBROUTINE ARGUMENTS
      EXTERNAL LLHDR
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  ARRAYS IN COMMON
      DOUBLE PRECISION DSTAK(12)
C
C  LOCAL SCALARS
      INTEGER
     +   ACC,C,FC,IFP,NALL0,NDIGIT,NNZW,PAR1,PARI,PVI,RED,RESI,
     +   RSDI,SDPVI,SDRESI,T,VCVI,WTI,WY,XMW
      LOGICAL
     +   PAGE,WIDE
C
C  LOCAL ARRAYS
      REAL
     +   RSTAK(12)
      INTEGER
     +   IPTOUT(4),ISTAK(12)
C
C  EXTERNAL FUNCTIONS
      INTEGER
     +   STKGET,STKST
      EXTERNAL STKGET,STKST
C
C  EXTERNAL SUBROUTINES
      EXTERNAL CPYMSS,LLER,LLSMN,PRTCNT,SCOPY,SETRV,STKCLR,STKSET
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
C     INTEGER ACC
C        THE STARTING LOCATION IN THE WORK AREA FOR
C        THE NUMBER OF ACCURATE DIGITS.
C     INTEGER C
C        *
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER FC
C        THE STARTING LOCATION IN THE WORK AREA FOR
C        THE ORTHONORMALIZATION MATRIX.
C     INTEGER IERR
C        THE INTEGER VALUE DESIGNATING WHETHER ANY ERRORS WERE
C        DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .EQ. 1, ERRORS WERE DETECTED.
C     INTEGER IFIT
C        THE INDICATOR VALUE DESIGNATING WHETHER THE FIT IS OF A
C        GENERAL MODEL (IFIT=3) OR A POLYNOMIAL MODEL (IFIT=1).
C     INTEGER IFP
C        AN INDICATOR FOR STACK ALLOCATION TYPE, WHERE IFP=3 INDICATES
C        SINGLE PRECISION AND IFP=4 INDICATES DOUBLE PRECISION.
C     INTEGER IPTOUT(4)
C        THE VARIABLE USED TO CONTROL PRINTED OUTPUT FOR EACH SECTION.
C     INTEGER ISTAK(12)
C        THE INTEGER VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER IVCV
C        THE FIRST DIMENSION OF THE MATRIX VCV.
C     INTEGER IXM
C        THE FIRST DIMENSION OF THE MATRIX XM.
C     INTEGER LDSTAK
C        THE LENGTH OF THE VECTOR DSTAK IN COMMON CSTAK.
C     EXTERNAL LLHDR
C        THE NAME OF THE ROUTINE THAT PRODUCED THE HEADING.
C     INTEGER LPAR
C        THE ACTUAL LENGTH OF THE VECTOR P.
C     INTEGER LPV
C        THE ACTUAL LENGTH OF THE VECTOR PV.
C     INTEGER LSDPV
C        THE ACTUAL LENGTH OF THE VECTOR SDPV.
C     INTEGER LSDRES
C        THE ACTUAL LENGTH OF THE VECTOR SDRES.
C     INTEGER LWT
C        THE ACTUAL LENGTH OF THE VECTOR WT.
C     INTEGER M
C        THE NUMBER OF INDEPENDENT VARIABLES.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS.
C     INTEGER NALL0
C        NUMBER OF ALLOCATIONS ON ENTRY.
C     INTEGER NDIGIT
C        THE NUMBER OF DIGITS IN THE PRINT CONTROL VALUE.
C     CHARACTER*1 NMSUB(6)
C        THE NAME OF THE SUBROUTINE CALLING THE ERROR CHECKING
C        SUBROUTINES.
C     INTEGER NNZW
C        THE NUMBER OF NON ZERO WEIGHTS.
C     INTEGER NPAR
C        THE NUMBER OF PARAMETERS.
C     INTEGER NPRT
C        THE INDICATOR VARIABLE USED TO DESIGNATE THE AMOUNT OF
C        PRINTED OUTPUT.
C     LOGICAL PAGE
C        THE VARIABLE USED TO INDICATE WHETHER A GIVEN SECTION OF
C        THE OUTPUT IS TO BEGIN ON A NEW PAGE (TRUE) OR NOT (FALSE).
C     REAL PAR(LPAR)
C        THE PARAMETERS TO BE ESTIMATED.
C     INTEGER PARI
C        THE STARTING LOCATION IN THE WORK AREA OF
C        THE PARAMETERS TO BE ESTIMATED.
C     INTEGER PAR1
C        THE STARTING LOCATION IN THE WORK AREA FOR
C        THE PARAMETERS TO BE ESTIMATED OMMITTING THE LAST
C        INDEPENDENT VARIABLE.
C     REAL PV(LPV)
C        THE PREDICTED VALUES.
C     INTEGER PVI
C        THE STARTING LOCATION IN THE WORK AREA FOR
C        THE PREDICTED VALUES.
C     INTEGER RED
C        THE STARTING LOCATION IN THE WORK AREA FOR
C        THE REDUCTION TO THE SUM OF SQUARES DUE TO EACH PARAMETER.
C     REAL RES(N)
C        THE RESIDUALS.
C     INTEGER RESI
C        THE STARTING LOCATION IN THE WORK AREA FOR
C        THE RESIDUALS.
C     REAL RSD
C        THE RESIDUAL STANDARD DEVIATION.
C     INTEGER RSDI
C        THE STARTING LOCATION IN THE WORK AREA FOR
C        THE RESIDUAL STANDARD DEVIATION.
C     REAL RSTAK(12)
C        THE REAL VERSION OF THE /CSTAK/ WORK AREA.
C     LOGICAL SAVE
C        THE VARIABLE USED TO INDICATE WHETHER ANY RESULTS OTHER THAN
C        THE RESIDUALS ARE TO VE SAVED (TRUE) OR NOT (FALSE).
C     REAL SDPV(LSDPV)
C        THE STANDARD DEVIATIONS OF THE PREDICTED VALUES.
C     INTEGER SDPVI
C        THE STARTING LOCATION IN THE WORK AREA FOR
C        THE STANDARD DEVIATIONS OF THE PREDICTED VALUES.
C     REAL SDRES(LSDRES)
C        THE STANDARDIZED RESIDUALS.
C     INTEGER SDRESI
C        THE STARTING LOCATION IN THE WORK AREA FOR
C        THE STANDARDIZED RESIDUALS.
C     INTEGER T
C        THE STARTING LOCATION IN THE WORK AREA FOR
C        THE TRIANGULAR MATRIX FROM THE DECOMPOSITION.
C     REAL VCV(IVCV,NPAR)
C        THE VARIANCE COVARIANCE MATRIX.
C     INTEGER VCVI
C        THE STARTING LOCATION IN THE WORK AREA FOR
C        THE VARIANCE COVARIANCE MATRIX.
C     LOGICAL WEIGHT
C        THE VARIABLE USED TO INDICATE WHETHER WEIGHTED ANALYSIS IS TO
C        BE PERFORMED (TRUE) OR NOT (FALSE).
C     LOGICAL WIDE
C        THE VARIABLE USED TO INDICATE WHETHER THE HEADING SHOULD
C        BE FULL WIDTH (TRUE) OR NOT (FALSE).
C     REAL WT(LWT)
C        THE WEIGHTS.
C     INTEGER WTI
C        THE STARTING LOCATION IN THE WORK AREA FOR
C        THE WEIGHTS.
C     INTEGER WY
C        THE STARTING LOCATION IN THE WORK AREA FOR
C        THE VECTOR CONTAINING SQRT(WT)*Y.
C     REAL XM(IXM,M)
C        THE INDEPENDENT VARIABLE.
C     INTEGER XMW
C        THE STARTING LOCATION IN THE WORK AREA FOR
C        THE MATRIX CONTAINING XM * SQRT(WT).
C     REAL Y(N)
C        THE DEPENDENT VARIABLE.
C
C
      WIDE = .TRUE.
      PAGE = .FALSE.
      NDIGIT = 4
C
      IFP = 3
C
C     SET PRINT CONTROL VALUES
C
      CALL PRTCNT(NPRT, NDIGIT, IPTOUT)
C
C     CHECK FOR ERRORS
C
      CALL LLER(NMSUB, IXM, IVCV, N, NPAR, LPAR, LDSTAK, WT, LWT,
     +   WEIGHT, NNZW, IFIT, SAVE)
      IF (IERR.NE.0) RETURN
C
      CALL STKSET(LDSTAK, 4)
      NALL0 = STKST(1)
C
C     SET UP SUBDIVISION OF WORK AREAS
C
      WTI = STKGET(N,IFP)
      RESI = STKGET(N,IFP)
      RSDI = STKGET(1,IFP)
      PARI = STKGET(NPAR,IFP)
      PVI = STKGET(N,IFP)
      SDPVI = STKGET(N,IFP)
      SDRESI = STKGET(N,IFP)
      VCVI = STKGET(NPAR*NPAR,IFP)
C
      WY = STKGET(N,IFP)
      XMW = STKGET(N*NPAR,IFP)
      RED = STKGET(NPAR,IFP)
      T = STKGET(NPAR*NPAR,IFP)
      PAR1 = STKGET(NPAR,IFP)
      ACC = STKGET(NPAR,IFP)
      C = STKGET(NPAR,IFP)
C
C     EQUIVALENCED LOCATIONS WITHIN SCRAT
C
      FC = XMW
C
C     SET UP WEIGHTS VECTOR
C
      IF (WEIGHT) THEN
         CALL SCOPY(N, WT, 1, RSTAK(WTI), 1)
      ELSE
         CALL SETRV(RSTAK(WTI), N, 1.0E0)
      END IF
C
      CALL LLSMN(Y, XM, RSTAK(WTI), N, M, NPAR, IXM, RSTAK(RESI),
     +   RSTAK(PARI), NNZW, RSTAK(RSDI), RSTAK(PVI), RSTAK(SDPVI),
     +   RSTAK(SDRESI), IPTOUT, RSTAK(WY), RSTAK(XMW), RSTAK(VCVI),
     +   RSTAK(FC), RSTAK(RED), RSTAK(T), RSTAK(PAR1), RSTAK(ACC), IFIT,
     +   WEIGHT, RSTAK(C), LLHDR, PAGE, WIDE)
C
      CALL SCOPY(N, RSTAK(RESI), 1, RES, 1)
C
      IF (SAVE) THEN
         RSD = RSTAK(RSDI)
         CALL SCOPY(NPAR, RSTAK(PARI), 1, PAR, 1)
         CALL SCOPY(N, RSTAK(PVI), 1, PV, 1)
         CALL SCOPY(N, RSTAK(SDPVI), 1, SDPV, 1)
         CALL SCOPY(N, RSTAK(SDRESI), 1, SDRES, 1)
         CALL CPYMSS(NPAR, NPAR, RSTAK(VCVI), NPAR, VCV, IVCV)
      END IF
      CALL STKCLR(NALL0)
C
      IF (IERR.EQ.3) IERR = 2
      IF (IERR.EQ.4) IERR = 3
C
      RETURN
C
      END
