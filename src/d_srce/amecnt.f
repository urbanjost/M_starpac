*AMECNT
      SUBROUTINE AMECNT(Y, WT, LWT, XM, N, M, IXM, MDL, NLDRV, APRXDV,
     +   DRV, PAR, NPAR, RES, IFIXED, LIFIXD, STP, LSTP, MIT, STOPSS,
     +   STOPP, SCALE, LSCALE, DELTA, IVAPRX, RSD, PV, LPV, SDPV,
     +   LSDPV, SDRES, LSDRES, VCV, IVCV, WEIGHT, SAVE, NNZW, NPARE,
     +   NLHDR, PAGE, WIDE, IPTOUT, NDIGIT, HLFRPT, NRESTS)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE CONTROLLING SUBROUTINE FOR NONLINEAR LEAST
C     SQUARES REGRESSION.
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  OCTOBER 3, 1983
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      DOUBLE PRECISION
     +   DELTA,RSD,STOPP,STOPSS
      INTEGER
     +   IVAPRX,IVCV,IXM,LIFIXD,LPV,LSCALE,LSDPV,LSDRES,LSTP,LWT,M,
     +   MIT,N,NDIGIT,NNZW,NPAR,NPARE,NRESTS
      LOGICAL
     +   APRXDV,HLFRPT,PAGE,SAVE,WEIGHT,WIDE
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   PAR(*),PV(*),RES(*),SCALE(*),SDPV(*),SDRES(*),STP(*),
     +   VCV(IVCV,*),WT(*),XM(IXM,*),Y(*)
      INTEGER
     +   IFIXED(*),IPTOUT(5)
C
C  SUBROUTINE ARGUMENTS
      EXTERNAL DRV,MDL,NLDRV,NLHDR
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
     +   D,IFIXD,IFP,IIWORK,IRWORK,IWORK,LVCVL,NALL0,PARE,PVI,
     +   RESTS,RWORK,SDPVI,SDRESI,VCVL
C
C  LOCAL ARRAYS
      DOUBLE PRECISION
     +   RSTAK(12)
      INTEGER
     +   ISTAK(12)
C
C  EXTERNAL FUNCTIONS
      INTEGER
     +   STKGET,STKST
      EXTERNAL STKGET,STKST
C
C  EXTERNAL SUBROUTINES
      EXTERNAL AMEMN,CPYASF,CPYVII,DCOPY,SETIV,STKCLR
C
C  COMMON BLOCKS
      COMMON /CSTAK/DSTAK
      COMMON /ERRCHK/IERR
C
C  EQUIVALENCES
      EQUIVALENCE (DSTAK(1),ISTAK(1))
      EQUIVALENCE (DSTAK(1),RSTAK(1))
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     LOGICAL APRXDV
C        THE VARIABLE USED TO INDICATE WHETHER NUMERICAL
C        APPROXIMATIONS TO THE DERIVATIVE WERE USED (TRUE) OR NOT
C        (FALSE).
C     INTEGER D
C        THE STARTING LOCATION IN RSTAK/DSTAK OF
C        THE ARRAY IN WHICH THE NUMERICAL DERIVATIVES WITH RESPECT TO
C        EACH PARAMETER ARE STORED.
C     DOUBLE PRECISION DELTA
C        THE MAXIMUM CHANGE ALLOWED IN THE MODEL PARAMETERS AT THE
C        FIRST ITERATION.
C     EXTERNAL DRV
C        THE NAME OF THE USER SUPPLIED SUBROUTINE WHICH COMPUTES THE
C        DERIVATIVE (JACOBIAN) MATRIX OF THE MODEL.
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     LOGICAL HLFRPT
C        THE VARIABLE WHICH INDICATES WHETHER THE DERIVATIVE
C        CHECKING ROUTINE HAS ALREADY PRINTED PART OF THE
C        INITIAL SUMMARY (TRUE) OR NOT (FALSE).
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .GE. 1, ERRORS WERE DETECTED.
C     INTEGER IFIXD
C        THE STARTING LOCATION IN ISTAK OF
C        THE INDICATOR VALUES USED TO DESIGNATE WHETHER THE
C        PARAMETERS ARE TO BE OPTIMIZED OR ARE TO BE HELD FIXED.
C        IF IFIXED(I).NE.0, THEN PAR(I) WILL BE HELD FIXED.
C        IF IFIXED(I).EQ.0, THEN PAR(I) WILL BE OPTIMIZED.
C     INTEGER IFIXED(LIFIXD)
C        THE INDICATOR VALUES USED TO DESIGNATE WHETHER THE
C        PARAMETERS ARE TO BE OPTIMIZED OR ARE TO BE HELD FIXED.  IF
C        IFIXED(I).NE.0, THEN PAR(I) WILL BE OPTIMIZED.  IF IFIXED(I).EQ
C        THEN PAR(I) WILL BE HELD FIXED.
C     INTEGER IFP
C        AN INDICATOR FOR STACK ALLOCATION TYPE, WHERE IFP=3 INDICATES
C        REAL AND IFP=4 INDICATES DOUBLE PRECISION.
C     INTEGER IIWORK
C        THE DIMENSION OF THE INTEGER WORK VECTOR IWORK.
C     INTEGER IPTOUT(5)
C        THE VARIABLE USED TO CONTROL PRINTED OUTPUT FOR EACH SECTION.
C     INTEGER IRWORK
C        THE DIMENSION OF THE DOUBLE PRECISION WORK VECTOR RWORK.
C     INTEGER ISTAK(12)
C        THE INTEGER VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER IVAPRX
C        AN INDICATOR VALUE USED TO DESIGNATE WHICH OPTION IS TO BE USED
C        TO COMPUTE THE VARIANCE COVARIANCE MATRIX (VCV), WHERE
C        IVAPRX LE 0 INDICATES THE THE DEFAULT OPTION WILL BE USED
C        IVAPRX EQ 1 INDICATES THE VCV IS TO BE COMPUTED BY
C                       INVERSE(TRANSPOSE(JACOBIAN)*JACOBIAN)
C                    USING BOTH THE MODEL SUBROUTINE THE USER SUPPLIED
C                    DERIVATIVE SUBROUTINE WHEN IT IS AVAILABLE
C        IVAPRX EQ 2 INDICATES THE VCV IS TO BE COMPUTED BY
C                       INVERSE(HESSIAN)
C                    USING BOTH THE MODEL SUBROUTINE THE USER SUPPLIED
C                    DERIVATIVE SUBROUTINE WHEN IT IS AVAILABLE
C        IVAPRX EQ 3 INDICATES THE VCV IS TO BE COMPUTED BY
C                       INVERSE(HESSIAN)*TRANSPOSE(JACOBIAN)*JACOBIAN
C                          *INVERSE(HESSIAN)
C                    USING BOTH THE MODEL SUBROUTINE THE USER SUPPLIED
C                    DERIVATIVE SUBROUTINE WHEN IT IS AVAILABLE
C        IVAPRX EQ 4 INDICATES THE VCV IS TO BE COMPUTED BY
C                       INVERSE(TRANSPOSE(JACOBIAN)*JACOBIAN)
C                    USING ONLY THE MODEL SUBROUTINE
C        IVAPRX EQ 5 INDICATES THE VCV IS TO BE COMPUTED BY
C                       INVERSE(HESSIAN)
C                    USING ONLY THE MODEL SUBROUTINE
C        IVAPRX EQ 6 INDICATES THE VCV IS TO BE COMPUTED BY
C                       INVERSE(HESSIAN)*TRANSPOSE(JACOBIAN)*JACOBIAN
C                          *INVERSE(HESSIAN)
C                    USING ONLY THE MODEL SUBROUTINE
C        IVAPRX GE 7 INDICATES THE DEFAULT OPTION WILL BE USED
C     INTEGER IVCV
C        THE FIRST DIMENSION OF THE VARIANCE COVARIANCE MATRIX VCV.
C     INTEGER IWORK
C        THE STARTING LOCATION IN ISTAK OF
C        THE INTEGER WORK SPACE VECTOR USED BY THE NL2 SUBROUTINES.
C     INTEGER IXM
C        THE FIRST DIMENSION OF THE INDEPENDENT VARIABLE ARRAY.
C     INTEGER LIFIXD
C        THE ACTUAL LENGTH OF THE VECTOR IFIXED.
C     INTEGER LPV
C        THE ACTUAL LENGTH OF THE VECTOR PV.
C     INTEGER LSCALE
C        THE ACTUAL LENGTH OF THE VECTOR SCALE.
C     INTEGER LSDPV
C        THE ACTUAL LENGTH OF THE VECTOR SDPV.
C     INTEGER LSDRES
C        THE ACTUAL LENGTH OF THE VECTOR SDRES.
C     INTEGER LSTP
C        THE ACTUAL LENGTH OF THE VECTOR STP.
C     INTEGER LVCVL
C        THE LENGTH OF THE VECTOR CONTAINING
C        THE LOWER HALF OF THE VCV MATRIX, STORED ROW WISE.
C     INTEGER LWT
C        THE ACTUAL LENGTH OF THE VECTOR WT.
C     INTEGER M
C        THE NUMBER OF INDEPENDENT VARIABLES.
C     INTEGER MIT
C        THE MAXIMUM NUMBER OF ITERATIONS ALLOWED.
C     EXTERNAL MDL
C        THE NAME OF THE USER SUPPLIED SUBROUTINE WHICH COMPUTES THE
C        PREDICTED VALUES BASED ON THE CURRENT PARAMETER ESTIMATE.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS.
C     INTEGER NALL0
C        NUMBER OF ALLOCATIONS ON ENTRY.
C     INTEGER NDIGIT
C        THE NUMBER OF DIGITS IN THE PRINT CONTROL VALUE.
C     EXTERNAL NLDRV
C        THE NAME OF THE ROUTINE WHICH CALCULATES THE DERIVATIVES.
C     EXTERNAL NLHDR
C        THE NAME OF THE ROUTINE WHICH PRODUCES THE HEADING.
C     INTEGER NNZW
C        THE NUMBER OF NON ZERO WEIGHTS.
C     INTEGER NPAR
C        THE NUMBER OF PARAMETERS IN THE MODEL.
C     INTEGER NPARE
C        THE NUMBER OF PARAMETERS TO BE ESTIMATED.
C     INTEGER NRESTS
C        THE MAXIMUM NUMBER OF RESIDUALS TO BE COMPUTED.
C     LOGICAL PAGE
C        THE VARIABLE USED TO INDICATE WHETHER A GIVEN SECTION OF
C        THE OUTPUT IS TO BEGIN ON A NEW PAGE (TRUE) OR NOT (FALSE).
C     DOUBLE PRECISION PAR(NPAR)
C        THE CURRENT ESTIMATES OF THE PARAMETERS.
C     INTEGER PARE
C        THE STARTING LOCATION IN RSTAK/DSTAK OF
C        THE CURRENT ESTIMATES OF THE PARAMETERS, BUT ONLY
C        THOSE TO BE OPTIMIZED (NOT THOSE WHOSE VALUES ARE FIXED).
C     DOUBLE PRECISION PV(LPV)
C        THE PREDICTED VALUE BASED ON THE CURRENT PARAMETER ESTIMATES
C     INTEGER PVI
C        THE STARTING LOCATION IN RSTAK/DSTAK OF
C        THE PREDICTED VALUES.
C     DOUBLE PRECISION RES(N)
C        THE RESIDUALS FROM THE FIT.
C     INTEGER RESTS
C        THE STARTING LOCATION IN RSTAK/DSTAK OF
C        THE RESIDUALS FROM THE ARIMA MODEL.
C     DOUBLE PRECISION RSD
C        THE VALUE OF THE RESIDUAL STANDARD DEVIATION AT THE SOLUTION.
C     DOUBLE PRECISION RSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER RWORK
C        THE STARTING LOCATION IN RSTAK/DSTAK OF
C        THE DOUBLE PRECISION WORK VECTOR USED BY THE NL2 SUBROUTINES.
C     LOGICAL SAVE
C        THE VARIABLE USED TO INDICATE WHETHER ANY RESULTS OTHER THAN
C        THE RESIDUALS AND PARAMETERS ARE TO BE SAVED (TRUE) OR NOT
C        (FALSE).
C     DOUBLE PRECISION SCALE(LSCALE)
C        A VALUE TO INDICATE USE OF THE DEFAULT VALUES OF
C        THE TYPICAL SIZE OF THE PARAMETERS.
C     DOUBLE PRECISION SDPV(LSDPV)
C        THE STANDARD DEVIATIONS OF THE PREDICTED VALUES.
C     INTEGER SDPVI
C        THE STARTING LOCATION IN RWORK OF
C        THE STANDARD DEVIATIONS OF THE PREDICTED VALUES.
C     DOUBLE PRECISION SDRES(LSDRES)
C        THE STANDARDIZED RESIDUALS.
C     INTEGER SDRESI
C        THE STARTING LOCATION IN RWORK OF THE
C        THE STANDARDIZED RESIDUALS.
C     DOUBLE PRECISION STOPP
C        THE STOPPING CRITERION FOR THE TEST BASED ON THE MAXIMUM SCALED
C        RELATIVE CHANGE IN THE ELEMENTS OF THE MODEL PARAMETER VECTOR
C     DOUBLE PRECISION STOPSS
C        THE STOPPING CRITERION FOR THE TEST BASED ON THE RATIO OF THE
C        PREDICTED DECREASE IN THE RESIDUAL SUM OF SQUARES (COMPUTED
C        BY STARPAC) TO THE CURRENT RESIDUAL SUM OF SQUARES ESTIMATE.
C     DOUBLE PRECISION STP(LSTP)
C        THE STEP SIZE ARRAY.
C     DOUBLE PRECISION VCV(IVCV,NPAR)
C        THE VARIANCE-COVARIANCE MATRIX.
C     INTEGER VCVL
C        THE STARTING LOCATION IN RWORK OF
C        THE VARIANCE-COVARIANCE MATRIX.
C     LOGICAL WEIGHT
C        THE VARIABLE USED TO INDICATE WHETHER WEIGHTED ANALYSIS IS TO
C        BE PERFORMED (TRUE) OR NOT (FALSE).
C     LOGICAL WIDE
C        THE VARIABLE USED TO INDICATE WHETHER THE HEADING SHOULD
C        BE FULL WIDTH (TRUE) OR NOT (FALSE).
C     DOUBLE PRECISION WT(LWT)
C        THE USER SUPPLIED WEIGHTS.
C     DOUBLE PRECISION XM(IXM,M)
C        THE ARRAY IN WHICH ONE ROW OF THE INDEPENDENT VARIABLE ARRAY
C        IS STORED.
C     DOUBLE PRECISION Y(N)
C        THE ARRAY OF THE DEPENDENT VARIABLE.
C
C     SET VARIOUS PROGRAM VALUES
C
      NALL0 = STKST(1)
C
      IFP = 4
C
      IERR = 0
C
C     SUBDIVIDE WORK AREA FOR LEAST SQUARES ANALYSIS
C
      IIWORK = NPARE + 60
      IRWORK = 94 + 2*NRESTS + NPARE*(3*NPARE+33)/2
C
      IFIXD = STKGET(NPAR,2)
      IWORK = STKGET(IIWORK,2)
C
      D = STKGET(NRESTS*NPAR,IFP)
      PARE = STKGET(NPARE,IFP)
      RESTS = STKGET(NRESTS,IFP)
      PVI = RESTS
      RWORK = STKGET(IRWORK,IFP)
C
      IF (IERR.EQ.1) RETURN
C
C     SET VALUES FOR IFIXD
C
      IF (IFIXED(1).GE.0) CALL CPYVII(NPAR, IFIXED, 1, ISTAK(IFIXD), 1)
      IF (IFIXED(1).LT.0) CALL SETIV(ISTAK(IFIXD), NPAR, 0)
C
      CALL AMEMN(Y, WEIGHT, NNZW, WT, LWT, XM, N, M, IXM, NRESTS,
     +   APRXDV, ISTAK(IFIXD), PAR, RSTAK(PARE), NPAR, RES, PAGE,
     +   WIDE, HLFRPT, STP, LSTP, MIT, STOPSS, STOPP, SCALE, LSCALE,
     +   DELTA, IVAPRX, IPTOUT, NDIGIT, RSD, RSTAK(RESTS), SDPVI,
     +   SDRESI, VCVL, LVCVL, RSTAK(D), ISTAK(IWORK), IIWORK,
     +   RSTAK(RWORK), IRWORK, NLHDR, NPARE, RSTAK(PVI))
C
      IF (.NOT.SAVE) GO TO 10
C
      SDPVI = RWORK + SDPVI - 1
      SDRESI = RWORK + SDRESI - 1
      VCVL = RWORK + VCVL - 1
C
      CALL DCOPY(N, RSTAK(PVI), 1, PV, 1)
      CALL DCOPY(N, RSTAK(SDPVI), 1, SDPV, 1)
      CALL DCOPY(N, RSTAK(SDRESI), 1, SDRES, 1)
      CALL CPYASF(NPARE, RSTAK(VCVL), LVCVL, VCV, IVCV)
C
   10 CALL STKCLR(NALL0)
C
      RETURN
C
      END
