*AMEDRV
      SUBROUTINE AMEDRV(Y, N, MSPEC, NFAC, PAR, NPAR,
     +   RES, LDSTAK, IFIXED, LIFIXD, STP, LSTP, MIT, STOPSS, STOPP,
     +   SCALE, LSCALE, DELTA, IVAPRX, NPRT, RSD, PV, LPV, SDPV, LSDPV,
     +   SDRES, LSDRES, VCV, IVCV, NMSUB, SAVE, NPARE)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE CONTROLLING SUBROUTINE FOR NONLINEAR LEAST
C     SQUARES REGRESSION USING NUMERICALLY APPROXIMATED DERIVATIVES.
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  DECEMBER 2, 1985
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      DOUBLE PRECISION
     +   DELTA,RSD,STOPP,STOPSS
      INTEGER
     +   IVAPRX,IVCV,LDSTAK,LIFIXD,LPV,LSCALE,LSDPV,LSDRES,LSTP,
     +   MIT,N,NFAC,NPAR,NPARE,NPRT
      LOGICAL
     +   SAVE
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   PAR(*),PV(*),RES(*),SCALE(*),SDPV(*),SDRES(*),STP(*),VCV(*),
     +   Y(*)
      INTEGER
     +   IFIXED(*),MSPEC(4,*)
      CHARACTER
     +   NMSUB(6)*1
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR,IFLAG,MBO,MBOL,MSPECT,NFACT,NPARAR,NPARDF,NPARMA,
     +   NRESTS,PARAR,PARDF,PARMA,T,TEMP
C
C  ARRAYS IN COMMON
      DOUBLE PRECISION DSTAK(12)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   EXMPT
      INTEGER
     +   IFP,IS,ISUBHD,IXM,LDSMIN,LWT,M,NALL0,NDIGIT,NETA,NNZW,STPT
      LOGICAL
     +   APRXDV,HLFRPT,PAGE,PRTFXD,WEIGHT,WIDE
C
C  LOCAL ARRAYS
      DOUBLE PRECISION
     +   RSTAK(12),WT(1)
      INTEGER
     +   IPTOUT(5),ISTAK(12)
C
C  EXTERNAL FUNCTIONS
      INTEGER
     +   ICNTI,STKGET,STKST
      EXTERNAL ICNTI,STKGET,STKST
C
C  EXTERNAL SUBROUTINES
      EXTERNAL AMECNT,AMEER,AMEHDR,AMESTP,BACKOP,CPYVII,
     +   DCOEF,DRV,LDSCMP,MDLTS1,MDLTS3,NLDRVN,PRTCNT,DCOPY,
     +   STKCLR,STKSET,STPAMO
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MAX
C
C  COMMON BLOCKS
      COMMON /CSTAK/DSTAK
      COMMON /ERRCHK/IERR
      COMMON /MDLTSC/MSPECT,NFACT,PARDF,NPARDF,PARAR,NPARAR,PARMA,
     +   NPARMA,MBO,MBOL,T,TEMP,NRESTS,IFLAG
C
C  EQUIVALENCES
      EQUIVALENCE (DSTAK(1),ISTAK(1))
      EQUIVALENCE (DSTAK(1),RSTAK(1))
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     EXTERNAL AMEHDR
C        THE ROUTINE USED TO PRINT THE HEADING
C     LOGICAL APRXDV
C        THE VARIABLE USED TO INDICATE WHETHER NUMERICAL
C        APPROXIMATIONS TO THE DERIVATIVE WERE USED (TRUE) OR NOT
C        (FALSE).
C     DOUBLE PRECISION DELTA
C        THE MAXIMUM CHANGE ALLOWED IN THE MODEL PARAMETERS AT THE
C        FIRST ITERATION.
C     EXTERNAL DRV
C        THE NAME OF THE USER SUPPLIED SUBROUTINE WHICH COMPUTES THE
C        DERIVATIVE (JACOBIAN) MATRIX OF THE MODEL.
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     DOUBLE PRECISION EXMPT
C        THE PROPORTION OF OBSERVATIONS FOR WHICH THE COMPUTED
C        NUMERICAL DERIVATIVES WRT A GIVEN PARAMETER ARE EXEMPTED
C        FROM MEETING THE DERIVATIVE ACCEPTANCE CRITERIA.
C     LOGICAL HLFRPT
C        THE VARIABLE WHICH INDICATES WHETHER THE DERIVATIVE
C        CHECKING ROUTINE HAS ALREADY PRINTED PART OF THE
C        INITIAL SUMMARY (TRUE) OR NOT (FALSE).
C     INTEGER IERR
C        THE VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .GE. 1, ERRORS WERE DETECTED.
C     INTEGER IFIXED(LIFIXD)
C        THE INDICATOR VALUES USED TO DESIGNATE WHETHER THE
C        PARAMETERS ARE TO BE OPTIMIZED OR ARE TO BE HELD FIXED.
C        IF IFIXED(I).NE.0, THEN PAR(I) WILL BE HELD FIXED.
C        IF IFIXED(I).EQ.0, THEN PAR(I) WILL BE OPTIMIZED.
C     INTEGER IFLAG
C        ...
C     INTEGER IFP
C        AN INDICATOR FOR THE PRECISION OF THE STACK ALLOCATION TYPE,
C        WHERE IFP=3 INDICATES SINGLE AND IFP=4 INDICATES DOUBLE.
C     INTEGER IPTOUT(5)
C        THE VARIABLE USED TO CONTROL PRINTED OUTPUT FOR EACH SECTION.
C     INTEGER IS
C        A VALUE USED TO DETERMINE THE AMOUNT OF WORK SPACE NEEDED
C        BASED ON WHETHER STEP SIZES ARE INPUT OR ARE TO BE CALCULATED.
C     INTEGER ISTAK(12)
C        THE INTEGER VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER ISUBHD
C        AN INDICATOR VALUE SPECIFYING SUBHEADINGS TO BE PRINTED.
C     INTEGER IVAPRX
C        AN INDICATOR VALUE USED TO DESIGNATE WHICH OPTION IS TO BE USED
C        TO COMPUTE THE VARIANCE COVARIANCE MATRIX (VCV), WHERE FOR
C        IVAPRX LE 0, VCV = THE DEFAULT OPTION
C        IVAPRX EQ 1, VCV = INVERSE(TRANSPOSE(J)*J)
C                     USING BOTH THE MODEL SUBROUTINE THE USER SUPPLIED
C                     DERIVATIVE SUBROUTINE WHEN IT IS AVAILABLE
C        IVAPRX EQ 2, VCV = INVERSE(H)
C                     USING BOTH THE MODEL SUBROUTINE THE USER SUPPLIED
C                     DERIVATIVE SUBROUTINE WHEN IT IS AVAILABLE
C        IVAPRX EQ 3, VCV = INVERSE(H)*TRANSPOSE(J)*JACOBIAN*INVERSE(H)
C                     USING BOTH THE MODEL SUBROUTINE THE USER SUPPLIED
C                     DERIVATIVE SUBROUTINE WHEN IT IS AVAILABLE
C        IVAPRX EQ 4, VCV = INVERSE(TRANSPOSE(J)*J)
C                     USING ONLY THE MODEL SUBROUTINE
C        IVAPRX EQ 5, VCV = INVERSE(H)
C                     USING ONLY THE MODEL SUBROUTINE
C        IVAPRX EQ 6, VCV = INVERSE(H)*TRANSPOSE(J)*JACOBIAN*INVERSE(H)
C                     USING ONLY THE MODEL SUBROUTINE
C        IVAPRX GE 7, VCV = THE DEFAULT OPTION
C        WITH J REPRESENTING THE JACOBIAN AND H THE HESSIAN.
C     INTEGER IVCV
C        THE FIRST DIMENSION OF MATRIX VCV.
C     INTEGER IXM
C        THE FIRST DIMENSION OF MATRIX XM.
C     INTEGER LDSMIN
C        THE MINIMUM LENGTH ALLOWED FOR THE ARRAY DSTAK.
C     INTEGER LDSTAK
C        THE LENGTH OF THE ARRAY DSTAK.
C     INTEGER LIFIXD
C        THE DIMENSION OF VECTOR IFIXED.
C     INTEGER LPV
C        THE DIMENSION OF VECTOR PV.
C     INTEGER LSCALE
C        THE DIMENSION OF VECTOR SCALE.
C     INTEGER LSDPV
C        THE DIMENSION OF VECTOR SDPV.
C     INTEGER LSDRES
C        THE DIMENSION OF VECTOR SDRES.
C     INTEGER LSTP
C        THE DIMENSION OF VECTOR STP.
C     INTEGER LWT
C        THE DIMENSION OF VECTOR WT.
C     INTEGER M
C        THE NUMBER OF INDEPENDENT VARIABLES.
C     INTEGER MBO
C        THE MAXIMUM BACK ORDER OPERATOR.
C     INTEGER MBOL
C        THE MAXIMUM BACK ORDER ON THE LEFT
C     EXTERNAL MDLTS1
C        THE STARPAC FORMAT SUBROUTINE FOR COMPUTING THE ARIMA MODEL
C        PREDICTED VALUES.
C     EXTERNAL MDLTS3
C        THE STARPAC FORMAT SUBROUTINE FOR COMPUTING THE ARIMA MODEL
C        RESIDUALS.
C     INTEGER MIT
C        THE MAXIMUM NUMBER OF ITERATIONS ALLOWED.
C     INTEGER MSPEC(4,NFAC)
C        THE ARRAY CONTAINING THE VALUES OF P, D, Q, AND S FOR EACH FACT
C     INTEGER MSPECT
C        THE STARTING LOCATION IN THE WORK SPACE FOR
C        THE ARRAY CONTAINING THE VALUES OF P, D, Q, AND S FOR EACH FACT
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS.
C     INTEGER NALL0
C        NUMBER OF STACK ALLOCATIONS OUTSTANDING.
C     INTEGER NDIGIT
C        THE NUMBER OF DIGITS IN THE PRINT CONTROL VALUE.
C     INTEGER NETA
C        THE NUMBER OF ACCURATE DIGITS IN THE MODEL RESULTS.
C     INTEGER NFAC
C        THE NUMBER OF FACTORS IN THE MODEL
C     INTEGER NFACT
C        THE NUMBER OF FACTORS IN THE MODEL
C     EXTERNAL NLDRVN
C        THE NAME OF THE ROUTINE WHICH CALCULATES THE DERIVATIVES.
C     CHARACTER*1 NMSUB(6)
C        THE NAME OF THE ROUTINE CALLING THE ERROR CHECKING ROUTINE
C     INTEGER NNZW
C        THE NUMBER OF NON ZERO WEIGHTS.
C     INTEGER NPAR
C        THE NUMBER OF PARAMETERS IN THE MODEL.
C     INTEGER NPARAR
C        THE NUMBER OF AUTOREGRESSIVE PARAMETERS
C     INTEGER NPARDF
C        THE ORDER OF THE EXPANDED DIFFERENCE FILTER.
C     INTEGER NPARE
C        THE NUMBER OF PARAMETERS TO BE OPTIMIZED.
C     INTEGER NPARMA
C        THE LENGTH OF THE VECTOR PARMA
C     INTEGER NPRT
C        THE PARAMETER USED TO INDICATE HOW MUCH PRINTED OUTPUT IS
C        TO BE PROVIDED.
C     INTEGER NRESTS
C        THE MAXIMUM NUMBER OF RESIDUALS TO BE COMPUTED.
C     LOGICAL PAGE
C        THE VARIABLE USED TO INDICATE WHETHER A GIVEN SECTION OF
C        THE OUTPUT IS TO BEGIN ON A NEW PAGE (TRUE) OR NOT (FALSE).
C     DOUBLE PRECISION PAR(NPAR)
C        THE CURRENT ESTIMATES OF THE PARAMETERS.
C     INTEGER PARAR
C        THE STARTING LOCATION IN THE WORK ARRAY FOR
C        THE AUTOREGRESSIVE PARAMETERS
C     INTEGER PARDF
C        THE STARTING LOCATION IN THE WORK SPACE FOR
C        THE VECTOR CONTAINING THE DIFFERENCE FILTER PARAMETERS
C     INTEGER PARMA
C        THE STARTING LOCATION IN THE WORK ARRAY FOR
C        THE MOVING AVERAGE PARAMETERS
C     LOGICAL PRTFXD
C        THE INDICATOR VALUE USED TO DESIGNATE WHETHER THE
C        OUTPUT IS TO INCLUDE INFORMATION ON WHETHER THE
C        PARAMETER IS FIXED (TRUE) OR NOT (FALSE).
C     DOUBLE PRECISION PV(LPV)
C        THE PREDICTED VALUE OF THE FIT.
C     DOUBLE PRECISION RES(N)
C        THE RESIDUALS FROM THE FIT.
C     DOUBLE PRECISION RSD
C        THE RESIDUAL STANDARD DEVIATION.
C     DOUBLE PRECISION RSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     LOGICAL SAVE
C        THE VARIABLE USED TO INDICATE WHETHER ANY RESULTS OTHER THAN
C        THE RESIDUALS AND PARAMETERS ARE TO BE SAVED (TRUE) OR NOT
C        (FALSE).
C     DOUBLE PRECISION SCALE(LSCALE)
C        THE TYPICAL SIZE OF THE PARAMETERS.
C     DOUBLE PRECISION SDPV(LSDPV)
C        THE STANDARD DEVIATIONS OF THE PREDICTED VALUES.
C     DOUBLE PRECISION SDRES(LSDRES)
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
C     EXTERNAL STPAMO
C        THE ROUTINE USED TO PRINT THE OUTPUT FROM THE STEP SIZE SELECTI
C        ROUTINES.
C     INTEGER STPT
C        THE STARTING LOCATION IN /CSTAK/ OF VECTOR STPT CONTAINING
C        THE STEP SIZE ARRAY.
C     INTEGER T
C        THE STARTING LOCATION IN THE WORK ARRAY FOR
C        A TEMPORARY WORK VECTOR.
C     INTEGER TEMP
C        THE STARTING LOCATION IN THE WORK ARRAY FOR
C        A TEMPORARY WORK VECTOR
C     DOUBLE PRECISION VCV(IVCV,NPAR)
C        THE VARIANCE-COVARIANCE MATRIX.
C     LOGICAL WEIGHT
C        THE VARIABLE USED TO INDICATE WHETHER WEIGHTED ANALYSIS IS TO
C        BE PERFORMED (TRUE) OR NOT (FALSE).
C     LOGICAL WIDE
C        THE VARIABLE USED TO INDICATE WHETHER THE HEADING SHOULD
C        BE FULL WIDTH (TRUE) OR NOT (FALSE).
C     DOUBLE PRECISION WT(1)
C        THE USER SUPPLIED WEIGHTS, UNUSED WHEN WEIGHT = FALSE.
C     DOUBLE PRECISION Y(N)
C        THE DEPENDENT VARIABLE.
C
C     SET VARIOUS PROGRAM VALUES
C
      WEIGHT = .FALSE.
      WT(1) = 1.0D0
      LWT = 1
C
      HLFRPT = .FALSE.
      APRXDV = .TRUE.
      PRTFXD = .TRUE.
      EXMPT = -1.0D0
      NETA = 0
C
      WIDE = .TRUE.
      PAGE = .FALSE.
C
      NDIGIT = 5
C
C     COMPUTE BACK OPERATORS
C
      CALL BACKOP(MSPEC, NFAC, NPARDF, MBOL, MBO, NPARMA, NPARAR)
      NNZW = N - NPARDF
C
C     SET UP FOR ERROR CHECKING
C
      IERR = 0
      NPARE = NPAR
      IF ((IFIXED(1).GE.0) .AND. (NPAR.GE.1)) NPARE =
     +   ICNTI(IFIXED,NPAR,0)
      IS = 0
      IF (STP(1).LE.0.0D0) IS = 1
C
      CALL LDSCMP(25, 0, MAX(IS*2*(N+NPAR),60+NPAR+NPARE) + 4*NFAC,
     +   0, 0, 0, 'D', 5*MBO +
     +   MAX(IS*(10*N+6*MBO+606),
     +        94+4*(N+MBO+101)+NPARE*(3*NPARE+35)/2),
     +   LDSMIN)
C
      CALL AMEER(NMSUB, N, NPAR, NPARE, LDSTAK,
     +   LDSMIN, STP, LSTP, SCALE, LSCALE, IVCV, SAVE, MSPEC, NFAC)
C
      IF (IERR.NE.0) RETURN
C
      CALL STKSET(LDSTAK, 4)
C
C     SET PRINT CONTROL VALUES
C
      CALL PRTCNT(NPRT, NDIGIT, IPTOUT)
C
C     SUBDIVIDE WORKSPACE FOR STEP SIZES
C
      NALL0 = STKST(1)
C
      IFP = 4
C
      STPT = STKGET(NPAR,IFP)
C
      PARDF = STKGET(MBO, IFP)
      PARAR = STKGET(MBO, IFP)
      PARMA = STKGET(MBO, IFP)
      T = STKGET(2*MBO, IFP)
C
      TEMP = T + MBO
C
      NFACT = NFAC
      MSPECT = STKGET(4*NFAC, 2)
C
C     SET UP FOR MODEL
C
      APRXDV = .TRUE.
      M = 1
      IXM = N
      NRESTS = MBO + 101 + N
C
      CALL CPYVII(NFAC, MSPEC(1,1), 4, ISTAK(MSPECT), 1)
      CALL CPYVII(NFAC, MSPEC(2,1), 4, ISTAK(MSPECT+NFAC), 1)
      CALL CPYVII(NFAC, MSPEC(3,1), 4, ISTAK(MSPECT+2*NFAC), 1)
      CALL CPYVII(NFAC, MSPEC(4,1), 4, ISTAK(MSPECT+3*NFAC), 1)
      CALL DCOEF (NFAC, ISTAK(MSPECT+NFAC), ISTAK(MSPECT+3*NFAC),
     +  NPARDF, RSTAK(PARDF), MBO, RSTAK(T))
C
C     COPY SUPPLIED STEP SIZES TO WORK SPACE
C
      CALL DCOPY(LSTP, STP, 1, RSTAK(STPT), 1)
C
      IF (IERR.NE.0) GO TO 10
C
C     SELECT STEP SIZES, IF DESIRED
C
      ISUBHD = 1
      IF (STP(1).LE.0.0D0) CALL AMESTP(Y, N, M, IXM, MDLTS3, PAR, NPAR,
     +  RSTAK(STPT), EXMPT, NETA, SCALE, LSCALE, IPTOUT(1), AMEHDR,
     +  PAGE, WIDE, ISUBHD, HLFRPT, PRTFXD, IFIXED, LIFIXD, STPAMO,
     +  NRESTS-N)
C
      CALL AMECNT(Y, WT, LWT, Y, N, M, IXM, MDLTS1, NLDRVN, APRXDV, DRV,
     +  PAR, NPAR, RES, IFIXED, LIFIXD, RSTAK(STPT), NPAR, MIT,
     +  STOPSS, STOPP, SCALE, LSCALE, DELTA, IVAPRX, RSD, PV, LPV,
     +  SDPV, LSDPV, SDRES, LSDRES, VCV, IVCV, WEIGHT, SAVE, NNZW,
     +  NPARE, AMEHDR, PAGE, WIDE, IPTOUT, NDIGIT, HLFRPT, NRESTS)
C
   10 CONTINUE
C
      CALL STKCLR(NALL0)
C
      RETURN
C
      END
