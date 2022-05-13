*NLSDC
      SUBROUTINE NLSDC(Y, XM, N, M, IXM, MDL, DRV, PAR, NPAR, RES,
     +   LDSTAK, IFIXED, IDRVCK, MIT, STOPSS, STOPP, SCALE, DELTA,
     +   IVAPRX, NPRT)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE USER CALLABLE SUBROUTINE FOR NONLINEAR LEAST
C     SQUARES REGRESSION USING ANALYTIC DERIVATIVES WITH USER
C     SUPPLIED CONTROL PARAMETERS.
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  APRIL 2, 1981
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      REAL
     +   DELTA,STOPP,STOPSS
      INTEGER
     +   IDRVCK,IVAPRX,IXM,LDSTAK,M,MIT,N,NPAR,NPRT
C
C  ARRAY ARGUMENTS
      REAL
     +   PAR(*),RES(*),SCALE(*),XM(*),Y(*)
      INTEGER
     +   IFIXED(*)
C
C  SUBROUTINE ARGUMENTS
      EXTERNAL DRV,MDL
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
     +   RSD
      INTEGER
     +   IPRT,IVCV,LIFIXD,LPV,LSCALE,LSDPV,LSDRES,LWT,NNZW,NPARE
      LOGICAL
     +   SAVE,WEIGHT
C
C  LOCAL ARRAYS
      REAL
     +   PV(1),SDPV(1),SDRES(1),VCV(1,1),WT(1)
      CHARACTER
     +   NMSUB(6)*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL IPRINT,NLCNTA
C
C  COMMON BLOCKS
      COMMON /CSTAK/DSTAK
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     REAL DELTA
C        THE MAXIMUM CHANGE ALLOWED IN THE MODEL PARAMETERS AT THE
C        FIRST ITERATION.
C     EXTERNAL DRV
C        THE NAME OF THE USER SUPPLIED SUBROUTINE WHICH COMPUTES THE
C        DERIVATIVE (JACOBIAN) MATRIX OF THE MODEL.
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER IDRVCK
C        THE VARIABLE USED TO INDICATE WHETHER OR NOT THE DERIVATIVES WE
C        CHECKED OR NOT.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .GE. 1, ERRORS WERE DETECTED.
C     INTEGER IFIXED(NPAR)
C        THE INDICATOR VALUES USED TO DESIGNATE WHETHER THE
C        PARAMETERS ARE TO BE OPTIMIZED OR ARE TO BE HELD FIXED.  IF
C        IFIXED(I).NE.0, THEN PAR(I) WILL BE OPTIMIZED.  IF
C        IFIXED(I).EQ.0, THEN PAR(I) WILL BE HELD FIXED.
C        IFIXED(1).LT.0, THEN ALL PAR(I),I=1,NPAR, WILL BE OPTIMIZED..
C     INTEGER IPRT
C        THE UNIT NUMBER FOR PRINTED OUTPUT.
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
C     INTEGER IXM
C        THE FIRST DIMENSION OF THE INDEPENDENT VARIABLE ARRAY.
C     INTEGER LDSTAK
C        THE LENGTH OF THE ARRAY DSTAK.
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
C     INTEGER LWT
C        THE ACTUAL LENGTH OF THE VECTOR WT.
C     INTEGER M
C        THE NUMBER OF INDEPENDENT VARIABLES.
C     INTEGER MIT
C        THE MAXIMUM NUMBER OF ITERATIONS ALLOWED.
C     EXTERNAL MDL
C        THE NAME OF THE USER SUPPLIED SUBROUTINE WHICH COMPUTES THE
C        PREDICTED VALUES BASED ON THE CURRENT PARAMETER ESTIMATE.   N
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS.
C     CHARACTER*1 NMSUB(6)
C        THE NAME OF THE SUBROUTINE CALLING THE ERROR CHECKING
C        SUBROUTINES.
C     INTEGER NNZW
C        THE NUMBER OF NON ZERO WEIGHTS.
C     INTEGER NPAR
C        THE NUMBER OF PARAMETERS IN THE MODEL.
C     INTEGER NPARE
C        THE NUMBER OF PARAMETERS TO BE ESTIMATED.
C     INTEGER NPRT
C        THE PARAMETER USED TO INDICATE HOW MUCH PRINTED OUTPUT IS
C        TO BE PROVIDED.
C     REAL PAR(NPAR)
C        THE CURRENT ESTIMATES OF THE PARAMETERS.
C     REAL PV(1)
C        A DUMMY ARRAY FOR
C        THE PREDICTED VALUE BASED ON THE CURRENT PARAMETER ESTIMATES
C     REAL RES(N)
C        THE RESIDUALS FROM THE FIT.
C     REAL RSD
C        THE VALUE OF THE RESIDUAL STANDARD DEVIATION AT THE SOLUTION.
C     LOGICAL SAVE
C        THE VARIABLE USED TO INDICATE WHETHER ANY RESULTS OTHER THAN
C        THE RESIDUALS AND PARAMETERS ARE TO BE SAVED (TRUE) OR NOT
C        (FALSE).
C     REAL SCALE(NPAR)
C        A VALUE TO INDICATE USE OF THE DEFAULT VALUES OF
C        THE TYPICAL SIZE OF THE PARAMETERS.
C     REAL SDPV(1)
C        A DUMMY ARRAY FOR
C        THE STANDARD DEVIATION OF THE PREDICTED VALUE.
C     REAL SDRES(1)
C        A DUMMY ARRAY FOR
C        THE STANDARD DEVIATIONS OF THE RESIDUALS.
C     REAL STOPP
C        THE STOPPING CRITERION FOR THE TEST BASED ON THE MAXIMUM SCALED
C        RELATIVE CHANGE IN THE ELEMENTS OF THE MODEL PARAMETER VECTOR
C     REAL STOPSS
C        THE STOPPING CRITERION FOR THE TEST BASED ON THE RATIO OF THE
C        PREDICTED DECREASE IN THE RESIDUAL SUM OF SQUARES (COMPUTED
C        BY STARPAC) TO THE CURRENT RESIDUAL SUM OF SQUARES ESTIMATE.
C     REAL VCV(1,1)
C        THE STARTING LOCATION IN RSTAK/DSTAK OF
C        THE COVARIANCE MATRIX.
C     LOGICAL WEIGHT
C        THE VARIABLE USED TO INDICATE WHETHER WEIGHTED ANALYSIS IS TO
C        BE PERFORMED (TRUE) OR NOT (FALSE).
C     REAL WT(1)
C        THE DUMMY ARRAY FOR THE WEIGHTS.
C     REAL XM(IXM,M)
C        THE ARRAY IN WHICH ONE ROW OF THE INDEPENDENT VARIABLE ARRAY
C        IS STORED.
C     REAL Y(N)
C        THE ARRAY OF THE DEPENDENT VARIABLE.
C
C     SET UP NAME ARRAYS
C
      DATA NMSUB(1), NMSUB(2), NMSUB(3), NMSUB(4), NMSUB(5), NMSUB(6) /
     +   'N','L','S','D','C',' '/
C
C     SET VARIOUS PROGRAM PARAMETERS
C
      WEIGHT = .FALSE.
      SAVE = .FALSE.
C
      WT(1) = 1.0E0
      LIFIXD = NPAR
      LPV = 1
      LSCALE = NPAR
      LSDPV = 1
      LSDRES = 1
      LWT = 1
      IVCV = 1
C
      CALL NLCNTA(Y, WT, LWT, XM, N, M, IXM, MDL, DRV, PAR, NPAR, RES,
     +   LDSTAK, IFIXED, LIFIXD, IDRVCK, MIT, STOPSS, STOPP, SCALE,
     +   LSCALE, DELTA, IVAPRX, NPRT, RSD, PV, LPV, SDPV, LSDPV, SDRES,
     +   LSDRES, VCV, IVCV, NMSUB, WEIGHT, SAVE, NNZW, NPARE)
C
      IF (IERR.NE.1) RETURN
C
C     PRINT PROPER CALL SEQUENCE
C
      CALL IPRINT(IPRT)
      WRITE (IPRT,1000)
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (/42H THE CORRECT FORM OF THE CALL STATEMENT IS//
     +   '       CALL NLSDC (Y, XM, N, M, IXM, NLSMDL, NLSDRV,'/
     +   '      +            PAR, NPAR, RES, LDSTAK,'/
     +   '      +            IFIXED, IDRVCK, MIT, STOPSS, STOPP,'/
     +   '      +            SCALE, DELTA, IVAPRX, NPRT)')
      END
