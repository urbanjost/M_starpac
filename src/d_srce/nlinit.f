*NLINIT
      SUBROUTINE NLINIT (N, IFIXD, PAR, NPAR, PARE, NPARE, MIT,
     +   STOPSS, STOPP, SCALE, LSCALE, DELTA, IVAPRX, APRXDV, IVCVPT,
     +   IWORK, IIWORK, RWORK, IRWORK, SCL)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE PERFORMS INITIALIZATION FOR THE NONLINEAR
C     LEAST SQUARES ROUTINES.
C
C     REFERENCES
C
C     DENNIS, J.E., GAY, D.M., AND WELSCH, R.E. (1979), AN ADAPTIVE
C             NONLINEAR LEAST-SQUARES ALGORITHM, (BEING REVISED).
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
     +   DELTA,STOPP,STOPSS
      INTEGER
     +   IIWORK,IRWORK,IVAPRX,IVCVPT,LSCALE,MIT,N,NPAR,NPARE,SCL
      LOGICAL
     +   APRXDV
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   PAR(NPAR),PARE(NPAR),RWORK(IRWORK),SCALE(LSCALE)
      INTEGER
     +   IFIXD(NPAR),IWORK(IIWORK)
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      INTEGER
     +   AFCTOL,CNVCOD,COVPRT,COVREQ,DINIT,DTYPE,ISCL,J,LMAX0,
     +   MXFCAL,MXITER,NITER,OUTLEV,PRUNIT,RFCTOL,SCLJ,SOLPRT,
     +   STATPR,X0PRT,XCTOL
C
C  EXTERNAL FUNCTIONS
      DOUBLE PRECISION
     +   RMDCON
      EXTERNAL RMDCON
C
C  EXTERNAL SUBROUTINES
      EXTERNAL DFAULT,NLSPK
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,IABS,MAX
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     INTEGER AFCTOL
C        THE LOCATION IN RWORK OF THE ABSOLUTE CONVERGENCE TOLERANCE.
C     LOGICAL APRXDV
C        THE VARIABLE USED TO INDICATE WHETHER NUMERICAL
C        APPROXIMATIONS TO THE DERIVATIVE WERE USED (TRUE) OR NOT
C        (FALSE).
C     INTEGER CNVCOD
C        A VALUE USED TO CONTROL THE PRINTING OF ITERATION REPORTS.
C     INTEGER COVPRT
C        THE LOCATION IN IWORK OF THE VARIABLE USED TO INDICATE WHETHER
C        THE COVARIANCE MATRIX IS TO BE PRINTED BY THE NL2 CODE, WHERE
C        IWORK(COVPRT) = 0 INDICATES IT IS NOT.
C     INTEGER COVREQ
C        THE LOCATION IN IWORK OF THE VARIABLE USED TO INDICATE HOW
C        THE COVARIANCE MATRIX IS TO BE COMPUTED BY THE NL2 CODE, WHERE
C        IWORK(COVREQ) = 3 INDICATES THE COVARIANCE MATRIX IS TO BE COMP
C        AS THE RESIDUAL VARIANCE TIMES THE INVERSE OF THE JACOBIAN MATR
C        TRANSPOSED TIMES THE JACOBIAN MATRIX .
C     DOUBLE PRECISION DELTA
C        THE MAXIMUM CHANGE ALLOWED IN THE MODEL PARAMETERS AT THE
C        FIRST ITERATION.
C     INTEGER DINIT
C        THE LOCATION IN IWORK OF THE VALUE USED TO INDICATE
C        WHETHER OR NOT USER SUPPLIED SCALE VALUES ARE TO BE
C        USED, WHERE THE (NL2) DEFAULT VALUE OF RWORK(DINIT) = 0.0D0
C        INIDCATES NO, AND THE VALUE RWORK(DINIT) = -1.0D0 INDICATES
C        YES.
C     INTEGER DTYPE
C        THE LOCATION IN IWORK OF THE VALUE INDICATING WHETHER THE
C        SCALE VALUES HAVE BEEN SUPPLIED BY THE USER (IWORK(DTYPE) .LE.
C        OR THE DEFAULT VALUES ARE TO BE USED (IWORK(DTYPE) .GT. 0).
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .GE. 1, ERRORS WERE DETECTED.
C     INTEGER IFIXD(NPAR)
C        THE INDICATOR VALUES USED TO DESIGNATE WHETHER THE
C        PARAMETERS ARE TO BE OPTIMIZED OR ARE TO BE HELD FIXED.
C        IF IFIXED(I).NE.0, THEN PAR(I) WILL BE HELD FIXED.
C        IF IFIXED(I).EQ.0, THEN PAR(I) WILL BE OPTIMIZED.
C     INTEGER IIWORK
C        THE DIMENSION OF THE INTEGER WORK VECTOR IWORK.
C     INTEGER IRWORK
C        THE DIMENSION OF THE DOUBLE PRECISION WORK VECTOR RWORK.
C     INTEGER ISCL
C        THE LOCATION IN IWORK INDICATING THE STARTING LOCATION IN
C         RWORK OF THE SCALE VECTOR.
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
C     INTEGER IVCVPT
C        AN INDICATOR VALUE USED TO DESIGNATE WHICH FORM OF THE
C        VARIANCE COVARIANCE MATRIX (VCV) IS BEING PRINTED, WHERE
C        IVCVPT = 1 INDICATES THE VCV WAS COMPUTED AS
C                   INVERSE(TRANSPOSE(JACOBIAN)*JACOBIAN)
C        IVCVPT = 2 INDICATES THE VCV WAS COMPUTED AS
C                   INVERSE(HESSIAN)
C        IVCVPT = 3 INDICATES THE VCV WAS COMPUTED AS
C                   INVERSE(HESSIAN)*TRANSPOSE(JACOBIAN)*JACOBIAN
C                       *INVERSE(HESSIAN)
C     INTEGER IWORK(IIWORK)
C        THE INTEGER WORK SPACE VECTOR USED BY THE NL2 SUBROUTINES.
C     INTEGER J
C        THE INDEX OF THE PARAMETER BEING EXAMINED.
C     INTEGER LMAX0
C        THE LOCATION IN RWORK OF THE VALUE INDICATING THE
C        MAXIMUM CHANGE ALLOWED IN THE MODEL PARAMETERS AT THE
C        FIRST ITERATION.
C     INTEGER MIT
C        THE MAXIMUM NUMBER OF ITERATIONS ALLOWED.
C     INTEGER MXFCAL
C        THE LOCATION IN IWORK OF THE VARIABLE DESIGNATING THE
C        MAXIMUM NUMBER OF FUNCTION CALLS ALLOWED, EXCLUDING
C        CALLS NECESSARY TO COMPUTE THE DERIVATIVES AND VARIANCE
C        COVARIANCE MATRIX.
C     INTEGER MXITER
C        THE LOCATION IN IWORK OF THE VARIABLE DESIGNATING THE
C        MAXIMUM NUMBER OF ITERATIONS ALLOWED.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS.
C     INTEGER NITER
C        THE LOCATION IN IWORK OF THE NUMBER OF THE CURRENT ITERATION.
C     INTEGER NPAR
C        THE NUMBER OF UNKNOWN PARAMETERS IN THE MODEL.
C     INTEGER NPARE
C        THE NUMBER OF UNKNOWN PARAMETERS TO BE OPTIMIZED.
C     INTEGER OUTLEV
C        THE LOCATION IN IWORK OF THE PARAMETER USED TO CONTROL THE
C        PRINTING OF THE ITERATION REPORTS BY NL2.
C     DOUBLE PRECISION PAR(NPAR)
C        THE ARRAY IN WHICH THE CURRENT ESTIMATES OF THE UNKNOWN
C        PARAMETERS ARE STORED.
C     DOUBLE PRECISION PARE(NPAR)
C        THE CURRENT ESTIMATES OF THE UNKNOWN PARAMETERS, BUT ONLY
C        THOSE TO BE OPTIMIZED (NOT THOSE WHOSE VALUES ARE FIXED).
C     INTEGER PRUNIT
C        THE LOCATION IN IWORK OF THE PARAMETER USED TO CONTROL
C        THE PRINT UNIT USED BY NL2.  IWORK(PRUNIT) = 0 MEANS
C        DONT PRINT ANYTHING.
C     INTEGER RFCTOL
C        THE LOCATION IN RWORK OF THE RELATIVE FUNCTION CONVERGENCE
C        TOLERANCE.
C     DOUBLE PRECISION RWORK(IRWORK)
C        THE DOUBLE PRECISION WORK VECTOR USED BY THE NL2 SUBROUTINES.
C     DOUBLE PRECISION SCALE(LSCALE)
C        THE TYPICAL SIZE OF THE UNKNOWN PARAMETERS.
C     INTEGER SCL
C        THE INDEX IN RWORK OF THE 1ST VALUE OF THE USER SUPPLIED SCALE
C        VALUE.
C     INTEGER SCLJ
C        THE INDEX IN RWORK OF THE JTH VALUE OF THE USER SUPPLIED SCALE
C        VALUE.
C     INTEGER SOLPRT
C        THE LOCATION IN IWORK OF THE PARAMETER USED TO CONTROL PRINTING
C        BY NL2 OF THE FINAL SOLUTION.
C     INTEGER STATPR
C        THE LOCATION IN IWORK OF THE PARAMETER USED TO CONTROL PRINTING
C        BY NL2 OF SUMMARY STATISTICS.
C     DOUBLE PRECISION STOPP
C        THE STOPPING CRITERION FOR THE TEST BASED ON THE MAXIMUM SCALED
C        RELATIVE CHANGE IN THE ELEMENTS OF THE MODEL PARAMETER VECTOR
C     DOUBLE PRECISION STOPSS
C        THE STOPPING CRITERION FOR THE TEST BASED ON THE RATIO OF THE
C        PREDICTED DECREASE IN THE RESIDUAL SUM OF SQUARES (COMPUTED
C        BY STARPAC) TO THE CURRENT RESIDUAL SUM OF SQUARES ESTIMATE.
C     INTEGER XCTOL
C        THE LOCATION IN RSTAK/DSTAK OF THE P CONVERGENCE TOLERANCE.
C     INTEGER X0PRT
C         THE LOCATION IN IWORK OF THE PARAMETER USED TO CONTROL PRINTIN
C        BY NL2 OF THE INITIAL PARAMETER AND SCALE VALUES.
C
C     IWORK SUBSCRIPT VALUES
C
      DATA CNVCOD /34/, COVPRT /14/, COVREQ /15/, DINIT /38/, DTYPE
     +   /16/, ISCL /27/, MXFCAL /17/, MXITER /18/,
     +   NITER /31/, OUTLEV /19/, PRUNIT /21/, SOLPRT /22/, STATPR
     +   /23/, X0PRT /24/
C
C     RWORK SUBSCRIPT VALUES
C
      DATA AFCTOL /31/, LMAX0 /35/, RFCTOL /32/, XCTOL /33/
C
C+++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++
C
C     PACK PARAMETERS INTO PARE
C
      CALL NLSPK(PAR, IFIXD, NPAR, PARE, NPAR)
C
C     SET NL2SOL DEFAULT VALUES
C
      CALL DFAULT(IWORK, RWORK)
C
C     SET NON NL2 DEFAULT VALUES
C
      IWORK(MXITER) = MIT
      IF (MIT.LE.0) IWORK(MXITER) = 21
C
      IWORK(MXFCAL) = 2*IWORK(MXITER)
C
C     SET STOPPING CRITERION
C
      RWORK(AFCTOL) = RMDCON(1)
      IF ((STOPSS.GE.RMDCON(3)) .AND. (STOPSS.LE.0.1)) RWORK(RFCTOL) =
     +   STOPSS
C
      IF ((STOPP.GE.0.0D0) .AND. (STOPP.LE.1.0D0))
     +   RWORK(XCTOL) = STOPP
C
C     SET SCALE VALUES
C
      SCL = 94 + 2*N + NPARE*(3*NPARE+31)/2
      IWORK(ISCL) = SCL
      IF (SCALE(1).GT.0.0D0) GO TO 40
C
      IWORK(DTYPE) = 1
C
C     INITIALIZE SCALE VALUES FOR FIRST ITERATION
C
      SCLJ = SCL - 1
      DO 30 J=1,NPAR
         IF (IFIXD(J).NE.0) GO TO 30
         SCLJ = SCLJ + 1
         IF (PAR(J).EQ.0.0D0) RWORK(SCLJ) = 1.0D0
         IF (PAR(J).NE.0.0D0) RWORK(SCLJ) = 1.0D0/ABS(PAR(J))
   30 CONTINUE
C
      GO TO 60
C
   40 IWORK(DTYPE) = 0
      RWORK(DINIT) = -1.0D0
      SCLJ = SCL - 1
      DO 50 J=1,NPAR
         IF (IFIXD(J).NE.0) GO TO 50
         SCLJ = SCLJ + 1
         RWORK(SCLJ) = 1.0D0/MAX(ABS(SCALE(J)),ABS(PAR(J)))
   50 CONTINUE
C
   60 IF (DELTA.LE.0.0D0) RWORK(LMAX0) = 100.0D0
      IF (DELTA.GT.0.0D0) RWORK(LMAX0) = DELTA
C
C     SET NL2 COVARIANCE COMPUTATION CONTROL PARAMETER
C
      IF ((IVAPRX.LE.1) .OR. (IVAPRX.EQ.4) .OR. (IVAPRX.GE.7))
     +   IWORK(COVREQ) = 3
      IF ((IVAPRX.EQ.2) .OR. (IVAPRX.EQ.5)) IWORK(COVREQ) = 2
      IF ((IVAPRX.EQ.3) .OR. (IVAPRX.EQ.6)) IWORK(COVREQ) = 1
      IF ((IVAPRX.GE.4) .AND. (IVAPRX.LE.6))
     +   IWORK(COVREQ) = -IWORK(COVREQ)
      IF (APRXDV) IWORK(COVREQ) = -IABS(IWORK(COVREQ))
      IF ((IVAPRX.LE.1) .OR. (IVAPRX.EQ.4) .OR. (IVAPRX.GE.7))
     +   IVCVPT = 1
      IF ((IVAPRX.EQ.2) .OR. (IVAPRX.EQ.5)) IVCVPT = 2
      IF ((IVAPRX.EQ.3) .OR. (IVAPRX.EQ.6)) IVCVPT = 3
C
C     INITIALIZE THE ITERATION COUNTER
C
      IWORK(NITER) = 0
C
C     SET NL2 PRINT CONTROL PARAMETERS
C
      IWORK(CNVCOD) = 0
      IWORK(COVPRT) = 0
      IWORK(OUTLEV) = 0
      IWORK(PRUNIT) = 0
      IWORK(SOLPRT) = 0
      IWORK(STATPR) = 0
      IWORK(X0PRT) = 0
C
      RETURN
C
      END
