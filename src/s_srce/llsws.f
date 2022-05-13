*LLSWS
      SUBROUTINE LLSWS(Y, WT, XM, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     CALL FOR GENERAL LINEAR MODEL LEAST SQUARES FIT
C     USER SUPPLIED WEIGHTS SPECIFIED
C     FULL STORAGE
C     USER CONTROL OF AUTOMATIC PRINTOUT
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
     +   IVCV,IXM,LDSTAK,N,NPAR,NPRT
C
C  ARRAY ARGUMENTS
      REAL
     +   PAR(*),PV(*),RES(*),SDPV(*),SDRES(*),VCV(*),WT(*),XM(*),Y(*)
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
     +   IPRT,LPAR,LPV,LSDPV,LSDRES,LWT
      LOGICAL
     +   SAVE,WEIGHT
C
C  LOCAL ARRAYS
      CHARACTER
     +   NMSUB(6)*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL IPRINT,LLCNTG
C
C  COMMON BLOCKS
      COMMON /CSTAK/DSTAK
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER IERR
C        THE INTEGER VALUE DESIGNATING WHETHER ANY ERRORS WERE
C        DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .EQ. 1, ERRORS WERE DETECTED.
C     INTEGER IPRT
C        THE UNIT NUMBER FOR PRINTED OUTPUT.
C     INTEGER IVCV
C        THE FIRST DIMENSION OF THE MATRIX VCV.
C     INTEGER IXM
C        THE FIRST DIMENSION OF THE MATRIX XM.
C     INTEGER LDSTAK
C        THE LENGTH OF THE VECTOR DSTAK IN COMMON CSTAK.
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
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS.
C     CHARACTER*1 NMSUB(6)
C        THE NAME OF THE SUBROUTINE CALLING THE ERROR CHECKING
C        SUBROUTINES.
C     INTEGER NPAR
C        THE NUMBER OF PARAMETERS.
C     INTEGER NPRT
C        THE INDICATOR VARIABLE USED TO DESIGNATE THE AMOUNT OF
C        PRINTED OUTPUT.
C     REAL PAR(NPAR)
C        THE PARAMETERS  TO BE ESTIMATED.
C     REAL PV(N)
C        THE PREDICTED VALUES.
C     REAL RES(N)
C        THE RESIDUALS.
C     REAL RSD
C        THE RESIDUAL STANDARD DEVIATION.
C     LOGICAL SAVE
C        THE VARIABLE USED TO INDICATE WHETHER ANY RESULTS OTHER THAN
C        THE RESIDUALS ARE TO VE SAVED (TRUE) OR NOT (FALSE).
C     REAL SDPV(N)
C        THE STANDARD DEVIATIONS OF THE PREDICTED VALUES.
C     REAL SDRES(N)
C        THE STANDARDIZED RESIDUALS.
C     REAL VCV(IVCV,NPAR)
C        THE VARIANCE COVARIANCE MATRIX.
C     LOGICAL WEIGHT
C        THE VARIABLE USED TO INDICATE WHETHER WEIGHTED ANALYSIS IS TO
C        BE PERFORMED (TRUE) OR NOT (FALSE).
C     REAL WT(N)
C        THE WEIGHTS (A DUMMY VECTOR IN THE UNWEIGHTED CASE).
C     REAL XM(IXM,NPAR)
C        THE INDEPENDENT VARIABLE.
C     REAL Y(N)
C        THE DEPENDENT VARIABLE.
C
C
      DATA NMSUB(1), NMSUB(2), NMSUB(3), NMSUB(4), NMSUB(5), NMSUB(6) /
     +   'L','L','S','W','S',' '/
C
C     SET PARAMETERS NECESSARY FOR THE COMPUTATIONS
C
      WEIGHT = .TRUE.
      SAVE = .TRUE.
      LPAR = NPAR
      LPV = N
      LSDPV = N
      LSDRES = N
      LWT = N
C
      CALL LLCNTG(Y, WT, LWT, XM, N, IXM, NPAR, RES, LDSTAK, NPRT,
     +   PAR, LPAR, RSD, PV, LPV, SDPV, LSDPV, SDRES, LSDRES, VCV, IVCV,
     +   NMSUB, WEIGHT, SAVE)
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
 1000 FORMAT (//42H THE CORRECT FORM OF THE CALL STATEMENT IS//
     + '       CALL LLSWS (Y, WT, XM, N, IXM, NPAR, RES, LSDTAK,'/
     + '      +            NPRT, PAR, RSD, PV, SDPV, SDRES, VCV, IVCV)')
      END
