*LLSPWS
      SUBROUTINE LLSPWS(Y, WT, XM, N, NDEG, RES, LDSTAK,
     +   NPRT, LPAR, PAR, NPAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     CALL FOR POLYNOMIAL MODEL LEAST SQUARES FIT
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
      DOUBLE PRECISION
     +   RSD
      INTEGER
     +   IVCV,LDSTAK,LPAR,N,NDEG,NPAR,NPRT
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
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
     +   IPRT,LPV,LSDPV,LSDRES,LWT
      LOGICAL
     +   SAVE,WEIGHT
C
C  LOCAL ARRAYS
      CHARACTER
     +   NMSUB(6)*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL IPRINT,LLCNTP
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
C     INTEGER NDEG
C        THE DEGREE OF THE MODEL.
C     CHARACTER*1 NMSUB(6)
C        THE NAME OF THE SUBROUTINE CALLING THE ERROR CHECKING
C        SUBROUTINES.
C     INTEGER NPAR
C        THE NUMBER OF PARAMETERS.
C     INTEGER NPRT
C        THE INDICATOR VARIABLE USED TO DESIGNATE THE AMOUNT OF
C        PRINTED OUTPUT.
C     DOUBLE PRECISION PAR(LPAR)
C        THE PARAMETERS  TO BE ESTIMATED.
C     DOUBLE PRECISION PV(N)
C        THE PREDICTED VALUES.
C     DOUBLE PRECISION RES(N)
C        THE RESIDUALS.
C     DOUBLE PRECISION RSD
C        THE RESIDUAL STANDARD DEVIATION.
C     LOGICAL SAVE
C        THE VARIABLE USED TO INDICATE WHETHER ANY RESULTS OTHER THAN
C        THE RESIDUALS ARE TO VE SAVED (TRUE) OR NOT (FALSE).
C     DOUBLE PRECISION SDPV(N)
C        THE STANDARD DEVIATIONS OF THE PREDICTED VALUES.
C     DOUBLE PRECISION SDRES(N)
C        THE STANDARDIZED RESIDUALS.
C     DOUBLE PRECISION VCV(IVCV,NPAR)
C        THE VARIANCE COVARIANCE MATRIX.
C     LOGICAL WEIGHT
C        THE VARIABLE USED TO INDICATE WHETHER WEIGHTED ANALYSIS IS TO
C        BE PERFORMED (TRUE) OR NOT (FALSE).
C     DOUBLE PRECISION WT(N)
C        THE WEIGHTS (A DUMMY VECTOR IN THE UNWEIGHTED CASE).
C     DOUBLE PRECISION XM(N,1)
C        THE INDEPENDENT VARIABLE.
C     DOUBLE PRECISION Y(N)
C        THE DEPENDENT VARIABLE.
C
C
      DATA NMSUB(1), NMSUB(2), NMSUB(3), NMSUB(4), NMSUB(5), NMSUB(6) /
     +   'L','L','S','P','W','S'/
C
C     SET PARAMETERS NECESSARY FOR THE COMPUTATIONS
C
      WEIGHT = .TRUE.
      SAVE = .TRUE.
      LPV = N
      LSDPV = N
      LSDRES = N
      LWT = N
C
      CALL LLCNTP(Y, WT, LWT, XM, N, NDEG, NPAR, RES, LDSTAK, NPRT,
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
     +   '       CALL LLSPWS (Y, WT, X, N, NDEG, RES, LSDTAK,'/
     +   '      +             NPRT, LPAR, PAR, NPAR, RSD, PV, SDPV,'/
     +   '      +             SDRES, VCV, IVCV)')
      END
