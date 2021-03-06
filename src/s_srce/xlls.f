*XLLS
      SUBROUTINE XLLS(LDS)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     TEST ROUTINES FOR LINEAR LEAST SQUARES SUBROUTINES.
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
      INTEGER
     +   LDS
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
     +   RSD,SUM,TERM
      INTEGER
     +   I,IPRT,IVCV,IXM,J,LDSMIN,LDSTAK,LPAR,N,NDEG,NPAR,NPRT
C
C  LOCAL ARRAYS
      REAL
     +   PAR(10),PV(50),RAND(1),RES(50),SDPV(50),SDRES(50),VCV(10,10),
     +   WT(50),X(50,9),XM(50,10),XM1(50,10),Y(50),Y1(50)
C
C  EXTERNAL SUBROUTINES
      EXTERNAL FITSXP,GENR,IPRINT,LDSCMP,LLS,LLSP,LLSPS,LLSPW,LLSPWS,
     +   LLSS,LLSW,LLSWS,NRAND,SETRV
C
C  COMMON BLOCKS
      COMMON /CSTAK/DSTAK
      COMMON /ERRCHK/IERR
C
C  EQUIVALENCES
      EQUIVALENCE (XM(1,2),X(1,1))
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER I
C        AN INDEX.
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
C        THE FIRST DIMENSION OF THE MATRIX X.
C     INTEGER J
C        AN INDEX.
C     INTEGER LDS
C       ..
C     INTEGER LDSMIN
C        THE MINIMUM LENGTH ALLOWED FOR THE ARRAY DSTAK.
C     INTEGER LDSTAK
C        THE LENGTH OF THE VECTOR DSTAK IN COMMON CSTAK.
C     INTEGER LPAR
C        THE ACTUAL LENGTH OF THE PARAMETER ARRAY.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS.
C     INTEGER NDEG
C        THE DEGREE OF THE POLYNOMIAL MODEL TO BE FIT.
C     INTEGER NPAR
C        THE NUMBER OF PARAMETERS.
C     INTEGER NPRT
C        THE INDICATOR VARIABLE USED TO DESIGNATE THE AMOUNT OF
C        PRINTED OUTPUT.
C     REAL PAR(10)
C        THE PARAMETERS  TO BE ESTIMATED.
C     REAL PV(50)
C        THE PREDICTED VALUES.
C     REAL RAND(1)
C        *
C     REAL RES(50)
C        THE RESIDUALS.
C     REAL RSD
C        THE RESIDUAL STANDARD DEVIATION.
C     REAL SDPV(50)
C        THE STANDARD DEVIATIONS OF THE PREDICTED VALUES.
C     REAL SDRES(50)
C        THE STANDARDIZED RESIDUALS.
C     REAL SUM, TERM
C        *
C     REAL VCV(10,10)
C        THE VARIANCE COVARIANCE MATRIX.
C     REAL WT(50)
C        THE WEIGHTS (A DUMMY VECTOR IN THE UNWEIGHTED CASE).
C     REAL X(50,9)
C        THE INDEPENDENT VARIABLE.
C     REAL XM(50,10)
C        THE INDEPENDENT VARIABLE.
C     REAL XM1(50,10)
C        THE INDEPENDENT VARIABLE.
C     REAL Y(50)
C        THE DEPENDENT VARIABLE.
C     REAL Y1(50)
C        THE DEPENDENT VARIABLE.
C
C
      DATA      XM(1,1),  XM(1,2),  XM(1,3),  XM(1,4)
     +    /      1.0E0, 42.2E0, 11.2E0, 31.9E0/
      DATA      XM(2,1),  XM(2,2),  XM(2,3),  XM(2,4)
     +    /      1.0E0, 48.6E0, 10.6E0, 13.2E0/
      DATA      XM(3,1),  XM(3,2),  XM(3,3),  XM(3,4)
     +    /      1.0E0, 42.6E0, 10.6E0, 28.7E0/
      DATA      XM(4,1),  XM(4,2),  XM(4,3),  XM(4,4)
     +    /      1.0E0, 39.0E0, 10.4E0, 26.1E0/
      DATA      XM(5,1),  XM(5,2),  XM(5,3),  XM(5,4)
     +    /      1.0E0, 34.7E0,  9.3E0, 30.1E0/
      DATA      XM(6,1),  XM(6,2),  XM(6,3),  XM(6,4)
     +    /      1.0E0, 44.5E0, 10.8E0,  8.5E0/
      DATA      XM(7,1),  XM(7,2),  XM(7,3),  XM(7,4)
     +    /      1.0E0, 39.1E0, 10.7E0, 24.3E0/
      DATA      XM(8,1),  XM(8,2),  XM(8,3),  XM(8,4)
     +    /      1.0E0, 40.1E0, 10.0E0, 18.6E0/
      DATA      XM(9,1),  XM(9,2),  XM(9,3),  XM(9,4)
     +    /      1.0E0, 45.9E0, 12.0E0, 20.4E0/
      DATA         Y(1),     Y(2),     Y(3)
     +    /    167.1E0,174.4E0,160.8E0/
      DATA         Y(4),     Y(5),     Y(6)
     +    /    162.0E0,140.8E0,174.6E0/
      DATA         Y(7),     Y(8),     Y(9)
     +    /    163.7E0,174.5E0,185.7E0/
C
C     SET PARAMETERS NECESSARY FOR THE COMPUTATIONS
C
      CALL IPRINT(IPRT)
      N = 9
      NPAR = 4
      NDEG = 3
      NPRT = 2
      LPAR = 10
      IVCV = 10
      IXM = 50
      LDSTAK = LDS
C
      CALL SETRV(WT, N, 1.0E0)
C
C     CHECK ERROR HANDLING
C
C        ERROR 1  -  NON POSITIVE NUMBER OF OBSERVATIONS AND PARAMETER
C                    NUMBER OF PARAMETERS GREATER THAN N
C                    IXM LESS THAN NUMBER OF OBSERVATIONS
C                    IVCV LESS THAN NUMBER OF PARAMETERS
C                    LPAR TOO SMALL
C
      N = -5
      NPAR = 0
      NDEG = -1
      IXM = -10
      LPAR = -1
      IVCV = -10
      NPRT = -1
      WRITE (IPRT,1200)
      WRITE (IPRT,1000)
      CALL LLS(Y, XM, N, IXM, NPAR, RES, LDSTAK)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1010)
      CALL LLSS(Y, XM, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1020)
      CALL LLSW(Y, WT, XM, N, IXM, NPAR, RES, LDSTAK)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1030)
      CALL LLSWS(Y, WT, XM, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV,SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1040)
      CALL LLSP(Y, X, N, NDEG, RES, LDSTAK)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1050)
      CALL LLSPS(Y, X, N, NDEG, RES, LDSTAK,
     +   NPRT, LPAR, PAR, NPAR, RSD, PV,SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1060)
      CALL LLSPW(Y, WT, X, N, NDEG, RES, LDSTAK)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1070)
      CALL LLSPWS(Y, WT, X, N, NDEG, RES, LDSTAK,
     +   NPRT, LPAR, PAR, NPAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      N = 9
      NPAR = 4
      NDEG = 3
      IXM = 50
      LPAR = -10
      IVCV = 10
C
C        ERROR 2  -  LDS TOO SMALL
C                    LPAR TOO SMALL
C
      LDSTAK = 0
      WRITE (IPRT,1220)
      WRITE (IPRT,1000)
      CALL LLS(Y, XM, N, IXM, NPAR, RES, LDSTAK)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1010)
      CALL LLSS(Y, XM, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1020)
      CALL LLSW(Y, WT, XM, N, IXM, NPAR, RES, LDSTAK)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1030)
      CALL LLSWS(Y, WT, XM, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV,SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1040)
      CALL LLSP(Y, X, N, NDEG, RES, LDSTAK)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1050)
      CALL LLSPS(Y, X, N, NDEG, RES, LDSTAK,
     +   NPRT, LPAR, PAR, NPAR, RSD, PV,SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1060)
      CALL LLSPW(Y, WT, X, N, NDEG, RES, LDSTAK)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1070)
      CALL LLSPWS(Y, WT, X, N, NDEG, RES, LDSTAK,
     +   NPRT, LPAR, PAR, NPAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      LDSTAK = LDS
      NPRT = 2
      LPAR = 10
C
C        ERROR 3  -  NEGATIVE WEIGHTS
C
      WT(1) = -1.0E0
      WRITE (IPRT,1240)
      WRITE (IPRT,1020)
      CALL LLSW(Y, WT, XM, N, IXM, NPAR, RES, LDSTAK)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1030)
      CALL LLSWS(Y, WT, XM, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV,SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1060)
      CALL LLSPW(Y, WT, X, N, NDEG, RES, LDSTAK)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1070)
      CALL LLSPWS(Y, WT, X, N, NDEG, RES, LDSTAK,
     +   NPRT, LPAR, PAR, NPAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      WT(1) = 1.0E0
C
C        ERROR 4  -  TOO FEW POSITIVE WEIGHTS
C
      CALL SETRV(WT(2), N-1, 0.0E0)
      WRITE (IPRT,1250)
      WRITE (IPRT,1020)
      CALL LLSW(Y, WT, XM, N, IXM, NPAR, RES, LDSTAK)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1030)
      CALL LLSWS(Y, WT, XM, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV,SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1060)
      CALL LLSPW(Y, WT, X, N, NDEG, RES, LDSTAK)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1070)
      CALL LLSPWS(Y, WT, X, N, NDEG, RES, LDSTAK,
     +   NPRT, LPAR, PAR, NPAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL SETRV(WT(2), N-1, 1.0E0)
C
C     CHECK RESULTS FROM VALID CALL
C
      WRITE (IPRT,1260)
      WRITE (IPRT,1000)
      CALL LLS(Y, XM, N, IXM, NPAR, RES, LDSTAK)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1430) (RES(I), I = 1, N)
      WRITE (IPRT,1260)
      WRITE (IPRT,1010)
      CALL LLSS(Y, XM, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
      WRITE (IPRT,1260)
      WRITE (IPRT,1020)
      CALL LLSW(Y, WT, XM, N, IXM, NPAR, RES, LDSTAK)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1430) (RES(I), I = 1, N)
      WRITE (IPRT,1260)
      WRITE (IPRT,1030)
      CALL LLSWS(Y, WT, XM, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV,SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
      WRITE (IPRT,1260)
      WRITE (IPRT,1040)
      CALL LLSP(Y, X, N, NDEG, RES, LDSTAK)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1430) (RES(I), I = 1, N)
      WRITE (IPRT,1260)
      WRITE (IPRT,1050)
      CALL LLSPS(Y, X, N, NDEG, RES, LDSTAK,
     +   NPRT, LPAR, PAR, NPAR, RSD, PV,SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
      WRITE (IPRT,1260)
      WRITE (IPRT,1060)
      CALL LLSPW(Y, WT, X, N, NDEG, RES, LDSTAK)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1430) (RES(I), I = 1, N)
      WRITE (IPRT,1260)
      WRITE (IPRT,1070)
      CALL LLSPWS(Y, WT, X, N, NDEG, RES, LDSTAK,
     +   NPRT, LPAR, PAR, NPAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
C
C     CHECK RESULTS FROM EXACT FIT
C
      N = NPAR
      NDEG = NPAR-1
C
      WRITE (IPRT,1270)
      WRITE (IPRT,1000)
      CALL LLS(Y, XM, N, IXM, NPAR, RES, LDSTAK)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1430) (RES(I), I = 1, N)
      WRITE (IPRT,1270)
      WRITE (IPRT,1010)
      CALL LLSS(Y, XM, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
      WRITE (IPRT,1270)
      WRITE (IPRT,1040)
      CALL LLSP(Y, X, N, NDEG, RES, LDSTAK)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1430) (RES(I), I = 1, N)
      WRITE (IPRT,1270)
      WRITE (IPRT,1050)
      CALL LLSPS(Y, X, N, NDEG, RES, LDSTAK,
     +   NPRT, LPAR, PAR, NPAR, RSD, PV,SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
C
      N = 9
C
      CALL SETRV(WT(NPAR+1), N-NPAR, 0.0E0)
C
      WRITE (IPRT,1270)
      WRITE (IPRT,1020)
      CALL LLSW(Y, WT, XM, N, IXM, NPAR, RES, LDSTAK)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1430) (RES(I), I = 1, N)
      WRITE (IPRT,1270)
      WRITE (IPRT,1030)
      CALL LLSWS(Y, WT, XM, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV,SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
      WRITE (IPRT,1270)
      WRITE (IPRT,1060)
      CALL LLSPW(Y, WT, X, N, NDEG, RES, LDSTAK)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1430) (RES(I), I = 1, N)
      WRITE (IPRT,1270)
      WRITE (IPRT,1070)
      CALL LLSPWS(Y, WT, X, N, NDEG, RES, LDSTAK,
     +   NPRT, LPAR, PAR, NPAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
C
      CALL SETRV(WT(NPAR+1), N-NPAR, 1.0E0)
C
C     CHECK RESULTS FROM RANK DEFICIENT FIT
C
      DO 10 I = 1, N
         XM(I,5) = XM(I,4)
   10 CONTINUE
      WRITE (IPRT,1280)
      WRITE (IPRT,1000)
      CALL LLS(Y, XM, N, IXM, NPAR+1, RES, LDSTAK)
      WRITE (IPRT,1500) IERR
C
C     CHECK RESULTS FROM A POORLY SCALED PROBLEM.
C
      DO 30 I = 1, N
         Y1(I) = Y(I) * 1.0E-8
         DO 20 J = 1, 4
            XM1(I,J) = XM(I,J)
   20    CONTINUE
         XM1(I,3) = XM1(I,3) * 1.0E+8
   30 CONTINUE
C
      WRITE (IPRT,1290)
      WRITE (IPRT,1000)
      CALL LLS(Y1, XM, N, IXM, NPAR, RES, LDSTAK)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1290)
      WRITE (IPRT,1000)
      CALL LLS(Y, XM1, N, IXM, NPAR, RES, LDSTAK)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1290)
      WRITE (IPRT,1000)
      CALL LLS(Y1, XM1, N, IXM, NPAR, RES, LDSTAK)
      WRITE (IPRT,1500) IERR
C
C     MINIMUM AMOUNT OF WORK AREA.
C
      CALL LDSCMP(15, 0, 0, 0, 0, 0, 'S',
     +            6*N + NPAR*(N+2*NPAR+5) + 1, LDSMIN)
C
      WRITE (IPRT,1300)
      WRITE (IPRT,1000)
      CALL LLS(Y, XM, N, IXM, NPAR, RES, LDSMIN)
      WRITE (IPRT,1500) IERR
      WRITE (IPRT,1430) (RES(I), I = 1, N)
C
C     CHECK RESULTS FOR WEIGHTED ANALYSIS
C
      NPRT = 1111
      CALL SETRV(WT, N, 100.0E0)
      WRITE (IPRT,1310)
      WRITE (IPRT,1030)
      CALL LLSWS(Y, WT, XM, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV,SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
C
      WT(1) = 0.0E0
      WT(5) = 0.0E0
      WT(9) = 0.0E0
C
      WRITE (IPRT,1310)
      WRITE (IPRT,1030)
      CALL LLSWS(Y, WT, XM, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV,SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
      CALL SETRV(WT, N, 100.0E0)
C
      CALL GENR(WT, N, 1.0E0, 1.0E0)
      WRITE (IPRT,1310)
      WRITE (IPRT,1030)
      CALL LLSWS(Y, WT, XM, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV,SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
      CALL SETRV(WT, N, 100.0E0)
C
C     CHECK PRINT CONTROL
C
      NPRT = 1000
      WRITE (IPRT,1320) NPRT
      WRITE (IPRT,1010)
      CALL LLSS(Y, XM, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
C
      NPRT = 2000
      WRITE (IPRT,1320) NPRT
      WRITE (IPRT,1010)
      CALL LLSS(Y, XM, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
C
      NPRT = 200
      WRITE (IPRT,1320) NPRT
      WRITE (IPRT,1010)
      CALL LLSS(Y, XM, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
C
      NPRT = 20
      WRITE (IPRT,1320) NPRT
      WRITE (IPRT,1010)
      CALL LLSS(Y, XM, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
C
      NPRT = 2
      WRITE (IPRT,1320) NPRT
      WRITE (IPRT,1010)
      CALL LLSS(Y, XM, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
C
      NPRT = 0
      WRITE (IPRT,1320) NPRT
      WRITE (IPRT,1010)
      CALL LLSS(Y, XM, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
C
C     CHECK RESULTS FOR N = 2, NPAR = ID+1 = 1
C
      NPRT = 2222
      N = 2
      NPAR = 1
      NDEG = 0
      WRITE (IPRT,1330)
      WRITE (IPRT,1010)
      CALL LLSS(Y, XM, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
      WRITE (IPRT,1330)
      WRITE (IPRT,1070)
      CALL LLSPWS(Y, WT, X, N, NDEG, RES, LDSTAK,
     +   NPRT, LPAR, PAR, NPAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
C
C     CHECK RESULTS FOR N = 1, NPAR = ID+1 = 1
C
      NPRT = 2222
      N = 1
      NPAR = 1
      NDEG = 0
      WRITE (IPRT,1330)
      WRITE (IPRT,1010)
      CALL LLSS(Y, XM, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
      WRITE (IPRT,1330)
      WRITE (IPRT,1070)
      CALL LLSPWS(Y, WT, X, N, NDEG, RES, LDSTAK,
     +   NPRT, LPAR, PAR, NPAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
      N = 9
      NPAR = 4
      NDEG = 3
C
C
C     ILL-CONDITIONED
C
      DO 40 I = 1, 50
         TERM = 1.0E0
         SUM = 0.0E0
         DO 35 J = 1, 6
            XM1(I,J) = TERM
            SUM = SUM + TERM
            TERM = (I-1)*TERM
   35    CONTINUE
         Y1(I) = SUM
   40 CONTINUE
C
      N = 21
      NPAR = 6
      NDEG = 5
      WRITE (IPRT,1340)
      WRITE (IPRT,1010)
      CALL LLSS(Y1, XM1, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
      WRITE (IPRT,1340)
      WRITE (IPRT,1050)
      CALL LLSPS(Y1, XM1(1,2), N, NDEG, RES, LDSTAK,
     +   NPRT, LPAR, PAR, NPAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
C
      N = 50
      NPAR = 6
      NDEG = 5
      CALL NRAND(RAND, 1, 223)
      DO 50 I = 1, N
         CALL NRAND(RAND, 1, 0)
         Y1(I) = Y1(I) + RAND(1)
   50 CONTINUE
      WRITE (IPRT,1340)
      WRITE (IPRT,1010)
      CALL LLSS(Y1, XM1, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
      WRITE (IPRT,1340)
      WRITE (IPRT,1050)
      CALL LLSPS(Y1, XM1(1,2), N, NDEG, RES, LDSTAK,
     +   NPRT, LPAR, PAR, NPAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
C
      NPRT = 1000
      WRITE (IPRT,1340)
      WRITE (IPRT,1010)
      CALL LLSS(Y1, XM1, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
      WRITE (IPRT,1340)
      WRITE (IPRT,1050)
      CALL LLSPS(Y1, XM1(1,2), N, NDEG, RES, LDSTAK,
     +   NPRT, LPAR, PAR, NPAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
C
      N = 45
      CALL SETRV(WT, N, 1.0E0)
      WRITE (IPRT,1340)
      WRITE (IPRT,1010)
      CALL LLSWS(Y1, WT, XM1, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
C
      N = 44
      WRITE (IPRT,1340)
      WRITE (IPRT,1010)
      CALL LLSWS(Y1, WT, XM1, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
C
      N = 41
      WRITE (IPRT,1340)
      WRITE (IPRT,1010)
      CALL LLSWS(Y1, WT, XM1, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
C
      N = 40
      WRITE (IPRT,1340)
      WRITE (IPRT,1010)
      CALL LLSWS(Y1, WT, XM1, N, IXM, NPAR, RES, LDSTAK,
     +   NPRT, PAR, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1500) IERR
      CALL FITSXP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV, RSD)
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (' CALL TO LLS   ')
 1010 FORMAT (' CALL TO LLSS  ')
 1020 FORMAT (' CALL TO LLSW  ')
 1030 FORMAT (' CALL TO LLSWS ')
 1040 FORMAT (' CALL TO LLSP  ')
 1050 FORMAT (' CALL TO LLSPS ')
 1060 FORMAT (' CALL TO LLSPW ')
 1070 FORMAT (' CALL TO LLSPWS')
 1200 FORMAT ('1MISCELLANEOUS ERRORS  -  TEST 1')
 1220 FORMAT ('1MISCELLANEOUS ERRORS  -  TEST 2')
 1240 FORMAT ('1NEGATIVE WEIGHTS')
 1250 FORMAT ('1TOO FEW POSITIVE WEIGHTS')
 1260 FORMAT ('1VALID PROBLEM')
 1270 FORMAT ('1ZERO RESIDUAL PROBLEM')
 1280 FORMAT ('1RANK DEFICIENT PROBLEM')
 1290 FORMAT ('1POORLY SCALED PROBLEM')
 1300 FORMAT ('1MINIMUM WORK AREA SIZE')
 1310 FORMAT ('1WEIGHTED ANALYSIS')
 1320 FORMAT ('1CHECK PRINT CONTROL  -  NPRT = ', I5)
 1330 FORMAT ('1CHECK MINIMUM PROBLEM SIZE')
 1340 FORMAT ('1ILL-CONDITIONED PROBLEM')
 1430 FORMAT (//4H RES/ (1X, E22.14))
 1500 FORMAT (/' IERR = ', I5)
C
      END
