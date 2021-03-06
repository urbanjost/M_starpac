*XAIMD
      SUBROUTINE XAIMD(LDSTAK)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     DEMONSTRATE THE USER CALLABLE ROUTINES IN THE ARIMA FAMILY.
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  SEPTEMBER 1, 1985
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   LDSTAK
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
     +   DELTA,RSD,STOPP,STOPSS
      INTEGER
     +   I,IFCST,IPRT,IVAPRX,IVCV,MIT,MXFAC,MXFC,MXFCO,MXN,MXPAR,N,
     +   NFAC,NFCST,NFCSTO,NPAR,NPARE,NPRT,NTEST
C
C  LOCAL ARRAYS
      DOUBLE PRECISION
     +   FCST(200,10),FCSTSD(200,10),PAR(10),PV(200),RES(200),
     +   SCALE(10),SDPV(200),SDRES(200),STP(10),VCV(10,10),Y(200)
      INTEGER
     +   IFCSTO(10),IFIXED(10),MSPEC(4,10)
C
C  EXTERNAL SUBROUTINES
      EXTERNAL AIME,AIMEC,AIMES,AIMF,AIMFS,AIMX1,FITXSP,IPRINT
C
C  INTRINSIC FUNCTIONS
      INTRINSIC LOG
C
C  COMMON BLOCKS
      COMMON /CSTAK/DSTAK
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION DELTA
C        THE MAXIMUM CHANGE ALLOWED IN THE MODEL PARAMETERS AT THE
C        FIRST ITERATION.
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     DOUBLE PRECISION FCST(200, 10)
C        THE FORECASTS.
C     DOUBLE PRECISION FCSTSD(200, 10)
C        THE STANDARD DEVIATIONS OF THE FORECASTS.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .GE. 1, ERRORS WERE DETECTED.
C     INTEGER IFCST
C        THE FIRST DIMENSION OF THE ARRAY FCST.
C     INTEGER IFCSTO(10)
C        THE INDICES OF THE ORIGINS FOR THE FORECASTS.
C     INTEGER IFIXED(10)
C        THE INDICATOR VALUES USED TO DESIGNATE WHETHER THE
C        PARAMETERS ARE TO BE OPTIMIZED OR ARE TO BE HELD FIXED.  IF
C        IFIXED(I).NE.0, THEN PAR(I) WILL BE OPTIMIZED.  IF
C        IFIXED(I).EQ.0, THEN PAR(I) WILL BE HELD FIXED.
C        IFIXED(I).LT.0, THEN ALL PAR(I),I=1,NPAR, WILL BE OPTIMIZED..
C     INTEGER IPRT
C        THE UNIT NUMBER FOR PRINTED OUTPUT.
C     INTEGER IVAPRX
C        AN INDICATOR VALUE USED TO DESIGNATE WHICH OPTION IS TO BE USED
C        TO COMPUTE THE VARIANCE COVARIANCE MATRIX (VCV), WHERE
C        IVAPRX LE 0 INDICATES THE THE DEFAULT OPTION WILL BE USED
C        IVAPRX EQ 1 INDICATES THE VCV IS TO BE COMPUTED BY
C                       INVERSE(HESSIAN)*TRANSPOSE(JACOBIAN)*JACOBIAN
C                          *INVERSE(HESSIAN)
C                    USING BOTH THE MODEL SUBROUTINE THE USER SUPPLIED
C                    DERIVATIVE SUBROUTINE WHEN IT IS AVAILABLE
C        IVAPRX EQ 2 INDICATES THE VCV IS TO BE COMPUTED BY
C                       INVERSE(HESSIAN)
C                    USING BOTH THE MODEL SUBROUTINE THE USER SUPPLIED
C                    DERIVATIVE SUBROUTINE WHEN IT IS AVAILABLE
C        IVAPRX EQ 3 INDICATES THE VCV IS TO BE COMPUTED BY
C                       INVERSE(TRANSPOSE(JACOBIAN)*JACOBIAN)
C                    USING BOTH THE MODEL SUBROUTINE THE USER SUPPLIED
C                    DERIVATIVE SUBROUTINE WHEN IT IS AVAILABLE
C        IVAPRX EQ 4 INDICATES THE VCV IS TO BE COMPUTED BY
C                       INVERSE(HESSIAN)*TRANSPOSE(JACOBIAN)*JACOBIAN
C                          *INVERSE(HESSIAN)
C                    USING ONLY THE MODEL SUBROUTINE
C        IVAPRX EQ 5 INDICATES THE VCV IS TO BE COMPUTED BY
C                       INVERSE(HESSIAN)
C                    USING ONLY THE MODEL SUBROUTINE
C        IVAPRX EQ 6 INDICATES THE VCV IS TO BE COMPUTED BY
C                       INVERSE(TRANSPOSE(JACOBIAN)*JACOBIAN)
C                    USING ONLY THE MODEL SUBROUTINE
C        IVAPRX GE 7 INDICATES THE DEFAULT OPTION WILL BE USED
C     INTEGER IVCV
C        THE FIRST DIMENSION OF THE VARIANCE COVARIANCE MATRIX VCV.
C     INTEGER LDSTAK
C        THE LENGTH OF THE ARRAY DSTAK.
C     INTEGER MIT
C        THE MAXIMUM NUMBER OF ITERATIONS ALLOWED.
C     INTEGER MSPEC(4,10)
C        THE ARRAY CONTAINING THE VALUES OF P, D, Q, AND S FOR EACH
C        FACTOR.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS.
C     INTEGER NFAC
C        THE NUMBER OF FACTORS IN THE MODEL
C     INTEGER NFCST
C        THE NUMBER OF FORECASTS.
C     INTEGER NFCSTO
C        THE NUMBER OF THE ORIGINS.
C     INTEGER NPAR
C        THE NUMBER OF UNKNOWN PARAMETERS IN THE MODEL.
C     INTEGER NPARE
C        THE NUMBER OF PARAMETERS ESTIMATED BY THE ROUTINE.
C     INTEGER NPRT
C        THE PARAMETER USED TO INDICATE HOW MUCH PRINTED OUTPUT IS
C        TO BE PROVIDED.
C     INTEGER NTEST
C        THE NUMBER OF THE CURRENT TEST.
C     DOUBLE PRECISION PAR(10)
C        THE ARRAY IN WHICH THE CURRENT ESTIMATES OF THE UNKNOWN
C        PARAMETERS ARE STORED.
C     DOUBLE PRECISION PV(200)
C        THE PREDICTED VALUE BASED ON THE CURRENT PARAMETER ESTIMATES
C     DOUBLE PRECISION RES(200)
C        THE RESIDUALS FROM THE FIT.
C     DOUBLE PRECISION RSD
C        THE VALUE OF THE RESIDUAL STANDARD DEVIATION AT THE SOLUTION.
C     DOUBLE PRECISION SCALE(10)
C        A VALUE TO INDICATE USE OF THE DEFAULT VALUES OF
C        THE TYPICAL SIZE OF THE UNKNOWN PARAMETERS.
C     DOUBLE PRECISION SDPV(200)
C        THE STANDARD DEVIATION OF THE PREDICTED VALUE.
C     DOUBLE PRECISION SDRES(200)
C        THE STANDARD DEVIATIONS OF THE RESIDUALS.
C     DOUBLE PRECISION STOPP
C        THE STOPPING CRITERION FOR THE TEST BASED ON THE MAXIMUM SCALED
C        RELATIVE CHANGE IN THE ELEMENTS OF THE MODEL PARAMETER VECTOR
C     DOUBLE PRECISION STOPSS
C        THE STOPPING CRITERION FOR THE TEST BASED ON THE RATIO OF THE
C        PREDICTED DECREASE IN THE RESIDUAL SUM OF SQUARES (COMPUTED
C        BY STARPAC) TO THE CURRENT RESIDUAL SUM OF SQUARES ESTIMATE.
C     DOUBLE PRECISION STP(10)
C        THE RCSTEP SIZE ARRAY.
C     DOUBLE PRECISION VCV(10,10)
C        THE COVARIANCE MATRIX.
C     DOUBLE PRECISION Y(200)
C        THE ARRAY OF THE DEPENDENT VARIABLE.
C
C     DEFINE CONSTANTS
C
      DATA    Y(  1),   Y(  2),   Y(  3),   Y(  4),   Y(  5),   Y(  6)
     +    / 112.0D0, 118.0D0, 132.0D0, 129.0D0, 121.0D0, 135.0D0/
      DATA    Y(  7),   Y(  8),   Y(  9),   Y( 10),   Y( 11),   Y( 12)
     +    / 148.0D0, 148.0D0, 136.0D0, 119.0D0, 104.0D0, 118.0D0/
      DATA    Y( 13),   Y( 14),   Y( 15),   Y( 16),   Y( 17),   Y( 18)
     +    / 115.0D0, 126.0D0, 141.0D0, 135.0D0, 125.0D0, 149.0D0/
      DATA    Y( 19),   Y( 20),   Y( 21),   Y( 22),   Y( 23),   Y( 24)
     +    / 170.0D0, 170.0D0, 158.0D0, 133.0D0, 114.0D0, 140.0D0/
      DATA    Y( 25),   Y( 26),   Y( 27),   Y( 28),   Y( 29),   Y( 30)
     +    / 145.0D0, 150.0D0, 178.0D0, 163.0D0, 172.0D0, 178.0D0/
      DATA    Y( 31),   Y( 32),   Y( 33),   Y( 34),   Y( 35),   Y( 36)
     +    / 199.0D0, 199.0D0, 184.0D0, 162.0D0, 146.0D0, 166.0D0/
      DATA    Y( 37),   Y( 38),   Y( 39),   Y( 40),   Y( 41),   Y( 42)
     +    / 171.0D0, 180.0D0, 193.0D0, 181.0D0, 183.0D0, 218.0D0/
      DATA    Y( 43),   Y( 44),   Y( 45),   Y( 46),   Y( 47),   Y( 48)
     +    / 230.0D0, 242.0D0, 209.0D0, 191.0D0, 172.0D0, 194.0D0/
      DATA    Y( 49),   Y( 50),   Y( 51),   Y( 52),   Y( 53),   Y( 54)
     +    / 196.0D0, 196.0D0, 236.0D0, 235.0D0, 229.0D0, 243.0D0/
      DATA    Y( 55),   Y( 56),   Y( 57),   Y( 58),   Y( 59),   Y( 60)
     +    / 264.0D0, 272.0D0, 237.0D0, 211.0D0, 180.0D0, 201.0D0/
      DATA    Y( 61),   Y( 62),   Y( 63),   Y( 64),   Y( 65),   Y( 66)
     +    / 204.0D0, 188.0D0, 235.0D0, 227.0D0, 234.0D0, 264.0D0/
      DATA    Y( 67),   Y( 68),   Y( 69),   Y( 70),   Y( 71),   Y( 72)
     +    / 302.0D0, 293.0D0, 259.0D0, 229.0D0, 203.0D0, 229.0D0/
      DATA    Y( 73),   Y( 74),   Y( 75),   Y( 76),   Y( 77),   Y( 78)
     +    / 242.0D0, 233.0D0, 267.0D0, 269.0D0, 270.0D0, 315.0D0/
      DATA    Y( 79),   Y( 80),   Y( 81),   Y( 82),   Y( 83),   Y( 84)
     +    / 364.0D0, 347.0D0, 312.0D0, 274.0D0, 237.0D0, 278.0D0/
      DATA    Y( 85),   Y( 86),   Y( 87),   Y( 88),   Y( 89),   Y( 90)
     +    / 284.0D0, 277.0D0, 317.0D0, 313.0D0, 318.0D0, 374.0D0/
      DATA    Y( 91),   Y( 92),   Y( 93),   Y( 94),   Y( 95),   Y( 96)
     +    / 413.0D0, 405.0D0, 355.0D0, 306.0D0, 271.0D0, 306.0D0/
      DATA    Y( 97),   Y( 98),   Y( 99),   Y(100),   Y(101),   Y(102)
     +    / 315.0D0, 301.0D0, 356.0D0, 348.0D0, 355.0D0, 422.0D0/
      DATA    Y(103),   Y(104),   Y(105),   Y(106),   Y(107),   Y(108)
     +    / 465.0D0, 467.0D0, 404.0D0, 347.0D0, 305.0D0, 336.0D0/
      DATA    Y(109),   Y(110),   Y(111),   Y(112),   Y(113),   Y(114)
     +    / 340.0D0, 318.0D0, 362.0D0, 348.0D0, 363.0D0, 435.0D0/
      DATA    Y(115),   Y(116),   Y(117),   Y(118),   Y(119),   Y(120)
     +    / 491.0D0, 505.0D0, 404.0D0, 359.0D0, 310.0D0, 337.0D0/
      DATA    Y(121),   Y(122),   Y(123),   Y(124),   Y(125),   Y(126)
     +    / 360.0D0, 342.0D0, 406.0D0, 396.0D0, 420.0D0, 472.0D0/
      DATA    Y(127),   Y(128),   Y(129),   Y(130),   Y(131),   Y(132)
     +    / 548.0D0, 559.0D0, 463.0D0, 407.0D0, 362.0D0, 405.0D0/
      DATA    Y(133),   Y(134),   Y(135),   Y(136),   Y(137),   Y(138)
     +    / 417.0D0, 391.0D0, 419.0D0, 461.0D0, 472.0D0, 535.0D0/
      DATA    Y(139),   Y(140),   Y(141),   Y(142),   Y(143),   Y(144)
     +    / 622.0D0, 606.0D0, 508.0D0, 461.0D0, 390.0D0, 432.0D0/
C
      CALL IPRINT(IPRT)
C
      DO 10 I = 1, 144
        Y(I) = LOG(Y(I))
   10 CONTINUE
C
C     SET DIMENSIONS
C
      MXN = 200
      MXPAR = 10
      MXFC = 200
      MXFCO = 10
      MXFAC = 10
C
C
      NTEST = 0
C
C
C     **TEST ON NORMAL STATEMENT**
C
C
      NTEST = NTEST + 1
      WRITE (IPRT,1330) NTEST
      WRITE (IPRT,1130)
      WRITE (IPRT,1000)
      CALL AIMX1(MXN, MXPAR, MXFC, MXFCO, MXFAC,
     +   1, N, MSPEC, NFAC, PAR, NPAR, RES,
     +   IFIXED, STP, MIT, STOPSS, STOPP, SCALE, DELTA, IVAPRX, NPRT,
     +   NPARE, RSD, PV, SDPV, SDRES, VCV, IVCV,
     +   NFCST, NFCSTO, IFCSTO, FCST, IFCST, FCSTSD)
      CALL AIME (Y, N, MSPEC, NFAC, PAR, NPAR, RES, LDSTAK)
      WRITE (IPRT,1120) IERR
      CALL FITXSP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV,
     +   N, NPARE, RSD)
C
      NTEST = NTEST + 1
      WRITE (IPRT,1330) NTEST
      WRITE (IPRT,1130)
      WRITE (IPRT,1010)
      WRITE (IPRT,1340) IFIXED(1), STP(1), MIT, STOPSS, STOPP,
     +   SCALE(1), DELTA, IVAPRX, NPRT
      CALL AIMX1(MXN, MXPAR, MXFC, MXFCO, MXFAC,
     +   1, N, MSPEC, NFAC, PAR, NPAR, RES,
     +   IFIXED, STP, MIT, STOPSS, STOPP, SCALE, DELTA, IVAPRX, NPRT,
     +   NPARE, RSD, PV, SDPV, SDRES, VCV, IVCV,
     +   NFCST, NFCSTO, IFCSTO, FCST, IFCST, FCSTSD)
      CALL AIMEC(Y, N, MSPEC, NFAC, PAR, NPAR, RES, LDSTAK,
     +   IFIXED, STP, MIT, STOPSS, STOPP, SCALE, DELTA, IVAPRX, NPRT)
      WRITE (IPRT,1350) IFIXED(1), STP(1), MIT, STOPSS, STOPP,
     +   SCALE(1), DELTA, IVAPRX, NPRT
      WRITE (IPRT,1120) IERR
      CALL FITXSP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV,
     +   N, NPARE, RSD)
C
      NTEST = NTEST + 1
      WRITE (IPRT,1330) NTEST
      WRITE (IPRT,1130)
      WRITE (IPRT,1020)
      WRITE (IPRT,1340) IFIXED(1), STP(1), MIT, STOPSS, STOPP,
     +   SCALE(1), DELTA, IVAPRX, NPRT
      CALL AIMX1(MXN, MXPAR, MXFC, MXFCO, MXFAC,
     +   1, N, MSPEC, NFAC, PAR, NPAR, RES,
     +   IFIXED, STP, MIT, STOPSS, STOPP, SCALE, DELTA, IVAPRX, NPRT,
     +   NPARE, RSD, PV, SDPV, SDRES, VCV, IVCV,
     +   NFCST, NFCSTO, IFCSTO, FCST, IFCST, FCSTSD)
      CALL AIMES(Y, N, MSPEC, NFAC, PAR, NPAR, RES, LDSTAK,
     +   IFIXED, STP, MIT, STOPSS, STOPP, SCALE, DELTA, IVAPRX, NPRT,
     +   NPARE, RSD, PV, SDPV, SDRES, VCV, IVCV)
      WRITE (IPRT,1350) IFIXED(1), STP(1), MIT, STOPSS, STOPP,
     +   SCALE(1), DELTA, IVAPRX, NPRT
      WRITE (IPRT,1120) IERR
      CALL FITXSP(PAR, PV, SDPV, RES, SDRES, VCV, N, NPAR, IVCV,
     +   N, NPARE, RSD)
C
      NTEST = NTEST + 1
      WRITE (IPRT,1330) NTEST
      WRITE (IPRT,1130)
      WRITE (IPRT,1030)
      CALL AIMX1(MXN, MXPAR, MXFC, MXFCO, MXFAC,
     +   1, N, MSPEC, NFAC, PAR, NPAR, RES,
     +   IFIXED, STP, MIT, STOPSS, STOPP, SCALE, DELTA, IVAPRX, NPRT,
     +   NPARE, RSD, PV, SDPV, SDRES, VCV, IVCV,
     +   NFCST, NFCSTO, IFCSTO, FCST, IFCST, FCSTSD)
      CALL AIMF (Y, N, MSPEC, NFAC, PAR, NPAR, LDSTAK)
      WRITE (IPRT,1120) IERR
C
      NTEST = NTEST + 1
      WRITE (IPRT,1330) NTEST
      WRITE (IPRT,1130)
      WRITE (IPRT,1040)
      WRITE (IPRT,1340) IFIXED(1), STP(1), MIT, STOPSS, STOPP,
     +   SCALE(1), DELTA, IVAPRX, NPRT
      CALL AIMX1(MXN, MXPAR, MXFC, MXFCO, MXFAC,
     +   1, N, MSPEC, NFAC, PAR, NPAR, RES,
     +   IFIXED, STP, MIT, STOPSS, STOPP, SCALE, DELTA, IVAPRX, NPRT,
     +   NPARE, RSD, PV, SDPV, SDRES, VCV, IVCV,
     +   NFCST, NFCSTO, IFCSTO, FCST, IFCST, FCSTSD)
      CALL AIMFS(Y, N, MSPEC, NFAC, PAR, NPAR, LDSTAK,
     +   NFCST, NFCSTO, IFCSTO, NPRT, FCST, IFCST, FCSTSD)
      WRITE (IPRT,1350) IFIXED(1), STP(1), MIT, STOPSS, STOPP,
     +   SCALE(1), DELTA, IVAPRX, NPRT
      WRITE (IPRT,1120) IERR
      CALL FITXSP(PAR, FCST(1,1), FCST(1,2), FCST(1,3), FCSTSD, VCV,
     +  N, NPAR, IVCV, N, NPARE, RSD)
C
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (14H TEST OF AIM  )
 1010 FORMAT (15H TEST OF AIMC  )
 1020 FORMAT (15H TEST OF AIMS  )
 1030 FORMAT (14H TEST OF AIMF )
 1040 FORMAT (15H TEST OF AIMFS )
 1120 FORMAT (/29H ***** RETURNED RESULTS *****, 5X, 15H (-1 INDICATES ,
     +   39HVALUE NOT CHANGED BY CALLED SUBROUTINE)//9H IERR IS , I3)
 1130 FORMAT (15H NORMAL PROBLEM)
 1330 FORMAT ('1ARIMA TEST NUMBER', I5)
 1340 FORMAT (//24H INPUT   -  IFIXED(1) = , I6, 9X, ', STP(1) = ',
     +   G15.8, ',    MIT = ',I5, ', STOPSS = ', G15.8, 10H, STOPP = ,
     +   G15.8/13X, 'SCALE(1) = ', G15.8, ',  DELTA = ', G15.8,
     +   ', IVAPRX = ', I5, ',   NPRT = ', I5//)
 1350 FORMAT (//24H OUTPUT  -  IFIXED(1) = , I6, 9X, ', STP(1) = ',
     +   G15.8, ',    MIT = ',I5, ', STOPSS = ', G15.8, 10H, STOPP = ,
     +   G15.8/13X, 'SCALE(1) = ', G15.8, ',  DELTA = ', G15.8,
     +   ', IVAPRX = ', I5, ',   NPRT = ', I5//)
      END
