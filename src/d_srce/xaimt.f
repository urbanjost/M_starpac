*XAIMT
      SUBROUTINE XAIMT(LDSTAK)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     ROUTINE TO TEST THE TIME SERIES MODEL ESTIMATION ROUTINES.
C
C     SERIES Y IS THE AIRLINE DATA LISTED ON PAGE 531 OF BOX AND JENKINS
C
C     WRITTEN BY - JANET R. DONALDSON
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  DECEMBER 2, 1985
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
     +   I,IFCST,IPRT,IVAPRX,IVCV,MIT,NFAC,NPAR,NPARE,NPRT,NY
C
C  LOCAL ARRAYS
      DOUBLE PRECISION
     +   FCST(50,5),FCSTSD(50),PAR(50),PV(200),RES(200),SCALE(50),
     +   SDPV(200),SDRES(200),STP(50),VCV(10,10),Y(200),YLOG(200),
     +   YT(200)
      INTEGER
     +   IFIXED(50),MSPEC(4,50)
C
C  EXTERNAL SUBROUTINES
      EXTERNAL AIME,AIMEC,AIMES,AIMF,AIMFS,IPRINT
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
C     DOUBLE PRECISION FCST(50,5)
C        THE FORECASTS.
C     DOUBLE PRECISION FCSTSD(50)
C        THE STANDARD DEVIATIONS OF THE FORECASTS.
C     INTEGER I
C        *
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .GE. 1, ERRORS WERE DETECTED.
C     INTEGER IFCST
C        THE FIRST DIMENSION OF THE ARRAY FCST.
C     INTEGER IFIXED(50)
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
C     INTEGER MSPEC(4,50)
C        THE ARRAY CONTAINING THE VALUES OF P, D, Q, AND S FOR EACH
C        FACTOR.
C     INTEGER NFAC
C        THE NUMBER OF FACTORS IN THE MODEL
C     INTEGER NPAR
C        THE NUMBER OF UNKNOWN PARAMETERS IN THE MODEL.
C     INTEGER NPARE
C        THE NUMBER OF PARAMETERS ESTIMATED BY THE ROUTINE.
C     INTEGER NPRT
C        THE PARAMETER USED TO INDICATE HOW MUCH PRINTED OUTPUT IS
C        TO BE PROVIDED.
C     INTEGER NY
C        THE NUMBER OF OBSERVATIONS.
C     DOUBLE PRECISION PAR(50)
C        THE ARRAY IN WHICH THE CURRENT ESTIMATES OF THE UNKNOWN
C        PARAMETERS ARE STORED.
C     DOUBLE PRECISION PV(200)
C        THE PREDICTED VALUE BASED ON THE CURRENT PARAMETER ESTIMATES
C     DOUBLE PRECISION RES(200)
C        THE RESIDUALS FROM THE FIT.
C     DOUBLE PRECISION RSD
C        THE VALUE OF THE RESIDUAL STANDARD DEVIATION AT THE SOLUTION.
C     DOUBLE PRECISION SCALE(50)
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
C     DOUBLE PRECISION STP(50)
C        THE RCSTEP SIZE ARRAY.
C     DOUBLE PRECISION VCV(10,10)
C        THE COVARIANCE MATRIX.
C     DOUBLE PRECISION Y(200),YLOG(200),YT(200)
C        THE ARRAY OF THE DEPENDENT VARIABLE.
C
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
C     COMMENCE BODY OF ROUTINE
C
C     TEST AGAINST PUBLISHED RESULTS
C
      NY = 144
      DO 10 I = 1, NY
         YLOG(I) = LOG(Y(I))
   10 CONTINUE
C
      NFAC = 2
      MSPEC(1,1) = 0
      MSPEC(2,1) = 1
      MSPEC(3,1) = 1
      MSPEC(4,1) = 1
C
      MSPEC(1,2) = 0
      MSPEC(2,2) = 1
      MSPEC(3,2) = 1
      MSPEC(4,2) = 12
C
      NPAR = 3
      PAR(1) = 0.0D0
      PAR(2) = 0.40D0
      PAR(3) = 0.60D0
C
      IFIXED(1) = 1
      IFIXED(2) = 0
      IFIXED(3) = 0
C
      STOPSS = -1.0D0
      STOPP = -1.0D0
      SCALE(1) = -1.0D0
      SCALE(2) = 1.0D-7
      SCALE(3) = 1.0D-7
      STP(1) = -1.0D0
      STP(2) = 1.0D-7
      STP(3) = 1.0D-7
      MIT = 0
      NPRT = -1
      DELTA = -1.0D0
      IVAPRX = -1
C
      WRITE(IPRT, 1000)
      CALL AIMEC (YLOG, NY, MSPEC, NFAC,
     +   PAR, NPAR, RES,LDSTAK, IFIXED,STP, MIT, STOPSS, STOPP,
     +   SCALE, DELTA, IVAPRX, NPRT)
C
      WRITE (IPRT, 1005)
      PAR(1) = 0.0D0
      PAR(2) = 0.395D0
      PAR(3) = 0.615D0
      CALL AIMFS (YLOG, NY, MSPEC, NFAC,
     +   PAR, NPAR, LDSTAK, NY/10+1, 1, NY, NPRT, FCST, 50, FCSTSD)
C
      SCALE(1) = 1.0D-7
      SCALE(2) = 1.0D-7
      SCALE(3) = 1.0D-7
C
      NFAC = 2
      MSPEC(1,1) = 0
      MSPEC(2,1) = 1
      MSPEC(3,1) = 1
      MSPEC(4,1) = 1
C
      MSPEC(1,2) = 0
      MSPEC(2,2) = 0
      MSPEC(3,2) = 1
      MSPEC(4,2) = 12
C
      WRITE (IPRT, 1000)
      CALL AIMEC (YLOG, NY, MSPEC, NFAC,
     +   PAR, NPAR, RES,LDSTAK, IFIXED,STP, MIT, STOPSS, STOPP,
     +   SCALE, DELTA, IVAPRX, NPRT)
      NY = 20
      WRITE (IPRT, 1000)
      CALL AIMEC (YLOG, NY, MSPEC, NFAC,
     +   PAR, NPAR, RES,LDSTAK, IFIXED,STP, MIT, STOPSS, STOPP,
     +   SCALE, DELTA, IVAPRX, NPRT)
C
      NFAC = 2
      MSPEC(1,1) = 0
      MSPEC(2,1) = 0
      MSPEC(3,1) = 1
      MSPEC(4,1) = 1
C
      MSPEC(1,2) = 0
      MSPEC(2,2) = 0
      MSPEC(3,2) = 1
      MSPEC(4,2) = 12
C
      NY = 144
      WRITE (IPRT, 1000)
      CALL AIMEC (YLOG, NY, MSPEC, NFAC,
     +   PAR, NPAR, RES,LDSTAK, IFIXED,STP, MIT, STOPSS, STOPP,
     +   SCALE, DELTA, IVAPRX, NPRT)
C
C     EXAMPLE FROM PAGE 212 OF BOX AND JENKINS (1970)
C     N.B. ADD PRINT STATEMENTS TO MDLTS2 TO CHECK COMPUTATIONS
C          AT FIRST CALL AGAINST THOSE LISTED ON PAGE 214.
C
      WRITE(IPRT, 1000)
      NY = 10
      YT(1) = 460.0D0
      YT(2) = 457.0D0
      YT(3) = 452.0D0
      YT(4) = 459.0D0
      YT(5) = 462.0D0
      YT(6) = 459.0D0
      YT(7) = 463.0D0
      YT(8) = 479.0D0
      YT(9) = 493.0D0
      YT(10) = 490.0D0
C
      NFAC = 1
      MSPEC(1,1) = 0
      MSPEC(2,1) = 1
      MSPEC(3,1) = 1
      MSPEC(4,1) = 1
C
      NPAR = 2
      PAR(1) = 0.0D0
      PAR(2) = 0.5D0
C
      IFIXED(1) = 1
      IFIXED(2) = 0
C
      CALL AIMEC (YT, NY, MSPEC, NFAC,
     +   PAR, NPAR, RES,LDSTAK, IFIXED,STP, MIT, STOPSS, STOPP,
     +   SCALE, DELTA, IVAPRX, NPRT)
C
C     EXAMPLE FROM PAGE 216 OF BOX AND JENKINS (1970)
C     N.B. ADD PRINT STATEMENTS TO MDLTS2 TO CHECK COMPUTATIONS
C          AT FIRST CALL AGAINST THOSE LISTED ON PAGE 218.
C
      WRITE(IPRT, 1000)
      NY = 12
      YT(1) = 2.0D0
      YT(2) = 0.8D0
      YT(3) = -0.3D0
      YT(4) = -0.3D0
      YT(5) = -1.9D0
      YT(6) = 0.3D0
      YT(7) = 3.2D0
      YT(8) = 1.6D0
      YT(9) = -0.7D0
      YT(10) = 3.0D0
      YT(11) = 4.3D0
      YT(12) = 1.1D0
C
      NFAC = 1
      MSPEC(1,1) = 1
      MSPEC(2,1) = 0
      MSPEC(3,1) = 1
      MSPEC(4,1) = 1
C
      NPAR = 3
      PAR(1) = 0.3D0
      PAR(2) = 0.0D0
      PAR(3) = 0.7D0
C
      IFIXED(1) = 0
      IFIXED(2) = 1
      IFIXED(3) = 0
C
      CALL AIMEC (YT, NY, MSPEC, NFAC,
     +   PAR, NPAR, RES, LDSTAK, IFIXED, STP, MIT, STOPSS, STOPP,
     +   SCALE, DELTA, IVAPRX, NPRT)
C
C     TEST ERROR MESSAGES
C
      WRITE (IPRT, 1010)
      NY = 0
      NFAC = 0
      CALL AIME (YT, NY, MSPEC, NFAC,
     +   PAR, NPAR, RES, LDSTAK)
      CALL AIMEC (YT, NY, MSPEC, NFAC,
     +   PAR, NPAR, RES, LDSTAK, IFIXED, STP, MIT, STOPSS, STOPP,
     +   SCALE, DELTA, IVAPRX, NPRT)
      CALL AIMES (YT, NY, MSPEC, NFAC,
     +   PAR, NPAR, RES, LDSTAK, IFIXED, STP, MIT, STOPSS, STOPP,
     +   SCALE, DELTA, IVAPRX, NPRT, NPARE, RSD, PV, SDPV, SDRES, VCV,
     +   IVCV)
      CALL AIMF (Y, NY, MSPEC, NFAC, PAR, NPAR, LDSTAK)
      CALL AIMFS (Y, NY, MSPEC, NFAC,
     +   PAR, NPAR, LDSTAK, NY/10+1, 1, NY, NPRT, FCST, 50, FCSTSD)
C
      NY = 144
      NFAC = 2
      MSPEC(1,1) = -1
      CALL AIME (YT, NY, MSPEC, NFAC,
     +   PAR, NPAR, RES, LDSTAK)
      CALL AIMEC (YT, NY, MSPEC, NFAC,
     +   PAR, NPAR, RES, LDSTAK, IFIXED, STP, MIT, STOPSS, STOPP,
     +   SCALE, DELTA, IVAPRX, NPRT)
      CALL AIMES (YT, NY, MSPEC, NFAC,
     +   PAR, NPAR, RES, LDSTAK, IFIXED, STP, MIT, STOPSS, STOPP,
     +   SCALE, DELTA, IVAPRX, NPRT, NPARE, RSD, PV, SDPV, SDRES, VCV,
     +   IVCV)
      CALL AIMF (Y, NY, MSPEC, NFAC, PAR, NPAR, LDSTAK)
      CALL AIMFS (Y, NY, MSPEC, NFAC,
     +   PAR, NPAR, LDSTAK, NY/10+1, 1, NY, NPRT, FCST, 50, FCSTSD)
      NY = 144
      NFAC = 2
      MSPEC(1,1) = 0
      NPAR = 1
      CALL AIME (YT, NY, MSPEC, NFAC,
     +   PAR, NPAR, RES, LDSTAK)
      CALL AIMEC (YT, NY, MSPEC, NFAC,
     +   PAR, NPAR, RES, LDSTAK, IFIXED, STP, MIT, STOPSS, STOPP,
     +   SCALE, DELTA, IVAPRX, NPRT)
      CALL AIMES (YT, NY, MSPEC, NFAC,
     +   PAR, NPAR, RES, LDSTAK, IFIXED, STP, MIT, STOPSS, STOPP,
     +   SCALE, DELTA, IVAPRX, NPRT, NPARE, RSD, PV, SDPV, SDRES, VCV,
     +   IVCV)
      CALL AIMF (Y, NY, MSPEC, NFAC, PAR, NPAR, LDSTAK)
      CALL AIMFS (Y, NY, MSPEC, NFAC,
     +   PAR, NPAR, LDSTAK, NY/10+1, 1, NY, NPRT, FCST, 50, FCSTSD)
      NY = 144
      NFAC = 2
      MSPEC(1,1) = 0
      NPAR = 3
      DO 20 I = 1, NPAR
        IFIXED(I) = 1
   20 CONTINUE
      IVCV = 0
      IFCST = 0
      CALL AIMEC (YT, NY, MSPEC, NFAC,
     +   PAR, NPAR, RES, LDSTAK, IFIXED, STP, MIT, STOPSS, STOPP,
     +   SCALE, DELTA, IVAPRX, NPRT)
      CALL AIMES (YT, NY, MSPEC, NFAC,
     +   PAR, NPAR, RES, LDSTAK, IFIXED, STP, MIT, STOPSS, STOPP,
     +   SCALE, DELTA, IVAPRX, NPRT, NPARE, RSD, PV, SDPV, SDRES, VCV,
     +   IVCV)
      CALL AIMFS (Y, NY, MSPEC, NFAC,
     +   PAR, NPAR, LDSTAK, NY/10+1, 1, NY, NPRT, FCST, IFCST, FCSTSD)
      DO 30 I = 1, NPAR
        IFIXED(I) = 1
   30 CONTINUE
      IVCV = 0
      STP(2) = -1.0D0
      SCALE(2) = -1.0D0
      CALL AIMEC (YT, NY, MSPEC, NFAC,
     +   PAR, NPAR, RES, LDSTAK, IFIXED, STP, MIT, STOPSS, STOPP,
     +   SCALE, DELTA, IVAPRX, NPRT)
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT ('1TEST OF ARIMA ESTIMATION ROUTINES')
 1005 FORMAT ('1TEST OF ARIMA FORECASTING ROUTINES')
 1010 FORMAT ('1TEST OF ERROR CHECKING FACILITIES')
      END
