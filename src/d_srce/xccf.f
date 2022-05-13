*XCCF
      SUBROUTINE XCCF(LDS)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     ROUTINE TO TEST THE TIME SERIES CORRELATION SUBROUTINES
C
C     SERIES Y1 AND Y2 ARE LISTED AS SERIES X1 AND X2 ON PAGE OF 361 OF
C     JENKINS AND WATTS.  CCF FOR SERIES Y1 AND Y2 ARE PLOTTED ON PAGE 3
C     AND LISTED ON PAGE 420.
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
      DOUBLE PRECISION
     +   CMISS,Y1MISS,Y2MISS,YMISS0
      INTEGER
     +   ICCOV,INLPPC,IPRT,ITEST,IYM,IYMFFT,JCCOV,JNLPPC,LAGMAX,
     +   LDSTAK,LYFFT,M,N,NLAG,NPRT,NYD
C
C  LOCAL ARRAYS
      DOUBLE PRECISION
     +   CCOV(30,5,5),Y1(100),Y2(100),YFFT1(150),YFFT2(150),YM(150,5),
     +   YMFFT(150,5),YMMISS(5)
      INTEGER
     +   NLPPC(30,5,5)
C
C  EXTERNAL SUBROUTINES
      EXTERNAL CCF,CCFF,CCFFS,CCFM,CCFMS,CCFS,CCFXP,IPRINT,DCOPY,SETRA,
     +   SETRV
C
C  COMMON BLOCKS
      COMMON /CSTAK/DSTAK
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION CCOV(30,5,5)
C        THE CROSS COVARIANCE ARRAY.
C     DOUBLE PRECISION CMISS
C        THE MISSING VALUE CODE FOR THE RETURNED CCVF ESTIMATES
C        (VECTOR CCOV).
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER ICCOV
C        THE FIRST DIMENSION OF THE ARRAY CCOV.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED
C     INTEGER INLPPC
C        THE FIRST DIMENSION OF THE ARRAY NLPPC.
C     INTEGER IPRT
C        THE UNIT NUMBER FOR PRINTED OUTPUT.
C     INTEGER ITEST
C        THE NUMBER OF THE TEST BEING PERFORMED
C     INTEGER IYM, IYMFFT
C        THE FIRST DIMENSION OF THE ARRAYS YM AND YMFFT, RESPECTIVELY.
C     INTEGER JCCOV, JNLPPC
C        THE SECOND DIMENSIONS OF THE ARRAYS CCOV AND NLPPC,
C        RESPECTIVELY.
C     INTEGER LAGMAX
C        THE MAXIMUM LAG VALUE REQUESTED.
C     INTEGER LDS, LDSTAK
C        THE LENGTH OF THE ARRAY DSTAK.
C     INTEGER LYFFT
C        THE LENGTH OF THE ARRAYS USED WHEN THE COMPUTATIONS ARE
C        PERFORMED BY THE FFT.
C     INTEGER M
C        THE NUMBER OF SERIES IN THE MULTIVARIATE TIME SERIES YM.
C     INTEGER N
C        THE INTEGER NUMBER OF OBSERVATIONS IN EACH SERIES
C     INTEGER NLAG
C        THE NUMBER OF LAGS AT WHICH THE ACVF WAS COMPUTED.
C     INTEGER NLPPC(30,5,5)
C        THE ARRAY CONTAINING THE NUMBER OF LAGGED PRODUCT PAIRS
C        USED TO COMPUTE EACH ACVF ESTIMATE.
C     INTEGER NPRT
C        THE INDICATOR VARIABLE USED TO SPECIFY WHETHER OR NOT
C        PRINTED OUTPUT IS TO BE GIVEN, WHERE IF THE VALUE OF
C        NPRT IS ZERO, NO OUTPUT IS MADE.
C     INTEGER NYD
C        THE NUMBER OF OBSERVATIONS IN THE SERIES TO BE DIFFERENCED.
C     DOUBLE PRECISION YFFT1(150), YFFT2(150)
C        THE VECTORS USED FOR STORING THE SERIES FOR THE ROUTINES
C        USING THE FFT.
C     DOUBLE PRECISION YM(150,5), YMFFT(150,5)
C        THE ARRAYS USED FOR MULTIVARIATE TIME SERIES.
C     DOUBLE PRECISION YMISS0, YMMISS(5)
C        THE MISSING VALUE CODES FOR SERIES Y AND YM.
C     DOUBLE PRECISION Y1(100), Y1MISS
C        THE FIRST SERIES, AND ITS MISSING VALUE CODE.
C     DOUBLE PRECISION Y2(100), Y2MISS
C        THE SECOND SERIES, AND ITS MISSING VALUE CODE.
C
C
      DATA   Y1(  1),  Y1(  2),  Y1(  3),  Y1(  4),  Y1(  5),  Y1(  6)
     +    /-0.88D0, -0.16D0, -1.87D0, -1.12D0,  1.38D0,  2.13D0/
      DATA   Y1(  7),  Y1(  8),  Y1(  9),  Y1( 10),  Y1( 11),  Y1( 12)
     +    / 2.76D0,  0.56D0, -0.69D0, -1.79D0, -3.82D0, -2.38D0/
      DATA   Y1( 13),  Y1( 14),  Y1( 15),  Y1( 16),  Y1( 17),  Y1( 18)
     +    / 1.00D0,  0.70D0, -0.15D0,  0.98D0,  0.11D0, -0.35D0/
      DATA   Y1( 19),  Y1( 20),  Y1( 21),  Y1( 22),  Y1( 23),  Y1( 24)
     +    /-0.73D0,  0.89D0, -1.63D0, -0.44D0, -1.37D0, -1.71D0/
      DATA   Y1( 25),  Y1( 26),  Y1( 27),  Y1( 28),  Y1( 29),  Y1( 30)
     +    /-1.22D0, -2.00D0, -0.22D0,  0.38D0,  1.31D0,  0.71D0/
      DATA   Y1( 31),  Y1( 32),  Y1( 33),  Y1( 34),  Y1( 35),  Y1( 36)
     +    / 0.32D0,  0.48D0, -1.88D0, -0.94D0, -1.54D0, -0.13D0/
      DATA   Y1( 37),  Y1( 38),  Y1( 39),  Y1( 40),  Y1( 41),  Y1( 42)
     +    / 1.02D0,  0.02D0, -0.77D0,  0.11D0, -0.60D0, -0.52D0/
      DATA   Y1( 43),  Y1( 44),  Y1( 45),  Y1( 46),  Y1( 47),  Y1( 48)
     +    /-0.09D0,  1.23D0,  1.46D0,  0.61D0,  0.42D0,  2.16D0/
      DATA   Y1( 49),  Y1( 50),  Y1( 51),  Y1( 52),  Y1( 53),  Y1( 54)
     +    / 3.18D0,  2.10D0,  0.37D0, -0.24D0,  0.57D0, -0.53D0/
      DATA   Y1( 55),  Y1( 56),  Y1( 57),  Y1( 58),  Y1( 59),  Y1( 60)
     +    / 2.44D0,  1.02D0, -0.53D0, -2.49D0, -2.12D0, -1.04D0/
      DATA   Y1( 61),  Y1( 62),  Y1( 63),  Y1( 64),  Y1( 65),  Y1( 66)
     +    /-0.12D0, -1.88D0, -1.50D0,  1.54D0,  3.33D0,  3.08D0/
      DATA   Y1( 67),  Y1( 68),  Y1( 69),  Y1( 70),  Y1( 71),  Y1( 72)
     +    / 1.71D0,  0.79D0,  1.55D0,  0.89D0, -0.89D0, -1.18D0/
      DATA   Y1( 73),  Y1( 74),  Y1( 75),  Y1( 76),  Y1( 77),  Y1( 78)
     +    / 0.89D0,  1.71D0,  3.05D0,  0.15D0, -1.04D0,  0.12D0/
      DATA   Y1( 79),  Y1( 80),  Y1( 81),  Y1( 82),  Y1( 83),  Y1( 84)
     +    / 0.08D0,  0.11D0, -2.62D0, -1.28D0,  1.07D0,  3.20D0/
      DATA   Y1( 85),  Y1( 86),  Y1( 87),  Y1( 88),  Y1( 89),  Y1( 90)
     +    / 1.92D0,  0.53D0, -1.08D0,  0.49D0, -0.58D0,  0.17D0/
      DATA   Y1( 91),  Y1( 92),  Y1( 93),  Y1( 94),  Y1( 95),  Y1( 96)
     +    / 1.15D0, -0.97D0, -1.63D0,  1.14D0, -0.67D0, -0.88D0/
      DATA   Y1( 97),  Y1( 98),  Y1( 99),  Y1(100)
     +    /-0.07D0,  0.24D0,  0.55D0, -2.16D0/
      DATA   Y2(  1),  Y2(  2),  Y2(  3),  Y2(  4),  Y2(  5),  Y2(  6)
     +    / 0.79D0,  1.12D0, -1.10D0, -2.39D0, -1.75D0, -0.82D0/
      DATA   Y2(  7),  Y2(  8),  Y2(  9),  Y2( 10),  Y2( 11),  Y2( 12)
     +    /-0.36D0,  1.27D0,  1.75D0,  2.44D0,  0.36D0, -2.10D0/
      DATA   Y2( 13),  Y2( 14),  Y2( 15),  Y2( 16),  Y2( 17),  Y2( 18)
     +    /-1.93D0, -1.30D0, -1.75D0, -0.34D0,  0.74D0,  0.49D0/
      DATA   Y2( 19),  Y2( 20),  Y2( 21),  Y2( 22),  Y2( 23),  Y2( 24)
     +    / 0.70D0,  0.71D0,  0.09D0,  0.59D0,  1.54D0,  0.14D0/
      DATA   Y2( 25),  Y2( 26),  Y2( 27),  Y2( 28),  Y2( 29),  Y2( 30)
     +    / 0.55D0, -1.40D0, -2.55D0, -1.66D0, -0.43D0,  0.58D0/
      DATA   Y2( 31),  Y2( 32),  Y2( 33),  Y2( 34),  Y2( 35),  Y2( 36)
     +    / 2.18D0, -0.24D0,  0.58D0, -0.18D0, -1.55D0, -0.64D0/
      DATA   Y2( 37),  Y2( 38),  Y2( 39),  Y2( 40),  Y2( 41),  Y2( 42)
     +    /-1.09D0,  0.90D0, -0.66D0, -0.35D0,  0.48D0,  0.50D0/
      DATA   Y2( 43),  Y2( 44),  Y2( 45),  Y2( 46),  Y2( 47),  Y2( 48)
     +    / 0.05D0, -0.68D0,  0.24D0,  0.58D0, -1.26D0, -0.25D0/
      DATA   Y2( 49),  Y2( 50),  Y2( 51),  Y2( 52),  Y2( 53),  Y2( 54)
     +    / 0.25D0,  2.18D0,  2.96D0,  1.56D0, -0.36D0, -0.59D0/
      DATA   Y2( 55),  Y2( 56),  Y2( 57),  Y2( 58),  Y2( 59),  Y2( 60)
     +    /-0.12D0,  3.03D0,  2.11D0,  0.78D0,  0.89D0, -1.45D0/
      DATA   Y2( 61),  Y2( 62),  Y2( 63),  Y2( 64),  Y2( 65),  Y2( 66)
     +    /-0.36D0, -0.37D0, -1.39D0, -4.19D0, -0.73D0, -0.98D0/
      DATA   Y2( 67),  Y2( 68),  Y2( 69),  Y2( 70),  Y2( 71),  Y2( 72)
     +    / 0.36D0,  0.06D0, -1.94D0, -0.08D0,  0.17D0,  1.00D0/
      DATA   Y2( 73),  Y2( 74),  Y2( 75),  Y2( 76),  Y2( 77),  Y2( 78)
     +    /-0.05D0,  0.43D0,  0.15D0,  2.69D0,  0.57D0,  0.29D0/
      DATA   Y2( 79),  Y2( 80),  Y2( 81),  Y2( 82),  Y2( 83),  Y2( 84)
     +    / 1.10D0,  0.48D0, -1.06D0, -2.28D0, -2.03D0, -0.75D0/
      DATA   Y2( 85),  Y2( 86),  Y2( 87),  Y2( 88),  Y2( 89),  Y2( 90)
     +    / 1.00D0,  1.71D0,  0.58D0,  1.97D0,  0.99D0,  1.94D0/
      DATA   Y2( 91),  Y2( 92),  Y2( 93),  Y2( 94),  Y2( 95),  Y2( 96)
     +    / 2.18D0,  3.14D0,  0.60D0,  0.51D0,  1.35D0,  0.56D0/
      DATA   Y2( 97),  Y2( 98),  Y2( 99),  Y2(100)
     +    / 0.11D0,  0.00D0,  2.34D0,  1.88D0/
C
      CALL IPRINT(IPRT)
      ITEST = 1
      LDSTAK = LDS
C
      N = 100
      LAGMAX = 20
      NLAG = 30
      NPRT = 1
      LYFFT = 150
      ICCOV = 30
      JCCOV = 5
      IYM = 150
      M = 4
      IYMFFT = 150
      INLPPC = 30
      JNLPPC = 5
      NYD = 144
      YMISS0 = 1.16D0
      Y1MISS = 0.89D0
      Y2MISS = 0.89D0
C
C     COPY DATA INTO YM FOR CCFS AND CCFMS
C
      CALL DCOPY(N, Y1, 1, YM(1,1), 1)
      CALL DCOPY(N, Y2, 1, YM(1,2), 1)
      CALL DCOPY(N, Y1, 1, YM(1,3), 1)
      CALL DCOPY(N, Y2, 1, YM(1,4), 1)
      CALL SETRV(YMMISS, 4, YMISS0)
C
C     TEST OF CCF
C
      WRITE (IPRT,1060)
      CALL CCF(Y1, Y2, N)
C
C     PRINT RETURNED RESULTS
C
      CALL CCFXP (.FALSE., LAGMAX, M, CCOV, ICCOV, JCCOV, .FALSE.,
     +   NLPPC,  INLPPC, JNLPPC, CMISS)
C
C     TEST OF CCFS
C
      WRITE (IPRT,1080)
      CALL CCFS(YM, N, M, IYM, LAGMAX, CCOV, ICCOV, JCCOV, NPRT,
     +   LDSTAK)
C
C     PRINT RETURNED RESULTS
C
      CALL CCFXP (.TRUE., LAGMAX, M, CCOV, ICCOV, JCCOV, .FALSE.,
     +   NLPPC,  INLPPC, JNLPPC, CMISS)
C
C     TEST OF CCFM WITHOUT MISSING VALUES
C
      WRITE (IPRT,1070)
      WRITE (IPRT, 1050)
      CALL CCFM(Y1, YMISS0, Y2, YMISS0, N)
C
C     PRINT RETURNED RESULTS
C
      CALL CCFXP (.FALSE., LAGMAX, M, CCOV, ICCOV, JCCOV, .TRUE.,
     +   NLPPC,  INLPPC, JNLPPC, CMISS)
C
C     TEST OF CCFMS WITHOUT MISSING VALUES
C
      WRITE (IPRT,1140)
      WRITE (IPRT, 1050)
      CALL CCFMS(YM, YMMISS, N, M, IYM, LAGMAX, CCOV, CMISS,
     +   ICCOV, JCCOV, NLPPC, INLPPC, JNLPPC, NPRT, LDSTAK)
C
C     PRINT RETURNED RESULTS
C
      CALL CCFXP (.TRUE., LAGMAX, M, CCOV, ICCOV, JCCOV, .TRUE.,
     +   NLPPC,  INLPPC, JNLPPC, CMISS)
C
C     COPY DATA INTO YFFT1, YFFT2 AND YMFFT FOR CCFF AND CCFFS
C
      CALL DCOPY(N, Y1, 1, YFFT1, 1)
      CALL DCOPY(N, Y2, 1, YFFT2, 1)
      CALL DCOPY(N, Y1, 1, YMFFT(1,1), 1)
      CALL DCOPY(N, Y2, 1, YMFFT(1,2), 1)
      CALL DCOPY(N, Y1, 1, YMFFT(1,3), 1)
      CALL DCOPY(N, Y2, 1, YMFFT(1,4), 1)
C
C     TEST OF CCFF
C
      WRITE (IPRT,1100)
      CALL CCFF(YFFT1, YFFT2, N, LYFFT, LDSTAK)
C
C     PRINT RETURNED RESULTS
C
      CALL CCFXP (.FALSE., LAGMAX, M, CCOV, ICCOV, JCCOV, .FALSE.,
     +   NLPPC,  INLPPC, JNLPPC, CMISS)
C
C     TEST OF CCFFS
C
      WRITE (IPRT,1150)
      CALL CCFFS(YMFFT, N, M, IYMFFT, LAGMAX, CCOV,
     +   ICCOV, JCCOV, NPRT, LDSTAK)
C
C     PRINT RETURNED RESULTS
C
      CALL CCFXP (.TRUE., LAGMAX, M, CCOV, ICCOV, JCCOV, .FALSE.,
     +   NLPPC,  INLPPC, JNLPPC, CMISS)
C
C     RESET YMMISS
C
      YMMISS(1) = Y1MISS
      YMMISS(2) = Y2MISS
      YMMISS(3) = Y1MISS
      YMMISS(4) = Y2MISS
C
C     TEST OF CCFM WITH MISSING VALUES
C
      WRITE (IPRT,1070)
      WRITE (IPRT, 1040)
      CALL CCFM(Y1, Y1MISS, Y2, Y2MISS, N)
C
C     PRINT RETURNED RESULTS
C
      CALL CCFXP (.FALSE., LAGMAX, M, CCOV, ICCOV, JCCOV, .TRUE.,
     +   NLPPC,  INLPPC, JNLPPC, CMISS)
C
C     TEST OF CCFMS WITH MISSING VALUES
C
      WRITE (IPRT,1140)
      WRITE (IPRT, 1040)
      CALL CCFMS(YM, YMMISS, N, M, IYM, LAGMAX, CCOV, CMISS,
     +   ICCOV, JCCOV, NLPPC, INLPPC, JNLPPC, NPRT, LDSTAK)
C
C     PRINT RETURNED RESULTS
C
      CALL CCFXP (.TRUE., LAGMAX, M, CCOV, ICCOV, JCCOV, .TRUE.,
     +   NLPPC,  INLPPC, JNLPPC, CMISS)
C
C     TEST PRINT CONTROL
C
      NPRT = 0
C
C     TEST OF CCFS
C
      WRITE (IPRT,1080)
      WRITE (IPRT, 1020)
      CALL CCFS(YM, N, M, LAGMAX, IYM, CCOV, ICCOV, JCCOV, NPRT,
     +   LDSTAK)
C
C     PRINT RETURNED RESULTS
C
      CALL CCFXP (.TRUE., LAGMAX, M, CCOV, ICCOV, JCCOV, .FALSE.,
     +   NLPPC,  INLPPC, JNLPPC, CMISS)
C
C     TEST OF CCFMS WITH MISSING VALUES
C
      WRITE (IPRT,1140)
      WRITE (IPRT, 1040)
      WRITE (IPRT, 1020)
      CALL CCFMS(YM, YMMISS, N, M, IYM, LAGMAX, CCOV, CMISS,
     +   ICCOV, JCCOV, NLPPC, INLPPC, JNLPPC, NPRT, LDSTAK)
C
C     PRINT RETURNED RESULTS
C
      CALL CCFXP (.TRUE., LAGMAX, M, CCOV, ICCOV, JCCOV, .TRUE.,
     +   NLPPC,  INLPPC, JNLPPC, CMISS)
C
C     COPY DATA INTO YMFFT FOR CCFFS
C
      CALL DCOPY(N, Y1, 1, YMFFT(1,1), 1)
      CALL DCOPY(N, Y2, 1, YMFFT(1,2), 1)
      CALL DCOPY(N, Y1, 1, YMFFT(1,3), 1)
      CALL DCOPY(N, Y2, 1, YMFFT(1,4), 1)
C
C     TEST OF CCFFS
C
      WRITE (IPRT,1150)
      WRITE (IPRT, 1020)
      CALL CCFFS(YMFFT, N, M, IYMFFT, LAGMAX, CCOV,
     +   ICCOV, JCCOV, NPRT, LDSTAK)
C
C     PRINT RETURNED RESULTS
C
      CALL CCFXP (.TRUE., LAGMAX, M, CCOV, ICCOV, JCCOV, .FALSE.,
     +   NLPPC,  INLPPC, JNLPPC, CMISS)
C
C     TEST LEAD/LAG MESSAGE
C
      NPRT = 1
C
      CALL SETRA(YMFFT, IYMFFT, M, N, 0.0D0)
      YMFFT(5,1) = 1.0D0
      YMFFT(15,2) = 1.0D0
      YMFFT(5,3) = YMFFT(5,1)
      YMFFT(15,4) = YMFFT(15,2)
C
C     TEST OF CCFFS
C
      WRITE (IPRT,1150)
      WRITE (IPRT, 1020)
      CALL CCFFS(YMFFT, N, M, IYMFFT, LAGMAX, CCOV,
     +   ICCOV, JCCOV, NPRT, LDSTAK)
C
C     PRINT RETURNED RESULTS
C
      CALL CCFXP (.TRUE., LAGMAX, M, CCOV, ICCOV, JCCOV, .FALSE.,
     +   NLPPC,  INLPPC, JNLPPC, CMISS)
C
      GO TO (100, 200, 300, 400), ITEST
C
C     TEST MINIMUM PROBLEM SIZE
C
  100 ITEST = ITEST + 1
      N = 3
      LAGMAX = 1
      LYFFT = 150
      ICCOV = 30
      JCCOV = 5
      IYM = 150
      M = 1
      IYMFFT = 150
      INLPPC = 30
      JNLPPC = 5
      NYD = 144
      YMISS0 = 1.16D0
      Y1MISS = 0.89D0
      Y2MISS = 0.89D0
C
C     TEST ERROR HANDLING
C
  200 ITEST = ITEST + 1
      N = 0
      LAGMAX = 1
      LYFFT = 0
      ICCOV = 0
      JCCOV = 0
      IYM = 0
      M = 0
      IYMFFT = 0
      INLPPC = 0
      JNLPPC = 0
      NYD = 0
C
C     TEST ERROR HANDLING
C
  300 ITEST = ITEST + 1
      N = 100
      LAGMAX = 100
      LYFFT = 0
      ICCOV = 0
      JCCOV = 0
      IYM = 0
      M = 0
      IYMFFT = 0
      INLPPC = 0
      JNLPPC = 0
      NYD = 144
      LDSTAK = 0
C
  400 RETURN
C
C     FORMAT STATEMENTS
C
 1020 FORMAT (18H OUTPUT SUPPRESSED)
 1040 FORMAT (20H WITH MISSING VALUES)
 1050 FORMAT (23H WITHOUT MISSING VALUES)
 1060 FORMAT ('1', 11HTEST OF CCF)
 1070 FORMAT ('1', 12HTEST OF CCFM)
 1080 FORMAT ('1', 12HTEST OF CCFS)
 1100 FORMAT ('1', 12HTEST OF CCFF)
 1140 FORMAT ('1', 13HTEST OF CCFMS)
 1150 FORMAT ('1', 13HTEST OF CCFFS)
      END
