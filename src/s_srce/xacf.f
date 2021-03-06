*XACF
      SUBROUTINE XACF(LDS)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     ROUTINE TO TEST THE TIME SERIES CORRELATION SUBROUTINES
C
C     SERIES Y IS LISTED AS SERIES X1 ON PAGE 362 IN JENKINS AND WATTS.
C
C     SERIES YD IS LISTED AS SERIES G ON PAGE 531 OF BOX AND JENKINS.
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
      REAL
     +   AMISS,YMISS
      INTEGER
     +   I,IAR,IPRT,ITEST,LACOV,LAGMAX,LDSTAK,LYFFT,N,NFAC,NPRT,NYD
C
C  LOCAL ARRAYS
      REAL
     +   ACOV(21),PHI(21),Y(100),YD(150),YFFT(150)
      INTEGER
     +   IOD(2),ND(2),NLPPA(21)
C
C  EXTERNAL SUBROUTINES
      EXTERNAL ACF,ACFD,ACFF,ACFFS,ACFM,ACFMS,ACFS,IPRINT
C
C  COMMON BLOCKS
      COMMON /CSTAK/DSTAK
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     REAL ACOV(21)
C        THE AUTOCOVARIANCE VECTOR.
C     REAL AMISS
C        THE MISSING VALUE CODE FOR THE RETURNED ACVF ESTIMATES
C        (VECTOR ACOV).
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER I
C        AN INDEXING VARIABLE.
C     INTEGER IAR
C        THE ORDER OF THE AUTOREGRESSIVE PROCESS CHOSEN.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED
C     INTEGER IOD(2)
C        THE ORDER OF EACH OF THE DIFFERENCE FACTORS.
C     INTEGER IPRT
C        THE UNIT NUMBER FOR PRINTED OUTPUT.
C     INTEGER ITEST
C        THE NUMBER OF THE TEST BEING RUN
C     INTEGER LACOV
C        THE LENGTH OF THE ACVF RELATED VECTORS.
C     INTEGER LAGMAX
C        THE MAXIMUM LAG VALUE REQUESTED.
C     INTEGER LDSTAK
C        THE LENGTH OF THE ARRAY DSTAK.
C     INTEGER LYFFT
C        THE LENGTH OF THE ARRAYS USED WHEN THE COMPUTATIONS ARE
C        PERFORMED BY THE FFT.
C     INTEGER N
C        THE INTEGER NUMBER OF OBSERVATIONS IN EACH SERIES
C     INTEGER ND(2)
C        THE ARRAY CONTAINING THE NUMBER OF TIMES THE DIFFERENCE
C        FACTORS ARE TO BE APPLIED.
C     INTEGER NFAC
C        THE NUMBER OF DIFFERENCE FACTORS.
C     INTEGER NLPPA(21)
C        THE ARRAY CONTAINING THE NUMBER OF LAGGED PRODUCT PAIRS
C        USED TO COMPUTE EACH ACVF ESTIMATE.
C     INTEGER NPRT
C        THE INDICATOR VARIABLE USED TO SPECIFY WHETHER OR NOT
C        PRINTED OUTPUT IS TO BE GIVEN, WHERE IF THE VALUE OF
C        NPRT IS ZERO, NO OUTPUT IS MADE.
C     INTEGER NYD
C        THE NUMBER OF OBSERVATIONS IN THE SERIES TO BE DIFFERENCED.
C     REAL PHI(21)
C        THE ARRAY OF AUTOREGRESSIVE COEFFICIENTS FOR THE SELECTED
C        ORDER.
C     REAL Y(100), YD(150)
C        THE VECTOR CONTAINING THE OBSERVED TIME SERIES
C     REAL YFFT(150)
C        THE VECTORS USED FOR STORING THE SERIES FOR THE ROUTINES
C        USING THE FFT.
C     REAL YMISS
C        THE MISSING VALUE CODES FOR SERIES Y AND YM.
C
      DATA    Y(  1),   Y(  2),   Y(  3),   Y(  4),   Y(  5),   Y(  6)
     +    / -2.07E0, -1.15E0,  0.69E0, -0.46E0, -1.49E0, -0.70E0/
      DATA    Y(  7),   Y(  8),   Y(  9),   Y( 10),   Y( 11),   Y( 12)
     +    / -1.07E0, -0.69E0, -0.68E0,  1.27E0, -1.05E0, -0.05E0/
      DATA    Y( 13),   Y( 14),   Y( 15),   Y( 16),   Y( 17),   Y( 18)
     +    / -0.84E0, -0.62E0, -0.49E0, -1.29E0, -0.49E0, -1.06E0/
      DATA    Y( 19),   Y( 20),   Y( 21),   Y( 22),   Y( 23),   Y( 24)
     +    / -0.38E0, -0.52E0, -0.13E0,  1.30E0, -1.51E0, -0.43E0/
      DATA    Y( 25),   Y( 26),   Y( 27),   Y( 28),   Y( 29),   Y( 30)
     +    / -1.33E0, -0.78E0,  0.31E0, -0.95E0, -0.90E0, -0.30E0/
      DATA    Y( 31),   Y( 32),   Y( 33),   Y( 34),   Y( 35),   Y( 36)
     +    / -1.02E0, -0.53E0,  0.15E0,  1.40E0,  1.22E0,  0.59E0/
      DATA    Y( 37),   Y( 38),   Y( 39),   Y( 40),   Y( 41),   Y( 42)
     +    /  0.70E0,  1.70E0,  2.78E0,  1.98E0,  1.39E0,  1.85E0/
      DATA    Y( 43),   Y( 44),   Y( 45),   Y( 46),   Y( 47),   Y( 48)
     +    /  2.60E0,  0.51E0,  2.77E0,  1.16E0,  1.07E0, -0.48E0/
      DATA    Y( 49),   Y( 50),   Y( 51),   Y( 52),   Y( 53),   Y( 54)
     +    / -0.52E0,  0.37E0,  0.00E0, -1.99E0, -1.75E0,  0.70E0/
      DATA    Y( 55),   Y( 56),   Y( 57),   Y( 58),   Y( 59),   Y( 60)
     +    /  0.73E0,  1.16E0,  0.06E0, -0.02E0,  1.10E0, -0.35E0/
      DATA    Y( 61),   Y( 62),   Y( 63),   Y( 64),   Y( 65),   Y( 66)
     +    / -1.67E0, -1.57E0,  1.16E0,  1.84E0,  3.35E0,  0.40E0/
      DATA    Y( 67),   Y( 68),   Y( 69),   Y( 70),   Y( 71),   Y( 72)
     +    /  0.45E0,  1.30E0,  0.93E0,  1.17E0, -1.74E0, -1.28E0/
      DATA    Y( 73),   Y( 74),   Y( 75),   Y( 76),   Y( 77),   Y( 78)
     +    / -0.07E0,  1.50E0,  0.53E0,  0.20E0, -0.42E0,  1.18E0/
      DATA    Y( 79),   Y( 80),   Y( 81),   Y( 82),   Y( 83),   Y( 84)
     +    /  0.82E0,  1.50E0,  2.92E0,  1.18E0,  1.23E0,  3.16E0/
      DATA    Y( 85),   Y( 86),   Y( 87),   Y( 88),   Y( 89),   Y( 90)
     +    /  0.79E0,  0.68E0,  1.14E0,  1.02E0,  1.02E0, -0.71E0/
      DATA    Y( 91),   Y( 92),   Y( 93),   Y( 94),   Y( 95),   Y( 96)
     +    / -0.17E0, -1.50E0, -0.26E0, -0.38E0,  0.93E0, -0.33E0/
      DATA    Y( 97),   Y( 98),   Y( 99),   Y(100)
     +    / -1.12E0, -2.95E0, -2.09E0, -1.11E0                    /
C
      DATA   YD(  1),  YD(  2),  YD(  3),  YD(  4),  YD(  5),  YD(  6)
     +    /  112.0E0, 118.0E0, 132.0E0, 129.0E0, 121.0E0, 135.0E0/
      DATA   YD(  7),  YD(  8),  YD(  9),  YD( 10),  YD( 11),  YD( 12)
     +    /  148.0E0, 148.0E0, 136.0E0, 119.0E0, 104.0E0, 118.0E0/
      DATA   YD( 13),  YD( 14),  YD( 15),  YD( 16),  YD( 17),  YD( 18)
     +    /  115.0E0, 126.0E0, 141.0E0, 135.0E0, 125.0E0, 149.0E0/
      DATA   YD( 19),  YD( 20),  YD( 21),  YD( 22),  YD( 23),  YD( 24)
     +    /  170.0E0, 170.0E0, 158.0E0, 133.0E0, 114.0E0, 140.0E0/
      DATA   YD( 25),  YD( 26),  YD( 27),  YD( 28),  YD( 29),  YD( 30)
     +    /  145.0E0, 150.0E0, 178.0E0, 163.0E0, 172.0E0, 178.0E0/
      DATA   YD( 31),  YD( 32),  YD( 33),  YD( 34),  YD( 35),  YD( 36)
     +    /  199.0E0, 199.0E0, 184.0E0, 162.0E0, 146.0E0, 166.0E0/
      DATA   YD( 37),  YD( 38),  YD( 39),  YD( 40),  YD( 41),  YD( 42)
     +    /  171.0E0, 180.0E0, 193.0E0, 181.0E0, 183.0E0, 218.0E0/
      DATA   YD( 43),  YD( 44),  YD( 45),  YD( 46),  YD( 47),  YD( 48)
     +    /  230.0E0, 242.0E0, 209.0E0, 191.0E0, 172.0E0, 194.0E0/
      DATA   YD( 49),  YD( 50),  YD( 51),  YD( 52),  YD( 53),  YD( 54)
     +    /  196.0E0, 196.0E0, 236.0E0, 235.0E0, 229.0E0, 243.0E0/
      DATA   YD( 55),  YD( 56),  YD( 57),  YD( 58),  YD( 59),  YD( 60)
     +    /  264.0E0, 272.0E0, 237.0E0, 211.0E0, 180.0E0, 201.0E0/
      DATA   YD( 61),  YD( 62),  YD( 63),  YD( 64),  YD( 65),  YD( 66)
     +    /  204.0E0, 188.0E0, 235.0E0, 227.0E0, 234.0E0, 264.0E0/
      DATA   YD( 67),  YD( 68),  YD( 69),  YD( 70),  YD( 71),  YD( 72)
     +    /  302.0E0, 293.0E0, 259.0E0, 229.0E0, 203.0E0, 229.0E0/
      DATA   YD( 73),  YD( 74),  YD( 75),  YD( 76),  YD( 77),  YD( 78)
     +    /  242.0E0, 233.0E0, 267.0E0, 269.0E0, 270.0E0, 315.0E0/
      DATA   YD( 79),  YD( 80),  YD( 81),  YD( 82),  YD( 83),  YD( 84)
     +    /  364.0E0, 347.0E0, 312.0E0, 274.0E0, 237.0E0, 278.0E0/
      DATA   YD( 85),  YD( 86),  YD( 87),  YD( 88),  YD( 89),  YD( 90)
     +    /  284.0E0, 277.0E0, 317.0E0, 313.0E0, 318.0E0, 374.0E0/
      DATA   YD( 91),  YD( 92),  YD( 93),  YD( 94),  YD( 95),  YD( 96)
     +    /  413.0E0, 405.0E0, 355.0E0, 306.0E0, 271.0E0, 306.0E0/
      DATA   YD( 97),  YD( 98),  YD( 99),  YD(100),  YD(101),  YD(102)
     +    /  315.0E0, 301.0E0, 356.0E0, 348.0E0, 355.0E0, 422.0E0/
      DATA   YD(103),  YD(104),  YD(105),  YD(106),  YD(107),  YD(108)
     +    /  465.0E0, 467.0E0, 404.0E0, 347.0E0, 305.0E0, 336.0E0/
      DATA   YD(109),  YD(110),  YD(111),  YD(112),  YD(113),  YD(114)
     +    /  340.0E0, 318.0E0, 362.0E0, 348.0E0, 363.0E0, 435.0E0/
      DATA   YD(115),  YD(116),  YD(117),  YD(118),  YD(119),  YD(120)
     +    /  491.0E0, 505.0E0, 404.0E0, 359.0E0, 310.0E0, 337.0E0/
      DATA   YD(121),  YD(122),  YD(123),  YD(124),  YD(125),  YD(126)
     +    /  360.0E0, 342.0E0, 406.0E0, 396.0E0, 420.0E0, 472.0E0/
      DATA   YD(127),  YD(128),  YD(129),  YD(130),  YD(131),  YD(132)
     +    /  548.0E0, 559.0E0, 463.0E0, 407.0E0, 362.0E0, 405.0E0/
      DATA   YD(133),  YD(134),  YD(135),  YD(136),  YD(137),  YD(138)
     +    /  417.0E0, 391.0E0, 419.0E0, 461.0E0, 472.0E0, 535.0E0/
      DATA   YD(139),  YD(140),  YD(141),  YD(142),  YD(143),  YD(144)
     +    /  622.0E0, 606.0E0, 508.0E0, 461.0E0, 390.0E0, 432.0E0/
C
      CALL IPRINT(IPRT)
      ITEST = 1
      LDSTAK = LDS
C
      N = 100
      LAGMAX = 20
      NPRT = 1
      LYFFT = 150
      LACOV = 21
      NYD = 144
      NFAC = 2
      ND(1) = 1
      ND(2) = 1
      IOD(1) = 12
      IOD(2) = 1
      YMISS = 1.16E0
C
C     TEST OF ACF
C
    5 WRITE (IPRT,1000)
      CALL ACF(Y, N)
      WRITE (IPRT,1010) IERR
C
C     TEST OF ACFS
C
      WRITE (IPRT,1020)
      CALL ACFS(Y, N, LAGMAX, LACOV, ACOV, IAR, PHI, NPRT, LDSTAK)
      WRITE (IPRT,1010) IERR
C
C     PRINT STORAGE FROM ACFS
C
      IF (IERR.EQ.0) THEN
        WRITE (IPRT,1030) (ACOV(I),I=1,LAGMAX+1)
        WRITE (IPRT,1030) (PHI(I),I=1,IAR)
      END IF
C
C     TEST OF ACFD
C
      WRITE (IPRT,1040)
      CALL ACFD(YD, NYD, LAGMAX, NFAC, ND, IOD, LDSTAK)
      WRITE (IPRT,1010) IERR
C
C     TEST OF ACFM
C
      WRITE (IPRT,1050)
      CALL ACFM(Y, YMISS, N)
      WRITE (IPRT,1010) IERR
C
C     TEST OF ACFMS
C
      WRITE (IPRT,1120)
      CALL ACFMS(Y, YMISS, N, LAGMAX, LACOV, ACOV, AMISS, NLPPA, NPRT,
     +   LDSTAK)
      WRITE (IPRT,1010) IERR
C
C     PRINT STORAGE FROM ACFMS
C
      IF (IERR.EQ.0) THEN
        WRITE (IPRT,1030) (ACOV(I),I=1,LAGMAX+1)
        WRITE (IPRT,1140) (NLPPA(I),I=1,LAGMAX+1)
      END IF
C
C     COPY DATA INTO YFFT FOR ACFF
C
      DO 10 I=1,N
         YFFT(I) = Y(I)
   10 CONTINUE
C
C     TEST OF ACFF
C
      WRITE (IPRT,1090)
      CALL ACFF(YFFT, N, LYFFT, LDSTAK)
      WRITE (IPRT,1010) IERR
C
C     COPY DATA INTO YFFT FOR ACFFS
C
      DO 20 I=1,N
         YFFT(I) = Y(I)
   20 CONTINUE
C
C     TEST OF ACFFS
C
      WRITE (IPRT,1130)
      CALL ACFFS(YFFT, N, LYFFT, LDSTAK, LAGMAX, LACOV, ACOV, IAR, PHI,
     +   NPRT)
      WRITE (IPRT,1010) IERR
C
C     PRINT STORAGE FROM ACFFS
C
      IF (IERR.EQ.0) THEN
        WRITE (IPRT,1030) (ACOV(I),I=1,LAGMAX+1)
        WRITE (IPRT,1030) (PHI(I),I=1,IAR)
      END IF
C
      GO TO (100, 200, 300, 400), ITEST
C
C     TEST MINIMUM PROBLEM SIZE
C
  100 ITEST = ITEST + 1
      N = 13
      LAGMAX = 1
      NFAC = 1
      ND(1) = 1
      IOD(1) = 1
      GO TO 5
C
C     CHECK ERROR HANDLING
C
  200 ITEST = ITEST + 1
      N = 0
      LAGMAX = 20
      LYFFT = 0
      LACOV = 0
      NYD = 0
      NFAC = 1
      ND(1) = 0
      IOD(1) = 0
      GO TO 5
C
C     CHECK ERROR HANDLING
C
  300 ITEST = ITEST + 1
      N = 100
      LAGMAX = 0
      LYFFT = 0
      LACOV = 0
      NYD = 144
      NFAC = 0
      LDSTAK = 0
      GO TO 5
C
  400 RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT ('1TEST OF ACF')
 1010 FORMAT (8H IERR IS, I5)
 1020 FORMAT ('1', 12HTEST OF ACFS)
 1030 FORMAT (9F10.5)
 1040 FORMAT ('1', 12HTEST OF ACFD)
 1050 FORMAT ('1', 12HTEST OF ACFM)
 1090 FORMAT ('1', 12HTEST OF ACFF)
 1120 FORMAT ('1', 13HTEST OF ACFMS)
 1130 FORMAT ('1', 13HTEST OF ACFFS)
 1140 FORMAT (9I10)
      END
