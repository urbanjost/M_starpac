*XXCH13
      SUBROUTINE XXCH13(LDSTAK)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS SUBPROGRAM PROVIDES A SIMPLE TEST OF
C     THE ARIMA MODELING AND FORECASTING FAMILY OF ROUTINES.
C
C     DATA IS THE AIRLINE DATA LISTED ON PAGE 531 OF BOX AND JENKINS.
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  AUGUST 3, 1987
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
      INTEGER
     +   I,IPRT,N,NFAC,NPAR
C
C  LOCAL ARRAYS
      DOUBLE PRECISION
     +   AIR(200),PAR(10),RES(200),Y(200)
      INTEGER
     +   MSPEC(4,10)
C
C  EXTERNAL SUBROUTINES
      EXTERNAL AIME,AIMF,IPRINT
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
C     DOUBLE PRECISION AIR(200)
C        THE AIRLINE DATA.
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER I
C        AN INDEX VARIABLE.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .GE. 1, ERRORS WERE DETECTED.
C     INTEGER IPRT
C        THE UNIT NUMBER FOR PRINTED OUTPUT.
C     INTEGER LDSTAK
C        THE LENGTH OF DSTAK IN COMMON /CSTAK/.
C     INTEGER MSPEC(4,10)
C        THE ARRAY CONTAINING THE VALUES OF P, D, Q, AND S FOR EACH
C        FACTOR.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS.
C     INTEGER NFAC
C        THE NUMBER OF FACTORS IN THE MODEL
C     INTEGER NPAR
C        THE NUMBER OF UNKNOWN PARAMETERS IN THE MODEL.
C     DOUBLE PRECISION PAR(10)
C        THE ARRAY IN WHICH THE CURRENT ESTIMATES OF THE UNKNOWN
C        PARAMETERS ARE STORED.
C     DOUBLE PRECISION RES(200)
C        THE RESIDUALS FROM THE FIT.
C     DOUBLE PRECISION Y(200)
C        THE ARRAY OF THE DEPENDENT VARIABLE.
C
C     DEFINE CONSTANTS
C
      DATA MSPEC(1,1), MSPEC(2,1), MSPEC(3,1), MSPEC(4,1)
     +   /          0,          1,          1,          1/
      DATA MSPEC(1,2), MSPEC(2,2), MSPEC(3,2), MSPEC(4,2)
     +   /          0,          1,          1,         12/
C
      DATA  AIR(  1), AIR(  2), AIR(  3), AIR(  4), AIR(  5), AIR(  6)
     +    / 112.0D0, 118.0D0, 132.0D0, 129.0D0, 121.0D0, 135.0D0/
      DATA  AIR(  7), AIR(  8), AIR(  9), AIR( 10), AIR( 11), AIR( 12)
     +    / 148.0D0, 148.0D0, 136.0D0, 119.0D0, 104.0D0, 118.0D0/
      DATA  AIR( 13), AIR( 14), AIR( 15), AIR( 16), AIR( 17), AIR( 18)
     +    / 115.0D0, 126.0D0, 141.0D0, 135.0D0, 125.0D0, 149.0D0/
      DATA  AIR( 19), AIR( 20), AIR( 21), AIR( 22), AIR( 23), AIR( 24)
     +    / 170.0D0, 170.0D0, 158.0D0, 133.0D0, 114.0D0, 140.0D0/
      DATA  AIR( 25), AIR( 26), AIR( 27), AIR( 28), AIR( 29), AIR( 30)
     +    / 145.0D0, 150.0D0, 178.0D0, 163.0D0, 172.0D0, 178.0D0/
      DATA  AIR( 31), AIR( 32), AIR( 33), AIR( 34), AIR( 35), AIR( 36)
     +    / 199.0D0, 199.0D0, 184.0D0, 162.0D0, 146.0D0, 166.0D0/
      DATA  AIR( 37), AIR( 38), AIR( 39), AIR( 40), AIR( 41), AIR( 42)
     +    / 171.0D0, 180.0D0, 193.0D0, 181.0D0, 183.0D0, 218.0D0/
      DATA  AIR( 43), AIR( 44), AIR( 45), AIR( 46), AIR( 47), AIR( 48)
     +    / 230.0D0, 242.0D0, 209.0D0, 191.0D0, 172.0D0, 194.0D0/
      DATA  AIR( 49), AIR( 50), AIR( 51), AIR( 52), AIR( 53), AIR( 54)
     +    / 196.0D0, 196.0D0, 236.0D0, 235.0D0, 229.0D0, 243.0D0/
      DATA  AIR( 55), AIR( 56), AIR( 57), AIR( 58), AIR( 59), AIR( 60)
     +    / 264.0D0, 272.0D0, 237.0D0, 211.0D0, 180.0D0, 201.0D0/
      DATA  AIR( 61), AIR( 62), AIR( 63), AIR( 64), AIR( 65), AIR( 66)
     +    / 204.0D0, 188.0D0, 235.0D0, 227.0D0, 234.0D0, 264.0D0/
      DATA  AIR( 67), AIR( 68), AIR( 69), AIR( 70), AIR( 71), AIR( 72)
     +    / 302.0D0, 293.0D0, 259.0D0, 229.0D0, 203.0D0, 229.0D0/
      DATA  AIR( 73), AIR( 74), AIR( 75), AIR( 76), AIR( 77), AIR( 78)
     +    / 242.0D0, 233.0D0, 267.0D0, 269.0D0, 270.0D0, 315.0D0/
      DATA  AIR( 79), AIR( 80), AIR( 81), AIR( 82), AIR( 83), AIR( 84)
     +    / 364.0D0, 347.0D0, 312.0D0, 274.0D0, 237.0D0, 278.0D0/
      DATA  AIR( 85), AIR( 86), AIR( 87), AIR( 88), AIR( 89), AIR( 90)
     +    / 284.0D0, 277.0D0, 317.0D0, 313.0D0, 318.0D0, 374.0D0/
      DATA  AIR( 91), AIR( 92), AIR( 93), AIR( 94), AIR( 95), AIR( 96)
     +    / 413.0D0, 405.0D0, 355.0D0, 306.0D0, 271.0D0, 306.0D0/
      DATA  AIR( 97), AIR( 98), AIR( 99), AIR(100), AIR(101), AIR(102)
     +    / 315.0D0, 301.0D0, 356.0D0, 348.0D0, 355.0D0, 422.0D0/
      DATA  AIR(103), AIR(104), AIR(105), AIR(106), AIR(107), AIR(108)
     +    / 465.0D0, 467.0D0, 404.0D0, 347.0D0, 305.0D0, 336.0D0/
      DATA  AIR(109), AIR(110), AIR(111), AIR(112), AIR(113), AIR(114)
     +    / 340.0D0, 318.0D0, 362.0D0, 348.0D0, 363.0D0, 435.0D0/
      DATA  AIR(115), AIR(116), AIR(117), AIR(118), AIR(119), AIR(120)
     +    / 491.0D0, 505.0D0, 404.0D0, 359.0D0, 310.0D0, 337.0D0/
      DATA  AIR(121), AIR(122), AIR(123), AIR(124), AIR(125), AIR(126)
     +    / 360.0D0, 342.0D0, 406.0D0, 396.0D0, 420.0D0, 472.0D0/
      DATA  AIR(127), AIR(128), AIR(129), AIR(130), AIR(131), AIR(132)
     +    / 548.0D0, 559.0D0, 463.0D0, 407.0D0, 362.0D0, 405.0D0/
      DATA  AIR(133), AIR(134), AIR(135), AIR(136), AIR(137), AIR(138)
     +    / 417.0D0, 391.0D0, 419.0D0, 461.0D0, 472.0D0, 535.0D0/
      DATA  AIR(139), AIR(140), AIR(141), AIR(142), AIR(143), AIR(144)
     +    / 622.0D0, 606.0D0, 508.0D0, 461.0D0, 390.0D0, 432.0D0/
C
      CALL IPRINT(IPRT)
C
      NFAC = 2
      N = 144
C
      NPAR = 3
      PAR(1) = 0.000
      PAR(2) = 0.395
      PAR(3) = 0.615
C
      DO 10 I = 1, N
        Y(I) = LOG(AIR(I))
   10 CONTINUE
C
C     RUN SIMPLE TEST OF AIME
C
      WRITE (IPRT,1000)
      WRITE (IPRT,1100)
      CALL AIME (Y, N, MSPEC, NFAC, PAR, NPAR, RES, LDSTAK)
      WRITE (IPRT,2000) IERR
C
C     RUN SIMPLE TEST OF AIMF
C
      WRITE (IPRT,1200)
      CALL AIMF (Y, N, MSPEC, NFAC, PAR, NPAR, LDSTAK)
      WRITE (IPRT,2000) IERR
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT ('1*CH13')
 1100 FORMAT (' SIMPLE TEST OF AIME')
 1200 FORMAT ('1SIMPLE TEST OF AIMF')
 2000 FORMAT (/' THE VALUE OF IERR IS ', I4)
      END
