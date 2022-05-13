*XVP
      SUBROUTINE XVP(LDSTAK)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     ROUTINE TO TEST THE PLOTTING SUBROUTINES
C
C     SERIES Y IS THE AIRLINE DATA LISTED ON PAGE 531 OF BOX AND
C     JENKINS.
C
C     WRITTEN BY - JANET R. DONALDSON
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  JANUARY 21, 1982
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
C  LOCAL SCALARS
      REAL
     +   XINC,XLB,YLB,YMISS,YUB
      INTEGER
     +   IBAR,ILOG,IPRT,IRLIN,ISIZE,ITEST,IYM,M,NOUT,NS,NY,NYM
C
C  LOCAL ARRAYS
      REAL
     +   AIR(144),Y(144),YM(12,12),YMMISS(144)
      INTEGER
     +   ISYM(144)
C
C  EXTERNAL SUBROUTINES
      EXTERNAL IPRINT,MVP,MVPC,MVPL,MVPM,MVPMC,MVPML,SCOPY,SETRV,SVP,
     +   SVPC,SVPL,SVPM,SVPMC,SVPML,VP,VPC,VPL,VPM,VPMC,VPML
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C  EQUIVALENCES
      EQUIVALENCE (Y(1),YM(1,1))
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     REAL AIR(144)
C        THE AIRLINE DATA.
C     INTEGER IBAR
C        THE VARIABLE USED TO DETERMINE IF SINGLE POINTS (IBAR .NE. 0)
C        OR BARS (IBAR .EQ. 0) ARE TO BE PLOTTED.
C     INTEGER IERR
C        A COMMON VARIABLE USED AS A FLAG TO INDICATE WHETHER
C        OR NOT THERE ARE ANY ERRORS, IF =0 THEN NO ERRORS.
C     INTEGER ILOG
C        THE TWO DIGIT INTEGER, PQ, USED TO SELECT AXIS SCALE, WHERE
C        P DESIGNATES THE X-AXIS AND Q DESIGNATES THE Y-AXIS.
C        IF P.EQ.0 (Q.EQ.0), THEN THE X-AXIS (Y-AXIS) IS LINEAR.
C        IF P.NE.0 (Q.NE.0), THEN THE X-AXIS (Y-AXIS) IS LOG.
C     INTEGER IPRT
C        OUTPUT LOGICAL UNIT NUMBER
C     INTEGER IRLIN
C        THE INDICATOR VARIABLE USED TO DESIGNATE WHETHER ZERO OR THE
C        SERIES MEAN IS TO BE PLOTTED AS A REFERENCE LINE, OR WHETHER
C        NO REFERENCE LINE IS TO BE PLOTTED.
C        IF IRLIN .LE. -1, NO REFERENCE LINE IS PLOTTED.
C        IF IRLIN .EQ.  0, ZERO IS PLOTTED AS THE REFERENCE LINE.
C        IF IRLIN .GE.  1, THE SERIES MEAN IS PLOTTED.
C     INTEGER ISIZE
C        THE TWO DIGIT INTEGER, PQ, USED TO SELECT AXIS SIZE, WHERE
C        P DESIGNATES THE X-AXIS AND Q DESIGNATES THE Y-AXIS.
C        IF P.EQ.0 (Q.EQ.0), THEN THE X-AXIS (Y-AXIS) IS THE MAXIMUM.
C        IF P.NE.0 (Q.NE.0), THEN THE X-AXIS (Y-AXIS) IS HALF THE MAXIMU
C     INTEGER ISYM(144)
C        VECTOR CONTAINING SYMBOL DESIGNATIONS FOR PLOTTING
C     INTEGER ITEST
C        THE NUMBER OF THE TEST.
C     INTEGER IYM
C        ACTUAL DIMENSION OF YM IN USERS MAIN PROGRAM
C     INTEGER LDSTAK
C        *
C     INTEGER M
C        THE NUMBER OF VECTORS IN YM
C     INTEGER NS
C        THE SAMPLING FREQUENCY,
C        WHERE IF NS .LE. 1, EVERY POINT IS PLOTTED,
C                       = 2, EVERY OTHER POINT IS PLOTTED,
C                       = 3, EVERY THIRD POINT IS PLOTTED, ETC.
C     INTEGER NY, NYM
C        THE NUMBER OF OBSERVATIONS IN ARRAYS Y AND YM, RESPECTIVELY.
C     REAL XINC
C        THE INCREMENT FOR THE X AXIS.
C     REAL XLB
C        THE LOWER BOUND FOR THE X-AXIS.
C     REAL Y(144)
C        VECTOR OF OBSERVATIONS FOR THE Y (VERTICAL) COORDINATES
C     REAL YLB
C        THE LOWER BOUND FOR THE Y-AXIS.  (YLB=YUB INDICATES LIMITS ARE
C        TO BE DETERMINED FROM THE RANGE OF THE DATA.)
C     REAL YM(12,12)
C        MULTIVARIATE OBSERVATIONS FOR THE Y (VERTICAL) COORDINATES.
C     REAL YMISS
C        THE MISSING VALUE CODE FOR THE Y-AXIS.
C     REAL YMMISS(144)
C        THE MISSING VALUE CODES FOR EACH COLUMN OF YM.
C     REAL YUB
C        THE UPPER BOUND FOR THE Y-AXIS.  (YLB=YUB INDICATES LIMITS ARE
C        TO BE DETERMINED FROM THE RANGE OF THE DATA.)
C
C
      DATA YMISS/180.0E0/
C
      DATA ISYM(  1),ISYM(  2),ISYM(  3),ISYM(  4),ISYM(  5),ISYM(  6)
     +    /    -5000,     6000,        7,        8,        9,       10/
      DATA ISYM(  7),ISYM(  8),ISYM(  9),ISYM( 10),ISYM( 11),ISYM( 12)
     +    /       11,       12,       13,       14,       15,       16/
      DATA ISYM( 13),ISYM( 14),ISYM( 15),ISYM( 16),ISYM( 17),ISYM( 18)
     +    /        5,        6,        7,        8,        9,       10/
      DATA ISYM( 19),ISYM( 20),ISYM( 21),ISYM( 22),ISYM( 23),ISYM( 24)
     +    /       11,       12,       13,       14,       15,       16/
      DATA ISYM( 25),ISYM( 26),ISYM( 27),ISYM( 28),ISYM( 29),ISYM( 30)
     +    /        5,        6,        7,        8,        9,       10/
      DATA ISYM( 31),ISYM( 32),ISYM( 33),ISYM( 34),ISYM( 35),ISYM( 36)
     +    /       11,       12,       13,       14,       15,       16/
      DATA ISYM( 37),ISYM( 38),ISYM( 39),ISYM( 40),ISYM( 41),ISYM( 42)
     +    /        5,        6,        7,        8,        9,       10/
      DATA ISYM( 43),ISYM( 44),ISYM( 45),ISYM( 46),ISYM( 47),ISYM( 48)
     +    /       11,       12,       13,       14,       15,       16/
      DATA ISYM( 49),ISYM( 50),ISYM( 51),ISYM( 52),ISYM( 53),ISYM( 54)
     +    /        5,        6,        7,        8,        9,       10/
      DATA ISYM( 55),ISYM( 56),ISYM( 57),ISYM( 58),ISYM( 59),ISYM( 60)
     +    /       11,       12,       13,       14,       15,       16/
      DATA ISYM( 61),ISYM( 62),ISYM( 63),ISYM( 64),ISYM( 65),ISYM( 66)
     +    /        5,        6,        7,        8,        9,       10/
      DATA ISYM( 67),ISYM( 68),ISYM( 69),ISYM( 70),ISYM( 71),ISYM( 72)
     +    /       11,       12,       13,       14,       15,       16/
      DATA ISYM( 73),ISYM( 74),ISYM( 75),ISYM( 76),ISYM( 77),ISYM( 78)
     +    /        5,        6,        7,        8,        9,       10/
      DATA ISYM( 79),ISYM( 80),ISYM( 81),ISYM( 82),ISYM( 83),ISYM( 84)
     +    /       11,       12,       13,       14,       15,       16/
      DATA ISYM( 85),ISYM( 86),ISYM( 87),ISYM( 88),ISYM( 89),ISYM( 90)
     +    /        5,        6,        7,        8,        9,       10/
      DATA ISYM( 91),ISYM( 92),ISYM( 93),ISYM( 94),ISYM( 95),ISYM( 96)
     +    /       11,       12,       13,       14,       15,       16/
      DATA ISYM( 97),ISYM( 98),ISYM( 99),ISYM(100),ISYM(101),ISYM(102)
     +    /        5,        6,        7,        8,        9,       10/
      DATA ISYM(103),ISYM(104),ISYM(105),ISYM(106),ISYM(107),ISYM(108)
     +    /       11,       12,       13,       14,       15,       16/
      DATA ISYM(109),ISYM(110),ISYM(111),ISYM(112),ISYM(113),ISYM(114)
     +    /        5,        6,        7,        8,        9,       10/
      DATA ISYM(115),ISYM(116),ISYM(117),ISYM(118),ISYM(119),ISYM(120)
     +    /       11,       12,       13,       14,       15,       16/
      DATA ISYM(121),ISYM(122),ISYM(123),ISYM(124),ISYM(125),ISYM(126)
     +    /        5,        6,        7,        8,        9,       10/
      DATA ISYM(127),ISYM(128),ISYM(129),ISYM(130),ISYM(131),ISYM(132)
     +    /       11,       12,       13,       14,       15,       16/
      DATA ISYM(133),ISYM(134),ISYM(135),ISYM(136),ISYM(137),ISYM(138)
     +    /        5,        6,        7,        8,        9,       10/
      DATA ISYM(139),ISYM(140),ISYM(141),ISYM(142),ISYM(143),ISYM(144)
     +    /       11,       12,       13,       14,       15,       16/
C
      DATA  AIR(  1), AIR(  2), AIR(  3), AIR(  4), AIR(  5), AIR(  6)
     +    / 112.0E0, 118.0E0, 132.0E0, 129.0E0, 121.0E0, 135.0E0/
      DATA  AIR(  7), AIR(  8), AIR(  9), AIR( 10), AIR( 11), AIR( 12)
     +    / 148.0E0, 148.0E0, 136.0E0, 119.0E0, 104.0E0, 118.0E0/
      DATA  AIR( 13), AIR( 14), AIR( 15), AIR( 16), AIR( 17), AIR( 18)
     +    / 115.0E0, 126.0E0, 141.0E0, 135.0E0, 125.0E0, 149.0E0/
      DATA  AIR( 19), AIR( 20), AIR( 21), AIR( 22), AIR( 23), AIR( 24)
     +    / 170.0E0, 170.0E0, 158.0E0, 133.0E0, 114.0E0, 140.0E0/
      DATA  AIR( 25), AIR( 26), AIR( 27), AIR( 28), AIR( 29), AIR( 30)
     +    / 145.0E0, 150.0E0, 178.0E0, 163.0E0, 172.0E0, 178.0E0/
      DATA  AIR( 31), AIR( 32), AIR( 33), AIR( 34), AIR( 35), AIR( 36)
     +    / 199.0E0, 199.0E0, 184.0E0, 162.0E0, 146.0E0, 166.0E0/
      DATA  AIR( 37), AIR( 38), AIR( 39), AIR( 40), AIR( 41), AIR( 42)
     +    / 171.0E0, 180.0E0, 193.0E0, 181.0E0, 183.0E0, 218.0E0/
      DATA  AIR( 43), AIR( 44), AIR( 45), AIR( 46), AIR( 47), AIR( 48)
     +    / 230.0E0, 242.0E0, 209.0E0, 191.0E0, 172.0E0, 194.0E0/
      DATA  AIR( 49), AIR( 50), AIR( 51), AIR( 52), AIR( 53), AIR( 54)
     +    / 196.0E0, 196.0E0, 236.0E0, 235.0E0, 229.0E0, 243.0E0/
      DATA  AIR( 55), AIR( 56), AIR( 57), AIR( 58), AIR( 59), AIR( 60)
     +    / 264.0E0, 272.0E0, 237.0E0, 211.0E0, 180.0E0, 201.0E0/
      DATA  AIR( 61), AIR( 62), AIR( 63), AIR( 64), AIR( 65), AIR( 66)
     +    / 204.0E0, 188.0E0, 235.0E0, 227.0E0, 234.0E0, 264.0E0/
      DATA  AIR( 67), AIR( 68), AIR( 69), AIR( 70), AIR( 71), AIR( 72)
     +    / 302.0E0, 293.0E0, 259.0E0, 229.0E0, 203.0E0, 229.0E0/
      DATA  AIR( 73), AIR( 74), AIR( 75), AIR( 76), AIR( 77), AIR( 78)
     +    / 242.0E0, 233.0E0, 267.0E0, 269.0E0, 270.0E0, 315.0E0/
      DATA  AIR( 79), AIR( 80), AIR( 81), AIR( 82), AIR( 83), AIR( 84)
     +    / 364.0E0, 347.0E0, 312.0E0, 274.0E0, 237.0E0, 278.0E0/
      DATA  AIR( 85), AIR( 86), AIR( 87), AIR( 88), AIR( 89), AIR( 90)
     +    / 284.0E0, 277.0E0, 317.0E0, 313.0E0, 318.0E0, 374.0E0/
      DATA  AIR( 91), AIR( 92), AIR( 93), AIR( 94), AIR( 95), AIR( 96)
     +    / 413.0E0, 405.0E0, 355.0E0, 306.0E0, 271.0E0, 306.0E0/
      DATA  AIR( 97), AIR( 98), AIR( 99), AIR(100), AIR(101), AIR(102)
     +    / 315.0E0, 301.0E0, 356.0E0, 348.0E0, 355.0E0, 422.0E0/
      DATA  AIR(103), AIR(104), AIR(105), AIR(106), AIR(107), AIR(108)
     +    / 465.0E0, 467.0E0, 404.0E0, 347.0E0, 305.0E0, 336.0E0/
      DATA  AIR(109), AIR(110), AIR(111), AIR(112), AIR(113), AIR(114)
     +    / 340.0E0, 318.0E0, 362.0E0, 348.0E0, 363.0E0, 435.0E0/
      DATA  AIR(115), AIR(116), AIR(117), AIR(118), AIR(119), AIR(120)
     +    / 491.0E0, 505.0E0, 404.0E0, 359.0E0, 310.0E0, 337.0E0/
      DATA  AIR(121), AIR(122), AIR(123), AIR(124), AIR(125), AIR(126)
     +    / 360.0E0, 342.0E0, 406.0E0, 396.0E0, 420.0E0, 472.0E0/
      DATA  AIR(127), AIR(128), AIR(129), AIR(130), AIR(131), AIR(132)
     +    / 548.0E0, 559.0E0, 463.0E0, 407.0E0, 362.0E0, 405.0E0/
      DATA  AIR(133), AIR(134), AIR(135), AIR(136), AIR(137), AIR(138)
     +    / 417.0E0, 391.0E0, 419.0E0, 461.0E0, 472.0E0, 535.0E0/
      DATA  AIR(139), AIR(140), AIR(141), AIR(142), AIR(143), AIR(144)
     +    / 622.0E0, 606.0E0, 508.0E0, 461.0E0, 390.0E0, 432.0E0/
C
      CALL SETRV(YMMISS, 144, YMISS)
      CALL SCOPY(144, AIR, 1, Y, 1)
C
C     DEFINE CONSTANTS
C
      CALL IPRINT(IPRT)
C
C     COMMENCE BODY OF ROUTINE
C
      ITEST = 0
C
C     SHORT CALLS
C
      NY = 144
      NYM = 12
      IYM = 12
      M = 12
      NS = 1
      ILOG = -1
      ISIZE = -1
      ISIZE = -1
      IRLIN = -1
      IBAR = -1
      YLB = 0.0E0
      YUB = 0.0E0
      XLB = 0.0E0
      XINC = 0.0E0
C
   10 CONTINUE
C
C     TEST OF VP
C
      WRITE(IPRT, 2000)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NY
      WRITE (IPRT, 3030) NS
      CALL VP(Y, NY, NS)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF VPM
C
      WRITE(IPRT, 2030)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NY
      WRITE (IPRT, 3030) NS
      CALL VPM (Y, YMISS, NY, NS)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF SVP
C
      WRITE(IPRT, 2120)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NY
      WRITE (IPRT, 3030) NS
      CALL SVP (Y, NY, NS, ISYM)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF SVPM
C
      WRITE(IPRT, 2150)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NY
      WRITE (IPRT, 3030) NS
      CALL SVPM (Y, YMISS, NY, NS, ISYM)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF MVP
C
      WRITE(IPRT, 2060)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NYM
      WRITE (IPRT, 3020) M, IYM
      WRITE (IPRT, 3030) NS
      CALL MVP (YM, NYM, M, IYM, NS)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF MVPM
C
      WRITE(IPRT, 2090)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NYM
      WRITE (IPRT, 3020) M, IYM
      WRITE (IPRT, 3030) NS
      CALL MVPM (YM, YMMISS, NYM, M, IYM, NS)
      WRITE (IPRT, 3000) IERR
C
C     LOG OPTION CALLS
C
   20 CONTINUE
C
C     TEST OF VPL
C
      WRITE(IPRT, 2010)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NY
      WRITE (IPRT, 3030) NS
      WRITE (IPRT, 3040) ILOG
      CALL VPL (Y, NY, NS, ILOG)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF VPML
C
      WRITE(IPRT, 2040)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NY
      WRITE (IPRT, 3030) NS
      WRITE (IPRT, 3040) ILOG
      CALL VPML (Y, YMISS, NY, NS, ILOG)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF SVPL
C
      WRITE(IPRT, 2130)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NY
      WRITE (IPRT, 3030) NS
      WRITE (IPRT, 3040) ILOG
      CALL SVPL (Y, NY, NS, ISYM, ILOG)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF SVPML
C
      WRITE(IPRT, 2160)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NY
      WRITE (IPRT, 3030) NS
      WRITE (IPRT, 3040) ILOG
      CALL SVPML (Y, YMISS, NY, NS, ISYM, ILOG)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF MVPL
C
      WRITE(IPRT, 2070)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NYM
      WRITE (IPRT, 3020) M, IYM
      WRITE (IPRT, 3030) NS
      WRITE (IPRT, 3040) ILOG
      CALL MVPL (YM, NYM, M, IYM, NS, ILOG)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF MVPML
C
      WRITE(IPRT, 2100)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NYM
      WRITE (IPRT, 3020) M, IYM
      WRITE (IPRT, 3030) NS
      WRITE (IPRT, 3040) ILOG
      CALL MVPML(YM, YMMISS, NYM, M, IYM, NS, ILOG)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF LONG CALLS
C
   30 CONTINUE
C
C     TEST OF VPC
C
      WRITE(IPRT, 2020)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NY
      WRITE (IPRT, 3030) NS
      WRITE (IPRT, 3040) ILOG
      WRITE (IPRT, 3060) ISIZE, IRLIN, IBAR
      WRITE (IPRT, 3070) YLB, YUB, XLB
      WRITE (IPRT, 3090) XINC
      CALL VPC (Y, NY, NS, ILOG, ISIZE, IRLIN, IBAR, YLB,
     +   YUB, XLB, XINC)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF VPMC
C
      WRITE(IPRT, 2050)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NY
      WRITE (IPRT, 3030) NS
      WRITE (IPRT, 3040) ILOG
      WRITE (IPRT, 3060) ISIZE, IRLIN, IBAR
      WRITE (IPRT, 3070) YLB, YUB, XLB
      WRITE (IPRT, 3090) XINC
      CALL VPMC (Y, YMISS, NY, NS, ILOG, ISIZE, IRLIN, IBAR, YLB,
     +   YUB, XLB, XINC)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF SVPC
C
      WRITE(IPRT, 2140)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NY
      WRITE (IPRT, 3030) NS
      WRITE (IPRT, 3040) ILOG
      WRITE (IPRT, 3060) ISIZE, IRLIN, IBAR
      WRITE (IPRT, 3070) YLB, YUB, XLB
      WRITE (IPRT, 3090) XINC
      CALL SVPC (Y, NY, NS, ISYM, ILOG, ISIZE, IRLIN, IBAR, YLB,
     +   YUB, XLB, XINC)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF SVPMC
C
      WRITE(IPRT, 2170)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NY
      WRITE (IPRT, 3030) NS
      WRITE (IPRT, 3040) ILOG
      WRITE (IPRT, 3060) ISIZE, IRLIN, IBAR
      WRITE (IPRT, 3070) YLB, YUB, XLB
      WRITE (IPRT, 3090) XINC
      CALL SVPMC(Y, YMISS, NY, NS, ISYM, ILOG, ISIZE, IRLIN, IBAR,
     +   YLB, YUB, XLB, XINC)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF MVPC
C
   40 WRITE(IPRT, 2080)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NYM
      WRITE (IPRT, 3020) M, IYM
      WRITE (IPRT, 3030) NS
      WRITE (IPRT, 3040) ILOG
      WRITE (IPRT, 3060) ISIZE, IRLIN, IBAR
      WRITE (IPRT, 3070) YLB, YUB, XLB
      WRITE (IPRT, 3090) XINC
      CALL MVPC(YM, NYM, M, IYM, NS, ILOG, ISIZE, YLB,
     +   YUB, XLB, XINC)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF MVPMC
C
   50 WRITE(IPRT, 2110)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NYM
      WRITE (IPRT, 3020) M, IYM
      WRITE (IPRT, 3030) NS
      WRITE (IPRT, 3040) ILOG
      WRITE (IPRT, 3060) ISIZE, IRLIN, IBAR
      WRITE (IPRT, 3070) YLB, YUB, XLB
      WRITE (IPRT, 3090) XINC
      CALL MVPMC(YM, YMMISS, NYM, M, IYM, NS, ILOG, ISIZE, YLB,
     +   YUB, XLB, XINC)
      WRITE (IPRT, 3000) IERR
C
      ITEST = ITEST + 1
C
      GO TO (110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 300),
     +  ITEST
C
C     TEST VALID OPTIONS
C
  110 ILOG = 0
      ISIZE = 0
      YLB = 100.0E0
      YUB = 700.0E0
      XLB = 4.0E0
      XINC = 16.0E0
      GO TO 20
C
  120 ILOG = 2
      ISIZE = 2
      NOUT = 5
      XINC = -1.0E0
      GO TO 20
C
  130 ILOG = 20
      ISIZE = 20
      NOUT = 55
      YUB = 300.0E0
      GO TO 30
C
  140 ILOG = 22
      ISIZE = 22
      GO TO 40
C
  150 NY = 1
      NYM = 1
      M = 144
      IYM = 1
      GO TO 40
C
  160 CALL SETRV(Y, 144, 1.0E0)
      NYM = 6
      IYM = 12
      M = 6
      NY = 36
      YLB = 0.0E0
      YUB = 0.0E0
      XLB = 0.0E0
      XINC = 0.0E0
      GO TO 30
C
C     TEST ERROR RESPONSE
C
  170 NY = 0
      NYM = 0
      M = 0
      IYM = -1
      GO TO 10
C
  180 NY = 144
      NYM = 12
      M = 12
      IYM = -1
      XLB = -1.0E0
      YLB = -1.0E0
      GO TO 40
C
  190 IYM = 12
      Y(1) = 0.0E0
      GO TO 50
C
  200 CALL SETRV(Y, 144, YMISS)
      XLB = XINC
      YLB = YUB
      GO TO 50
C
  300 CONTINUE
C
      RETURN
C
C     FORMAT STATEMENTS
C
 2000 FORMAT ('1', 10HTEST OF VP)
 2010 FORMAT ('1', 11HTEST OF VPL)
 2020 FORMAT ('1', 11HTEST OF VPC)
 2030 FORMAT ('1', 11HTEST OF VPM)
 2040 FORMAT ('1', 12HTEST OF VPML)
 2050 FORMAT ('1', 12HTEST OF VPMC)
 2060 FORMAT ('1', 11HTEST OF MVP)
 2070 FORMAT ('1', 12HTEST OF MVPL)
 2080 FORMAT ('1', 12HTEST OF MVPC)
 2090 FORMAT ('1', 12HTEST OF MVPM)
 2100 FORMAT ('1', 13HTEST OF MVPML)
 2110 FORMAT ('1', 13HTEST OF MVPMC)
 2120 FORMAT ('1', 11HTEST OF SVP)
 2130 FORMAT ('1', 12HTEST OF SVPL)
 2140 FORMAT ('1', 12HTEST OF SVPC)
 2150 FORMAT ('1', 12HTEST OF SVPM)
 2160 FORMAT ('1', 13HTEST OF SVPML)
 2170 FORMAT ('1', 13HTEST OF SVPMC)
 3000 FORMAT (/8H IERR = , I4)
 3010 FORMAT (' ', 5X, 10H   N     =, I5)
 3020 FORMAT ('+', 20X, 10H / M     =, I5, 10H / IYM   =, I5)
 3030 FORMAT ('+', 50X, 10H / NS    =, I5)
 3040 FORMAT ('+', 65X, 10H / ILOG  =, I5)
 3060 FORMAT (' ',  5X, '   ISIZE=', I5, ' / IRLIN=', I5,
     +   10H / IBAR  =, I5)
 3070 FORMAT ('+', 50X, 10H / YLB   =, F10.4, 10H / YUB   =, F10.4,
     +   10H / XLB   =, F10.4)
 3090 FORMAT ('+', 110X, 10H / XINC  =, F10.4)
 3100 FORMAT (' ', 13H TEST NUMBER , I5)
      END
