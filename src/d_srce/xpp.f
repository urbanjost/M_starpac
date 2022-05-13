*XPP
      SUBROUTINE XPP(LDSTAK)
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
      DOUBLE PRECISION
     +   XLB,XMISS,XUB,YLB,YMISS,YUB
      INTEGER
     +   ILOG,IPRT,ISIZE,ITEST,IYM,M,NOUT,NY,NYM
C
C  LOCAL ARRAYS
      DOUBLE PRECISION
     +   AIR(144),TIME(144),X(144),Y(144),YM(12,12),YMMISS(144)
      INTEGER
     +   ISYM(144)
C
C  EXTERNAL SUBROUTINES
      EXTERNAL IPRINT,MPP,MPPC,MPPL,MPPM,MPPMC,MPPML,PP,PPC,PPL,PPM,
     +   PPMC,PPML,DCOPY,SETRV,SPP,SPPC,SPPL,SPPM,SPPMC,SPPML
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C  EQUIVALENCES
      EQUIVALENCE (Y(1),YM(1,1))
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION AIR(144)
C        THE AIRLINE DATA.
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
C     INTEGER NOUT
C        USED TO INDICATE HOW MANY OF THE POINTS OUTSIDE THE BOUNDS
C        OF THE PLOT ARE TO BE LISTED.
C     INTEGER NY, NYM
C        THE NUMBER OF OBSERVATIONS IN ARRAYS Y AND YM, RESPECTIVELY.
C     DOUBLE PRECISION TIME(144)
C        THE TIME VALUES FOR THE AIRLINE DATA.
C     DOUBLE PRECISION X(144)
C        VECTOR OF OBSERVATIONS FOR X(HORIZONTAL) COORDINATES
C     DOUBLE PRECISION XLB
C        THE LOWER BOUND FOR THE X-AXIS.  (XLB=XUB INDICATES LIMITS ARE
C        TO BE DETERMINED FROM THE RANGE OF THE DATA.)
C     DOUBLE PRECISION XMISS
C        THE MISSING VALUE CODE FOR THE X-AXIS.
C     DOUBLE PRECISION XUB
C        THE UPPER BOUND FOR THE X-AXIS.  (XLB=XUB INDICATES LIMITS ARE
C        TO BE DETERMINED FROM THE RANGE OF THE DATA.)
C     DOUBLE PRECISION Y(144)
C        VECTOR OF OBSERVATIONS FOR THE Y (VERTICAL) COORDINATES
C     DOUBLE PRECISION YLB
C        THE LOWER BOUND FOR THE Y-AXIS.  (YLB=YUB INDICATES LIMITS ARE
C        TO BE DETERMINED FROM THE RANGE OF THE DATA.)
C     DOUBLE PRECISION YM(12,12)
C        MULTIVARIATE OBSERVATIONS FOR THE Y (VERTICAL) COORDINATES.
C     DOUBLE PRECISION YMISS
C        THE MISSING VALUE CODE FOR THE Y-AXIS.
C     DOUBLE PRECISION YMMISS(144)
C        THE MISSING VALUE CODES FOR EACH COLUMN OF YM.
C     DOUBLE PRECISION YUB
C        THE UPPER BOUND FOR THE Y-AXIS.  (YLB=YUB INDICATES LIMITS ARE
C        TO BE DETERMINED FROM THE RANGE OF THE DATA.)
C
C
      DATA     XMISS,    YMISS
     +    /      7.0D0,    180.0D0/
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
      DATA TIME(  1),TIME(  2),TIME(  3),TIME(  4),TIME(  5),TIME(  6)
     +    /   1.0D0,    2.0D0,    3.0D0,    4.0D0,    5.0D0,    6.0D0/
      DATA TIME(  7),TIME(  8),TIME(  9),TIME( 10),TIME( 11),TIME( 12)
     +    /   7.0D0,    8.0D0,    9.0D0,   10.0D0,   11.0D0,   12.0D0/
      DATA TIME( 13),TIME( 14),TIME( 15),TIME( 16),TIME( 17),TIME( 18)
     +    /  13.0D0,   14.0D0,   15.0D0,   16.0D0,   17.0D0,   18.0D0/
      DATA TIME( 19),TIME( 20),TIME( 21),TIME( 22),TIME( 23),TIME( 24)
     +    /  19.0D0,   20.0D0,   21.0D0,   22.0D0,   23.0D0,   24.0D0/
      DATA TIME( 25),TIME( 26),TIME( 27),TIME( 28),TIME( 29),TIME( 30)
     +    /  25.0D0,   26.0D0,   27.0D0,   28.0D0,   29.0D0,   30.0D0/
      DATA TIME( 31),TIME( 32),TIME( 33),TIME( 34),TIME( 35),TIME( 36)
     +    /  31.0D0,   32.0D0,   33.0D0,   34.0D0,   35.0D0,   36.0D0/
      DATA TIME( 37),TIME( 38),TIME( 39),TIME( 40),TIME( 41),TIME( 42)
     +    /  37.0D0,   38.0D0,   39.0D0,   40.0D0,   41.0D0,   42.0D0/
      DATA TIME( 43),TIME( 44),TIME( 45),TIME( 46),TIME( 47),TIME( 48)
     +    /  43.0D0,   44.0D0,   45.0D0,   46.0D0,   47.0D0,   48.0D0/
      DATA TIME( 49),TIME( 50),TIME( 51),TIME( 52),TIME( 53),TIME( 54)
     +    /  49.0D0,   50.0D0,   51.0D0,   52.0D0,   53.0D0,   54.0D0/
      DATA TIME( 55),TIME( 56),TIME( 57),TIME( 58),TIME( 59),TIME( 60)
     +    /  55.0D0,   56.0D0,   57.0D0,   58.0D0,   59.0D0,   60.0D0/
      DATA TIME( 61),TIME( 62),TIME( 63),TIME( 64),TIME( 65),TIME( 66)
     +    /  61.0D0,   62.0D0,   63.0D0,   64.0D0,   65.0D0,   66.0D0/
      DATA TIME( 67),TIME( 68),TIME( 69),TIME( 70),TIME( 71),TIME( 72)
     +    /  67.0D0,   68.0D0,   69.0D0,   70.0D0,   71.0D0,   72.0D0/
      DATA TIME( 73),TIME( 74),TIME( 75),TIME( 76),TIME( 77),TIME( 78)
     +    /  73.0D0,   74.0D0,   75.0D0,   76.0D0,   77.0D0,   78.0D0/
      DATA TIME( 79),TIME( 80),TIME( 81),TIME( 82),TIME( 83),TIME( 84)
     +    /  79.0D0,   80.0D0,   81.0D0,   82.0D0,   83.0D0,   84.0D0/
      DATA TIME( 85),TIME( 86),TIME( 87),TIME( 88),TIME( 89),TIME( 90)
     +    /  85.0D0,   86.0D0,   87.0D0,   88.0D0,   89.0D0,   90.0D0/
      DATA TIME( 91),TIME( 92),TIME( 93),TIME( 94),TIME( 95),TIME( 96)
     +    /  91.0D0,   92.0D0,   93.0D0,   94.0D0,   95.0D0,   96.0D0/
      DATA TIME( 97),TIME( 98),TIME( 99),TIME(100),TIME(101),TIME(102)
     +    /  97.0D0,   98.0D0,   99.0D0,  100.0D0,  101.0D0,  102.0D0/
      DATA TIME(103),TIME(104),TIME(105),TIME(106),TIME(107),TIME(108)
     +    / 103.0D0,  104.0D0,  105.0D0,  106.0D0,  107.0D0,  108.0D0/
      DATA TIME(109),TIME(110),TIME(111),TIME(112),TIME(113),TIME(114)
     +    / 109.0D0,  110.0D0,  111.0D0,  112.0D0,  113.0D0,  114.0D0/
      DATA TIME(115),TIME(116),TIME(117),TIME(118),TIME(119),TIME(120)
     +    / 115.0D0,  116.0D0,  117.0D0,  118.0D0,  119.0D0,  120.0D0/
      DATA TIME(121),TIME(122),TIME(123),TIME(124),TIME(125),TIME(126)
     +    / 121.0D0,  122.0D0,  123.0D0,  124.0D0,  125.0D0,  126.0D0/
      DATA TIME(127),TIME(128),TIME(129),TIME(130),TIME(131),TIME(132)
     +    / 127.0D0,  128.0D0,  129.0D0,  130.0D0,  131.0D0,  132.0D0/
      DATA TIME(133),TIME(134),TIME(135),TIME(136),TIME(137),TIME(138)
     +    / 133.0D0,  134.0D0,  135.0D0,  136.0D0,  137.0D0,  138.0D0/
      DATA TIME(139),TIME(140),TIME(141),TIME(142),TIME(143),TIME(144)
     +    / 139.0D0,  140.0D0,  141.0D0,  142.0D0,  143.0D0,  144.0D0/
C
      DATA  AIR(  1), AIR(  2), AIR(  3), AIR(  4), AIR(  5), AIR(  6)
     +    / 112.0D0,  118.0D0,  132.0D0,  129.0D0,  121.0D0,  135.0D0/
      DATA  AIR(  7), AIR(  8), AIR(  9), AIR( 10), AIR( 11), AIR( 12)
     +    / 148.0D0,  148.0D0,  136.0D0,  119.0D0,  104.0D0,  118.0D0/
      DATA  AIR( 13), AIR( 14), AIR( 15), AIR( 16), AIR( 17), AIR( 18)
     +    / 115.0D0,  126.0D0,  141.0D0,  135.0D0,  125.0D0,  149.0D0/
      DATA  AIR( 19), AIR( 20), AIR( 21), AIR( 22), AIR( 23), AIR( 24)
     +    / 170.0D0,  170.0D0,  158.0D0,  133.0D0,  114.0D0,  140.0D0/
      DATA  AIR( 25), AIR( 26), AIR( 27), AIR( 28), AIR( 29), AIR( 30)
     +    / 145.0D0,  150.0D0,  178.0D0,  163.0D0,  172.0D0,  178.0D0/
      DATA  AIR( 31), AIR( 32), AIR( 33), AIR( 34), AIR( 35), AIR( 36)
     +    / 199.0D0,  199.0D0,  184.0D0,  162.0D0,  146.0D0,  166.0D0/
      DATA  AIR( 37), AIR( 38), AIR( 39), AIR( 40), AIR( 41), AIR( 42)
     +    / 171.0D0,  180.0D0,  193.0D0,  181.0D0,  183.0D0,  218.0D0/
      DATA  AIR( 43), AIR( 44), AIR( 45), AIR( 46), AIR( 47), AIR( 48)
     +    / 230.0D0,  242.0D0,  209.0D0,  191.0D0,  172.0D0,  194.0D0/
      DATA  AIR( 49), AIR( 50), AIR( 51), AIR( 52), AIR( 53), AIR( 54)
     +    / 196.0D0,  196.0D0,  236.0D0,  235.0D0,  229.0D0,  243.0D0/
      DATA  AIR( 55), AIR( 56), AIR( 57), AIR( 58), AIR( 59), AIR( 60)
     +    / 264.0D0,  272.0D0,  237.0D0,  211.0D0,  180.0D0,  201.0D0/
      DATA  AIR( 61), AIR( 62), AIR( 63), AIR( 64), AIR( 65), AIR( 66)
     +    / 204.0D0,  188.0D0,  235.0D0,  227.0D0,  234.0D0,  264.0D0/
      DATA  AIR( 67), AIR( 68), AIR( 69), AIR( 70), AIR( 71), AIR( 72)
     +    / 302.0D0,  293.0D0,  259.0D0,  229.0D0,  203.0D0,  229.0D0/
      DATA  AIR( 73), AIR( 74), AIR( 75), AIR( 76), AIR( 77), AIR( 78)
     +    / 242.0D0,  233.0D0,  267.0D0,  269.0D0,  270.0D0,  315.0D0/
      DATA  AIR( 79), AIR( 80), AIR( 81), AIR( 82), AIR( 83), AIR( 84)
     +    / 364.0D0,  347.0D0,  312.0D0,  274.0D0,  237.0D0,  278.0D0/
      DATA  AIR( 85), AIR( 86), AIR( 87), AIR( 88), AIR( 89), AIR( 90)
     +    / 284.0D0,  277.0D0,  317.0D0,  313.0D0,  318.0D0,  374.0D0/
      DATA  AIR( 91), AIR( 92), AIR( 93), AIR( 94), AIR( 95), AIR( 96)
     +    / 413.0D0,  405.0D0,  355.0D0,  306.0D0,  271.0D0,  306.0D0/
      DATA  AIR( 97), AIR( 98), AIR( 99), AIR(100), AIR(101), AIR(102)
     +    / 315.0D0,  301.0D0,  356.0D0,  348.0D0,  355.0D0,  422.0D0/
      DATA  AIR(103), AIR(104), AIR(105), AIR(106), AIR(107), AIR(108)
     +    / 465.0D0,  467.0D0,  404.0D0,  347.0D0,  305.0D0,  336.0D0/
      DATA  AIR(109), AIR(110), AIR(111), AIR(112), AIR(113), AIR(114)
     +    / 340.0D0,  318.0D0,  362.0D0,  348.0D0,  363.0D0,  435.0D0/
      DATA  AIR(115), AIR(116), AIR(117), AIR(118), AIR(119), AIR(120)
     +    / 491.0D0,  505.0D0,  404.0D0,  359.0D0,  310.0D0,  337.0D0/
      DATA  AIR(121), AIR(122), AIR(123), AIR(124), AIR(125), AIR(126)
     +    / 360.0D0,  342.0D0,  406.0D0,  396.0D0,  420.0D0,  472.0D0/
      DATA  AIR(127), AIR(128), AIR(129), AIR(130), AIR(131), AIR(132)
     +    / 548.0D0,  559.0D0,  463.0D0,  407.0D0,  362.0D0,  405.0D0/
      DATA  AIR(133), AIR(134), AIR(135), AIR(136), AIR(137), AIR(138)
     +    / 417.0D0,  391.0D0,  419.0D0,  461.0D0,  472.0D0,  535.0D0/
      DATA  AIR(139), AIR(140), AIR(141), AIR(142), AIR(143), AIR(144)
     +    / 622.0D0,  606.0D0,  508.0D0,  461.0D0,  390.0D0,  432.0D0/
C
      CALL SETRV(YMMISS, 144, YMISS)
      CALL DCOPY(144, AIR, 1, Y, 1)
      CALL DCOPY(144, TIME, 1, X, 1)
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
      ILOG = -1
      ISIZE = -1
      NOUT = -1
      YLB = 0.0D0
      YUB = 0.0D0
      XLB = 0.0D0
      XUB = 0.0D0
C
   10 CONTINUE
C
C     TEST OF PP
C
      WRITE(IPRT, 1000)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NY
      CALL PP(Y, X, NY)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF PPM
C
      WRITE(IPRT, 1030)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NY
      CALL PPM(Y, YMISS, X, XMISS, NY)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF SPP
C
      WRITE(IPRT, 1120)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NY
      CALL SPP(Y, X, NY, ISYM)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF SPPM
C
      WRITE(IPRT, 1150)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NY
      CALL SPPM(Y, YMISS, X, XMISS, NY, ISYM)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF MPP
C
      WRITE(IPRT, 1060)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NYM
      WRITE (IPRT, 3020) M, IYM
      CALL MPP(YM, X, NYM, M, IYM)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF MPPM
C
      WRITE(IPRT, 1090)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NYM
      WRITE (IPRT, 3020) M, IYM
      CALL MPPM(YM, YMMISS, X, XMISS, NYM, M, IYM)
      WRITE (IPRT, 3000) IERR
C
C
C     LOG OPTION CALLS
C
C
   20 CONTINUE
C
C     TEST OF PPL
C
      WRITE(IPRT, 1010)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NY
      WRITE (IPRT, 3040) ILOG
      CALL PPL(Y, X, NY, ILOG)
C
      WRITE (IPRT, 3000) IERR
C
C     TEST OF PPML
C
      WRITE(IPRT, 1040)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NY
      WRITE (IPRT, 3040) ILOG
      CALL PPML(Y, YMISS, X, XMISS, NY, ILOG)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF SPPL
C
      WRITE(IPRT, 1130)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NY
      WRITE (IPRT, 3040) ILOG
      CALL SPPL(Y, X, NY, ISYM, ILOG)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF SPPML
C
      WRITE(IPRT, 1160)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NY
      WRITE (IPRT, 3040) ILOG
      CALL SPPML(Y, YMISS, X, XMISS, NY, ISYM, ILOG)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF MPPL
C
      WRITE(IPRT, 1070)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NYM
      WRITE (IPRT, 3020) M, IYM
      WRITE (IPRT, 3040) ILOG
      CALL MPPL(YM, X, NYM, M, IYM, ILOG)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF MPPML
C
      WRITE(IPRT, 1100)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NYM
      WRITE (IPRT, 3020) M, IYM
      WRITE (IPRT, 3040) ILOG
      CALL MPPML(YM, YMMISS, X, XMISS, NYM, M, IYM, ILOG)
      WRITE (IPRT, 3000) IERR
C
C
C     TEST OF LONG CALLS
C
C
   30 CONTINUE
C
C     TEST OF PPC
C
      WRITE(IPRT, 1020)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NY
      WRITE (IPRT, 3040) ILOG
      WRITE (IPRT, 3050) ISIZE, NOUT
      WRITE (IPRT, 3070) YLB, YUB, XLB
      WRITE (IPRT, 3080) XUB
      CALL PPC(Y, X, NY, ILOG, ISIZE, NOUT, YLB, YUB, XLB, XUB)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF PPMC
C
      WRITE(IPRT, 1050)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NY
      WRITE (IPRT, 3040) ILOG
      WRITE (IPRT, 3050) ISIZE, NOUT
      WRITE (IPRT, 3070) YLB, YUB, XLB
      WRITE (IPRT, 3080) XUB
      CALL PPMC(Y, YMISS, X, XMISS, NY, ILOG, ISIZE, NOUT, YLB, YUB,
     +   XLB, XUB)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF SPPC
C
      WRITE(IPRT, 1140)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NY
      WRITE (IPRT, 3040) ILOG
      WRITE (IPRT, 3050) ISIZE, NOUT
      WRITE (IPRT, 3070) YLB, YUB, XLB
      WRITE (IPRT, 3080) XUB
      CALL SPPC(Y, X, NY, ISYM, ILOG, ISIZE, NOUT, YLB, YUB, XLB,
     +   XUB)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF SPPMC
C
      WRITE(IPRT, 1170)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NY
      WRITE (IPRT, 3040) ILOG
      WRITE (IPRT, 3050) ISIZE, NOUT
      WRITE (IPRT, 3070) YLB, YUB, XLB
      WRITE (IPRT, 3080) XUB
      CALL SPPMC(Y, YMISS, X, XMISS, NY, ISYM, ILOG, ISIZE, NOUT,
     +   YLB, YUB, XLB, XUB)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF MPPC
C
   40 WRITE(IPRT, 1080)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NYM
      WRITE (IPRT, 3020) M, IYM
      WRITE (IPRT, 3040) ILOG
      WRITE (IPRT, 3050) ISIZE, NOUT
      WRITE (IPRT, 3070) YLB, YUB, XLB
      WRITE (IPRT, 3080) XUB
      CALL MPPC(YM, X, NYM, M, IYM, ILOG, ISIZE, NOUT, YLB, YUB,
     +   XLB, XUB)
      WRITE (IPRT, 3000) IERR
C
C     TEST OF MPPMC
C
   50 WRITE(IPRT, 1110)
      WRITE (IPRT, 3100) ITEST
      WRITE (IPRT, 3010) NYM
      WRITE (IPRT, 3020) M, IYM
      WRITE (IPRT, 3040) ILOG
      WRITE (IPRT, 3050) ISIZE, NOUT
      WRITE (IPRT, 3070) YLB, YUB, XLB
      WRITE (IPRT, 3080) XUB
      CALL MPPMC(YM, YMMISS, X, XMISS, NYM, M, IYM, ILOG, ISIZE, NOUT,
     +   YLB, YUB, XLB, XUB)
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
      NOUT = 0
      YLB = 100.0D0
      YUB = 700.0D0
      XLB = 4.0D0
      XUB = 16.0D0
      GO TO 20
C
  120 ILOG = 2
      ISIZE = 2
      NOUT = 5
      GO TO 20
C
  130 ILOG = 20
      ISIZE = 20
      NOUT = 55
      YUB = 300.0D0
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
      X(1) = 10.0D0
      GO TO 40
C
  160 CALL SETRV(Y, 144, 1.0D0)
      CALL SETRV(X, 144, 1.0D0)
      NYM = 6
      IYM = 12
      M = 6
      NY = 36
      YLB = 0.0D0
      YUB = 0.0D0
      XLB = 0.0D0
      XUB = 0.0D0
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
      XLB = -1.0D0
      YLB = -1.0D0
      GO TO 40
C
  190 IYM = 12
      X(1) = 0.0D0
      Y(1) = 0.0D0
      GO TO 50
C
  200 CALL SETRV(X, 144, XMISS)
      CALL SETRV(Y, 144, YMISS)
      XLB = XUB
      YLB = YUB
      GO TO 50
C
  300 CONTINUE
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT ('1', 10HTEST OF PP)
 1010 FORMAT ('1', 11HTEST OF PPL)
 1020 FORMAT ('1', 11HTEST OF PPC)
 1030 FORMAT ('1', 11HTEST OF PPM)
 1040 FORMAT ('1', 12HTEST OF PPML)
 1050 FORMAT ('1', 12HTEST OF PPMC)
 1060 FORMAT ('1', 11HTEST OF MPP)
 1070 FORMAT ('1', 12HTEST OF MPPL)
 1080 FORMAT ('1', 12HTEST OF MPPC)
 1090 FORMAT ('1', 12HTEST OF MPPM)
 1100 FORMAT ('1', 13HTEST OF MPPML)
 1110 FORMAT ('1', 13HTEST OF MPPMC)
 1120 FORMAT ('1', 11HTEST OF SPP)
 1130 FORMAT ('1', 12HTEST OF SPPL)
 1140 FORMAT ('1', 12HTEST OF SPPC)
 1150 FORMAT ('1', 12HTEST OF SPPM)
 1160 FORMAT ('1', 13HTEST OF SPPML)
 1170 FORMAT ('1', 13HTEST OF SPPMC)
 3000 FORMAT (/8H IERR = , I4)
 3010 FORMAT (' ', 5X, 10H   N     =, I5)
 3020 FORMAT ('+', 20X, 10H / M     =, I5, 10H / IYM   =, I5)
 3040 FORMAT ('+', 65X, 10H / ILOG  =, I5)
 3050 FORMAT (' ',  5X, 10H   ISIZE =, I5, 10H / NOUT  =, I5)
 3070 FORMAT ('+', 50X, 10H / YLB   =, F10.4, 10H / YUB   =, F10.4,
     +   10H / XLB   =, F10.4)
 3080 FORMAT ('+', 110X, 10H / XUB   =, F10.4)
 3100 FORMAT (' ', 13H TEST NUMBER , I5)
      END
