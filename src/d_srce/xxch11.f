*XXCH11
      SUBROUTINE XXCH11(LDSTAK)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS SUBPROGRAM PROVIDES A SIMPLE TEST OF
C     THE COMPLEX DEMODULATION FAMILY OF ROUTINES.
C
C     DATA IS THE WOLF SUMSPOT NUMBERS FOR THE YEARS 1700 TO 1960 AS
C     TABULATED BY WALDMEIER [1961].
C
C     WRITTEN BY - JANET R. DONALDSON
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
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
      DOUBLE PRECISION
     +   FC,FD
      INTEGER
     +   IPRT,K,N
C
C  LOCAL ARRAYS
      DOUBLE PRECISION
     +   Y(300)
C
C  EXTERNAL SUBROUTINES
      EXTERNAL DEMOD,IPRINT
C
C  COMMON BLOCKS
      COMMON /CSTAK/DSTAK
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     DOUBLE PRECISION FC
C        THE CUTOFF FREQUENCY USED FOR THE LOW PASS FILTER.
C     DOUBLE PRECISION FD
C        THE DEMODULATION FREQUENCY.
C     INTEGER IERR
C        A COMMON VARIABLE USED AS A FLAG TO INDICATE WHETHER
C        OR NOT THERE ARE ANY ERRORS, IF =0 THEN NO ERRORS.
C     INTEGER IPRT
C        LOGICAL OUTPUT UNIT.
C     INTEGER K
C        THE NUMBER OF TERMS IN THE SYMETRIC LINEAR FILTER.
C     INTEGER LDSTAK
C        THE LENGTH OF DSTAK IN COMMON /CSTAK/.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS.
C     DOUBLE PRECISION Y(300)
C        THE LOG OF THE AIRLINE DATA.
C
C
      DATA   Y(  1),  Y(  2),  Y(  3),  Y(  4),  Y(  5),  Y(  6)
     +    /     5.0D0, 11.0D0, 16.0D0, 23.0D0, 36.0D0, 58.0D0/
      DATA   Y(  7),  Y(  8),  Y(  9),  Y( 10),  Y( 11),  Y( 12)
     +    /    29.0D0, 20.0D0, 10.0D0,  8.0D0,  3.0D0,  0.0D0/
      DATA   Y( 13),  Y( 14),  Y( 15),  Y( 16),  Y( 17),  Y( 18)
     +    /     0.0D0, 2.0D0, 11.0D0, 27.0D0, 47.0D0, 63.0D0/
      DATA   Y( 19),  Y( 20),  Y( 21),  Y( 22),  Y( 23),  Y( 24)
     +    /    60.0D0, 39.0D0, 28.0D0, 26.0D0, 22.0D0, 11.0D0/
      DATA   Y( 25),  Y( 26),  Y( 27),  Y( 28),  Y( 29),  Y( 30)
     +    /    21.0D0, 40.0D0, 78.0D0,122.0D0,103.0D0, 73.0D0/
      DATA   Y( 31),  Y( 32),  Y( 33),  Y( 34),  Y( 35),  Y( 36)
     +    /    47.0D0, 35.0D0, 11.0D0,  5.0D0, 16.0D0, 34.0D0/
      DATA   Y( 37),  Y( 38),  Y( 39),  Y( 40),  Y( 41),  Y( 42)
     +    /    70.0D0, 81.0D0,111.0D0,101.0D0, 73.0D0, 40.0D0/
      DATA   Y( 43),  Y( 44),  Y( 45),  Y( 46),  Y( 47),  Y( 48)
     +    /    20.0D0, 16.0D0,  5.0D0, 11.0D0, 22.0D0, 40.0D0/
      DATA   Y( 49),  Y( 50),  Y( 51),  Y( 52),  Y( 53),  Y( 54)
     +    /    60.0D0, 80.9D0, 83.4D0, 47.7D0, 47.8D0, 30.7D0/
      DATA   Y( 55),  Y( 56),  Y( 57),  Y( 58),  Y( 59),  Y( 60)
     +    /    12.2D0,  9.6D0, 10.2D0, 32.4D0, 47.6D0, 54.0D0/
      DATA   Y( 61),  Y( 62),  Y( 63),  Y( 64),  Y( 65),  Y( 66)
     +    /    62.9D0, 85.9D0, 61.2D0, 45.1D0, 36.4D0, 20.9D0/
      DATA   Y( 67),  Y( 68),  Y( 69),  Y( 70),  Y( 71),  Y( 72)
     +    /    11.4D0, 37.8D0, 69.8D0,106.1D0,100.8D0, 81.6D0/
      DATA   Y( 73),  Y( 74),  Y( 75),  Y( 76),  Y( 77),  Y( 78)
     +    /    66.5D0, 34.8D0, 30.6D0,  7.0D0, 19.8D0, 92.5D0/
      DATA   Y( 79),  Y( 80),  Y( 81),  Y( 82),  Y( 83),  Y( 84)
     +    /   154.4D0,125.9D0, 84.8D0, 68.1D0, 38.5D0, 22.8D0/
      DATA   Y( 85),  Y( 86),  Y( 87),  Y( 88),  Y( 89),  Y( 90)
     +    /    10.2D0, 24.1D0, 82.9D0,132.0D0,130.9D0,118.1D0/
      DATA   Y( 91),  Y( 92),  Y( 93),  Y( 94),  Y( 95),  Y( 96)
     +    /    89.9D0, 66.6D0, 60.0D0, 46.9D0, 41.0D0, 21.3D0/
      DATA   Y( 97),  Y( 98),  Y( 99),  Y(100),  Y(101),  Y(102)
     +    /    16.0D0,  6.4D0,  4.1D0,  6.8D0, 14.5D0, 34.0D0/
      DATA   Y(103),  Y(104),  Y(105),  Y(106),  Y(107),  Y(108)
     +    /    45.0D0, 43.1D0, 47.5D0, 42.2D0, 28.1D0, 10.1D0/
      DATA   Y(109),  Y(110),  Y(111),  Y(112),  Y(113),  Y(114)
     +    /     8.1D0,  2.5D0,  0.0D0,  1.4D0,  5.0D0, 12.2D0/
      DATA   Y(115),  Y(116),  Y(117),  Y(118),  Y(119),  Y(120)
     +    /    13.9D0, 35.4D0, 45.8D0, 41.1D0, 30.1D0, 23.9D0/
      DATA   Y(121),  Y(122),  Y(123),  Y(124),  Y(125),  Y(126)
     +    /    15.6D0,  6.6D0,  4.0D0,  1.8D0,  8.5D0, 16.6D0/
      DATA   Y(127),  Y(128),  Y(129),  Y(130),  Y(131),  Y(132)
     +    /    36.3D0, 49.6D0, 64.2D0, 67.0D0, 70.9D0, 47.8D0/
      DATA   Y(133),  Y(134),  Y(135),  Y(136),  Y(137),  Y(138)
     +    /    27.5D0,  8.5D0, 13.2D0, 56.9D0,121.5D0,138.3D0/
      DATA   Y(139),  Y(140),  Y(141),  Y(142),  Y(143),  Y(144)
     +    /   103.2D0, 85.7D0, 64.6D0, 36.7D0, 24.2D0, 10.7D0/
      DATA   Y(145),  Y(146),  Y(147),  Y(148),  Y(149),  Y(150)
     +    /    15.0D0, 40.1D0, 61.5D0, 98.5D0,124.7D0, 96.3D0/
      DATA   Y(151),  Y(152),  Y(153),  Y(154),  Y(155),  Y(156)
     +    /    66.6D0, 64.5D0, 54.1D0, 39.0D0, 20.6D0,  6.7D0/
      DATA   Y(157),  Y(158),  Y(159),  Y(160),  Y(161),  Y(162)
     +    /     4.3D0, 22.7D0, 54.8D0, 93.8D0, 95.8D0, 77.2D0/
      DATA   Y(163),  Y(164),  Y(165),  Y(166),  Y(167),  Y(168)
     +    /    59.1D0, 44.0D0, 47.0D0, 30.5D0, 16.3D0,  7.3D0/
      DATA   Y(169),  Y(170),  Y(171),  Y(172),  Y(173),  Y(174)
     +    /    37.6D0, 74.0D0,139.0D0,111.2D0,101.6D0, 66.2D0/
      DATA   Y(175),  Y(176),  Y(177),  Y(178),  Y(179),  Y(180)
     +    /    44.7D0, 17.0D0, 11.3D0, 12.4D0,  3.4D0,  6.0D0/
      DATA   Y(181),  Y(182),  Y(183),  Y(184),  Y(185),  Y(186)
     +    /    32.3D0, 54.3D0, 59.7D0, 63.7D0, 63.5D0, 52.2D0/
      DATA   Y(187),  Y(188),  Y(189),  Y(190),  Y(191),  Y(192)
     +    /    25.4D0, 13.1D0,  6.8D0,  6.3D0,  7.1D0, 35.6D0/
      DATA   Y(193),  Y(194),  Y(195),  Y(196),  Y(197),  Y(198)
     +    /    73.0D0, 85.1D0, 78.0D0, 64.0D0, 41.8D0, 26.2D0/
      DATA   Y(199),  Y(200),  Y(201),  Y(202),  Y(203),  Y(204)
     +    /    26.7D0, 12.1D0,  9.5D0,  2.7D0,  5.0D0, 24.4D0/
      DATA   Y(205),  Y(206),  Y(207),  Y(208),  Y(209),  Y(210)
     +    /    42.0D0, 63.5D0, 53.8D0, 62.0D0, 48.5D0, 43.9D0/
      DATA   Y(211),  Y(212),  Y(213),  Y(214),  Y(215),  Y(216)
     +    /    18.6D0,  5.7D0,  3.6D0,  1.4D0,  9.6D0, 47.4D0/
      DATA   Y(217),  Y(218),  Y(219),  Y(220),  Y(221),  Y(222)
     +    /    57.1D0,103.9D0, 80.6D0, 63.6D0, 37.6D0, 26.1D0/
      DATA   Y(223),  Y(224),  Y(225),  Y(226),  Y(227),  Y(228)
     +    /    14.2D0,  5.8D0, 16.7D0, 44.3D0, 63.9D0, 69.0D0/
      DATA   Y(229),  Y(230),  Y(231),  Y(232),  Y(233),  Y(234)
     +    /    77.8D0, 64.9D0, 35.7D0, 21.2D0, 11.1D0,  5.7D0/
      DATA   Y(235),  Y(236),  Y(237),  Y(238),  Y(239),  Y(240)
     +    /     8.7D0, 36.1D0, 79.7D0,114.4D0,109.6D0, 88.8D0/
      DATA   Y(241),  Y(242),  Y(243),  Y(244),  Y(245),  Y(246)
     +    /    67.8D0, 47.5D0, 30.6D0, 16.3D0,  9.6D0, 33.2D0/
      DATA   Y(247),  Y(248),  Y(249),  Y(250),  Y(251),  Y(252)
     +    /    92.6D0,151.6D0,136.3D0,134.7D0, 83.9D0, 69.4D0/
      DATA   Y(253),  Y(254),  Y(255),  Y(256),  Y(257),  Y(258)
     +    /    31.5D0, 13.9D0,  4.4D0, 38.0D0,141.7D0,190.2D0/
      DATA   Y(259),  Y(260),  Y(261)
     +    /   184.8D0,159.0D0,112.3D0/
C
C     DEFINE CONSTANTS
C
      CALL IPRINT(IPRT)
      N = 261
      FD = 1.0/11.0
      FC = 1.0/22.0
      K = 41
C
C     WRITE HEADER
C
      WRITE(IPRT, 1000)
C
C     RUN SIMPLE TEST OF DIF
C
      WRITE(IPRT, 1100)
      CALL DEMOD (Y, N, FD, FC, K, LDSTAK)
      WRITE (IPRT,2000) IERR
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT ('1*CH11')
 1100 FORMAT (' SIMPLE TEST OF DEMOD')
 2000 FORMAT (/' THE VALUE OF IERR IS ', I4)
      END
