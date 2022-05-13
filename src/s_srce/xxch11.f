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
      REAL
     +   FC,FD
      INTEGER
     +   IPRT,K,N
C
C  LOCAL ARRAYS
      REAL
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
C     REAL FC
C        THE CUTOFF FREQUENCY USED FOR THE LOW PASS FILTER.
C     REAL FD
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
C     REAL Y(300)
C        THE LOG OF THE AIRLINE DATA.
C
C
      DATA   Y(  1),  Y(  2),  Y(  3),  Y(  4),  Y(  5),  Y(  6)
     +    /     5.0E0, 11.0E0, 16.0E0, 23.0E0, 36.0E0, 58.0E0/
      DATA   Y(  7),  Y(  8),  Y(  9),  Y( 10),  Y( 11),  Y( 12)
     +    /    29.0E0, 20.0E0, 10.0E0,  8.0E0,  3.0E0,  0.0E0/
      DATA   Y( 13),  Y( 14),  Y( 15),  Y( 16),  Y( 17),  Y( 18)
     +    /     0.0E0, 2.0E0, 11.0E0, 27.0E0, 47.0E0, 63.0E0/
      DATA   Y( 19),  Y( 20),  Y( 21),  Y( 22),  Y( 23),  Y( 24)
     +    /    60.0E0, 39.0E0, 28.0E0, 26.0E0, 22.0E0, 11.0E0/
      DATA   Y( 25),  Y( 26),  Y( 27),  Y( 28),  Y( 29),  Y( 30)
     +    /    21.0E0, 40.0E0, 78.0E0,122.0E0,103.0E0, 73.0E0/
      DATA   Y( 31),  Y( 32),  Y( 33),  Y( 34),  Y( 35),  Y( 36)
     +    /    47.0E0, 35.0E0, 11.0E0,  5.0E0, 16.0E0, 34.0E0/
      DATA   Y( 37),  Y( 38),  Y( 39),  Y( 40),  Y( 41),  Y( 42)
     +    /    70.0E0, 81.0E0,111.0E0,101.0E0, 73.0E0, 40.0E0/
      DATA   Y( 43),  Y( 44),  Y( 45),  Y( 46),  Y( 47),  Y( 48)
     +    /    20.0E0, 16.0E0,  5.0E0, 11.0E0, 22.0E0, 40.0E0/
      DATA   Y( 49),  Y( 50),  Y( 51),  Y( 52),  Y( 53),  Y( 54)
     +    /    60.0E0, 80.9E0, 83.4E0, 47.7E0, 47.8E0, 30.7E0/
      DATA   Y( 55),  Y( 56),  Y( 57),  Y( 58),  Y( 59),  Y( 60)
     +    /    12.2E0,  9.6E0, 10.2E0, 32.4E0, 47.6E0, 54.0E0/
      DATA   Y( 61),  Y( 62),  Y( 63),  Y( 64),  Y( 65),  Y( 66)
     +    /    62.9E0, 85.9E0, 61.2E0, 45.1E0, 36.4E0, 20.9E0/
      DATA   Y( 67),  Y( 68),  Y( 69),  Y( 70),  Y( 71),  Y( 72)
     +    /    11.4E0, 37.8E0, 69.8E0,106.1E0,100.8E0, 81.6E0/
      DATA   Y( 73),  Y( 74),  Y( 75),  Y( 76),  Y( 77),  Y( 78)
     +    /    66.5E0, 34.8E0, 30.6E0,  7.0E0, 19.8E0, 92.5E0/
      DATA   Y( 79),  Y( 80),  Y( 81),  Y( 82),  Y( 83),  Y( 84)
     +    /   154.4E0,125.9E0, 84.8E0, 68.1E0, 38.5E0, 22.8E0/
      DATA   Y( 85),  Y( 86),  Y( 87),  Y( 88),  Y( 89),  Y( 90)
     +    /    10.2E0, 24.1E0, 82.9E0,132.0E0,130.9E0,118.1E0/
      DATA   Y( 91),  Y( 92),  Y( 93),  Y( 94),  Y( 95),  Y( 96)
     +    /    89.9E0, 66.6E0, 60.0E0, 46.9E0, 41.0E0, 21.3E0/
      DATA   Y( 97),  Y( 98),  Y( 99),  Y(100),  Y(101),  Y(102)
     +    /    16.0E0,  6.4E0,  4.1E0,  6.8E0, 14.5E0, 34.0E0/
      DATA   Y(103),  Y(104),  Y(105),  Y(106),  Y(107),  Y(108)
     +    /    45.0E0, 43.1E0, 47.5E0, 42.2E0, 28.1E0, 10.1E0/
      DATA   Y(109),  Y(110),  Y(111),  Y(112),  Y(113),  Y(114)
     +    /     8.1E0,  2.5E0,  0.0E0,  1.4E0,  5.0E0, 12.2E0/
      DATA   Y(115),  Y(116),  Y(117),  Y(118),  Y(119),  Y(120)
     +    /    13.9E0, 35.4E0, 45.8E0, 41.1E0, 30.1E0, 23.9E0/
      DATA   Y(121),  Y(122),  Y(123),  Y(124),  Y(125),  Y(126)
     +    /    15.6E0,  6.6E0,  4.0E0,  1.8E0,  8.5E0, 16.6E0/
      DATA   Y(127),  Y(128),  Y(129),  Y(130),  Y(131),  Y(132)
     +    /    36.3E0, 49.6E0, 64.2E0, 67.0E0, 70.9E0, 47.8E0/
      DATA   Y(133),  Y(134),  Y(135),  Y(136),  Y(137),  Y(138)
     +    /    27.5E0,  8.5E0, 13.2E0, 56.9E0,121.5E0,138.3E0/
      DATA   Y(139),  Y(140),  Y(141),  Y(142),  Y(143),  Y(144)
     +    /   103.2E0, 85.7E0, 64.6E0, 36.7E0, 24.2E0, 10.7E0/
      DATA   Y(145),  Y(146),  Y(147),  Y(148),  Y(149),  Y(150)
     +    /    15.0E0, 40.1E0, 61.5E0, 98.5E0,124.7E0, 96.3E0/
      DATA   Y(151),  Y(152),  Y(153),  Y(154),  Y(155),  Y(156)
     +    /    66.6E0, 64.5E0, 54.1E0, 39.0E0, 20.6E0,  6.7E0/
      DATA   Y(157),  Y(158),  Y(159),  Y(160),  Y(161),  Y(162)
     +    /     4.3E0, 22.7E0, 54.8E0, 93.8E0, 95.8E0, 77.2E0/
      DATA   Y(163),  Y(164),  Y(165),  Y(166),  Y(167),  Y(168)
     +    /    59.1E0, 44.0E0, 47.0E0, 30.5E0, 16.3E0,  7.3E0/
      DATA   Y(169),  Y(170),  Y(171),  Y(172),  Y(173),  Y(174)
     +    /    37.6E0, 74.0E0,139.0E0,111.2E0,101.6E0, 66.2E0/
      DATA   Y(175),  Y(176),  Y(177),  Y(178),  Y(179),  Y(180)
     +    /    44.7E0, 17.0E0, 11.3E0, 12.4E0,  3.4E0,  6.0E0/
      DATA   Y(181),  Y(182),  Y(183),  Y(184),  Y(185),  Y(186)
     +    /    32.3E0, 54.3E0, 59.7E0, 63.7E0, 63.5E0, 52.2E0/
      DATA   Y(187),  Y(188),  Y(189),  Y(190),  Y(191),  Y(192)
     +    /    25.4E0, 13.1E0,  6.8E0,  6.3E0,  7.1E0, 35.6E0/
      DATA   Y(193),  Y(194),  Y(195),  Y(196),  Y(197),  Y(198)
     +    /    73.0E0, 85.1E0, 78.0E0, 64.0E0, 41.8E0, 26.2E0/
      DATA   Y(199),  Y(200),  Y(201),  Y(202),  Y(203),  Y(204)
     +    /    26.7E0, 12.1E0,  9.5E0,  2.7E0,  5.0E0, 24.4E0/
      DATA   Y(205),  Y(206),  Y(207),  Y(208),  Y(209),  Y(210)
     +    /    42.0E0, 63.5E0, 53.8E0, 62.0E0, 48.5E0, 43.9E0/
      DATA   Y(211),  Y(212),  Y(213),  Y(214),  Y(215),  Y(216)
     +    /    18.6E0,  5.7E0,  3.6E0,  1.4E0,  9.6E0, 47.4E0/
      DATA   Y(217),  Y(218),  Y(219),  Y(220),  Y(221),  Y(222)
     +    /    57.1E0,103.9E0, 80.6E0, 63.6E0, 37.6E0, 26.1E0/
      DATA   Y(223),  Y(224),  Y(225),  Y(226),  Y(227),  Y(228)
     +    /    14.2E0,  5.8E0, 16.7E0, 44.3E0, 63.9E0, 69.0E0/
      DATA   Y(229),  Y(230),  Y(231),  Y(232),  Y(233),  Y(234)
     +    /    77.8E0, 64.9E0, 35.7E0, 21.2E0, 11.1E0,  5.7E0/
      DATA   Y(235),  Y(236),  Y(237),  Y(238),  Y(239),  Y(240)
     +    /     8.7E0, 36.1E0, 79.7E0,114.4E0,109.6E0, 88.8E0/
      DATA   Y(241),  Y(242),  Y(243),  Y(244),  Y(245),  Y(246)
     +    /    67.8E0, 47.5E0, 30.6E0, 16.3E0,  9.6E0, 33.2E0/
      DATA   Y(247),  Y(248),  Y(249),  Y(250),  Y(251),  Y(252)
     +    /    92.6E0,151.6E0,136.3E0,134.7E0, 83.9E0, 69.4E0/
      DATA   Y(253),  Y(254),  Y(255),  Y(256),  Y(257),  Y(258)
     +    /    31.5E0, 13.9E0,  4.4E0, 38.0E0,141.7E0,190.2E0/
      DATA   Y(259),  Y(260),  Y(261)
     +    /   184.8E0,159.0E0,112.3E0/
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
