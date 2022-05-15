subroutine xdemod ( lds ,typ)
use M_starpac, only : demod,demods,iprint
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XDEMOD tests the time series complex demodulation routines.
!
!  Discussion:
!
!    Series Y is the Wolf sunspot data from 1700 to 1960 as
!    tabulated by Waldmeier.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Reference:
!
!    Max Waldmeier,
!    The Sunspot-Activity in the Years 1610-1960,
!    Shulthess, Zurich, 1961.
!
!  Parameters:
!
!     real ampl(300)
!        the array in which the amplitudes are stored.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     real fc
!        the cutoff frequency used for the low pass filter.
!     real fd
!        the demodulation frequency.
!     integer i
!        an indexing variable.
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr == 1, errors have been detected.
!     integer itest
!        the number of the test being run
!     integer k
!        the number of terms in the symetric linear filter.
!     integer lds, ldstak
!        the length of the array dstak.
!     integer n
!        the number of observations in the input series.
!     integer ndem
!        the number of values in the demodulated series, i. e., it
!        is the number of values in the amplitude and phase arrays.
!     integer nprt
!        a code used to specify the type of plot, where if
!        nprt == 0 the plot is suppressed
!        nprt .ne. 1 the plot is provided
!     real phas(300)
!        the array in which the primary phase estimates are returned.
!     real y(300)
!        the vector containing the observed time series.
!

  integer &
     lds
!
!  scalars in common
  integer &
     ierr
!
!  arrays in common
  double precision dstak(12)
!
!  local scalars
  real &
     fc,fd
  integer &
     i,itest,k,ldstak,n,ndem,nprt
!
!  local arrays
  real &
     ampl(300),phas(300),y(300)
!
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr

  data   y(  1),  y(  2),  y(  3),  y(  4),  y(  5),  y(  6) &
      /     5.0e0, 11.0e0, 16.0e0, 23.0e0, 36.0e0, 58.0e0/
  data   y(  7),  y(  8),  y(  9),  y( 10),  y( 11),  y( 12) &
      /    29.0e0, 20.0e0, 10.0e0,  8.0e0,  3.0e0,  0.0e0/
  data   y( 13),  y( 14),  y( 15),  y( 16),  y( 17),  y( 18) &
      /     0.0e0, 2.0e0, 11.0e0, 27.0e0, 47.0e0, 63.0e0/
  data   y( 19),  y( 20),  y( 21),  y( 22),  y( 23),  y( 24) &
      /    60.0e0, 39.0e0, 28.0e0, 26.0e0, 22.0e0, 11.0e0/
  data   y( 25),  y( 26),  y( 27),  y( 28),  y( 29),  y( 30) &
      /    21.0e0, 40.0e0, 78.0e0,122.0e0,103.0e0, 73.0e0/
  data   y( 31),  y( 32),  y( 33),  y( 34),  y( 35),  y( 36) &
      /    47.0e0, 35.0e0, 11.0e0,  5.0e0, 16.0e0, 34.0e0/
  data   y( 37),  y( 38),  y( 39),  y( 40),  y( 41),  y( 42) &
      /    70.0e0, 81.0e0,111.0e0,101.0e0, 73.0e0, 40.0e0/
  data   y( 43),  y( 44),  y( 45),  y( 46),  y( 47),  y( 48) &
      /    20.0e0, 16.0e0,  5.0e0, 11.0e0, 22.0e0, 40.0e0/
  data   y( 49),  y( 50),  y( 51),  y( 52),  y( 53),  y( 54) &
      /    60.0e0, 80.9e0, 83.4e0, 47.7e0, 47.8e0, 30.7e0/
  data   y( 55),  y( 56),  y( 57),  y( 58),  y( 59),  y( 60) &
      /    12.2e0,  9.6e0, 10.2e0, 32.4e0, 47.6e0, 54.0e0/
  data   y( 61),  y( 62),  y( 63),  y( 64),  y( 65),  y( 66) &
      /    62.9e0, 85.9e0, 61.2e0, 45.1e0, 36.4e0, 20.9e0/
  data   y( 67),  y( 68),  y( 69),  y( 70),  y( 71),  y( 72) &
      /    11.4e0, 37.8e0, 69.8e0,106.1e0,100.8e0, 81.6e0/
  data   y( 73),  y( 74),  y( 75),  y( 76),  y( 77),  y( 78) &
      /    66.5e0, 34.8e0, 30.6e0,  7.0e0, 19.8e0, 92.5e0/
  data   y( 79),  y( 80),  y( 81),  y( 82),  y( 83),  y( 84) &
      /   154.4e0,125.9e0, 84.8e0, 68.1e0, 38.5e0, 22.8e0/
  data   y( 85),  y( 86),  y( 87),  y( 88),  y( 89),  y( 90) &
      /    10.2e0, 24.1e0, 82.9e0,132.0e0,130.9e0,118.1e0/
  data   y( 91),  y( 92),  y( 93),  y( 94),  y( 95),  y( 96) &
      /    89.9e0, 66.6e0, 60.0e0, 46.9e0, 41.0e0, 21.3e0/
  data   y( 97),  y( 98),  y( 99),  y(100),  y(101),  y(102) &
      /    16.0e0,  6.4e0,  4.1e0,  6.8e0, 14.5e0, 34.0e0/
  data   y(103),  y(104),  y(105),  y(106),  y(107),  y(108) &
      /    45.0e0, 43.1e0, 47.5e0, 42.2e0, 28.1e0, 10.1e0/
  data   y(109),  y(110),  y(111),  y(112),  y(113),  y(114) &
      /     8.1e0,  2.5e0,  0.0e0,  1.4e0,  5.0e0, 12.2e0/
  data   y(115),  y(116),  y(117),  y(118),  y(119),  y(120) &
      /    13.9e0, 35.4e0, 45.8e0, 41.1e0, 30.1e0, 23.9e0/
  data   y(121),  y(122),  y(123),  y(124),  y(125),  y(126) &
      /    15.6e0,  6.6e0,  4.0e0,  1.8e0,  8.5e0, 16.6e0/
  data   y(127),  y(128),  y(129),  y(130),  y(131),  y(132) &
      /    36.3e0, 49.6e0, 64.2e0, 67.0e0, 70.9e0, 47.8e0/
  data   y(133),  y(134),  y(135),  y(136),  y(137),  y(138) &
      /    27.5e0,  8.5e0, 13.2e0, 56.9e0,121.5e0,138.3e0/
  data   y(139),  y(140),  y(141),  y(142),  y(143),  y(144) &
      /   103.2e0, 85.7e0, 64.6e0, 36.7e0, 24.2e0, 10.7e0/
  data   y(145),  y(146),  y(147),  y(148),  y(149),  y(150) &
      /    15.0e0, 40.1e0, 61.5e0, 98.5e0,124.7e0, 96.3e0/
  data   y(151),  y(152),  y(153),  y(154),  y(155),  y(156) &
      /    66.6e0, 64.5e0, 54.1e0, 39.0e0, 20.6e0,  6.7e0/
  data   y(157),  y(158),  y(159),  y(160),  y(161),  y(162) &
      /     4.3e0, 22.7e0, 54.8e0, 93.8e0, 95.8e0, 77.2e0/
  data   y(163),  y(164),  y(165),  y(166),  y(167),  y(168) &
      /    59.1e0, 44.0e0, 47.0e0, 30.5e0, 16.3e0,  7.3e0/
  data   y(169),  y(170),  y(171),  y(172),  y(173),  y(174) &
      /    37.6e0, 74.0e0,139.0e0,111.2e0,101.6e0, 66.2e0/
  data   y(175),  y(176),  y(177),  y(178),  y(179),  y(180) &
      /    44.7e0, 17.0e0, 11.3e0, 12.4e0,  3.4e0,  6.0e0/
  data   y(181),  y(182),  y(183),  y(184),  y(185),  y(186) &
      /    32.3e0, 54.3e0, 59.7e0, 63.7e0, 63.5e0, 52.2e0/
  data   y(187),  y(188),  y(189),  y(190),  y(191),  y(192) &
      /    25.4e0, 13.1e0,  6.8e0,  6.3e0,  7.1e0, 35.6e0/
  data   y(193),  y(194),  y(195),  y(196),  y(197),  y(198) &
      /    73.0e0, 85.1e0, 78.0e0, 64.0e0, 41.8e0, 26.2e0/
  data   y(199),  y(200),  y(201),  y(202),  y(203),  y(204) &
      /    26.7e0, 12.1e0,  9.5e0,  2.7e0,  5.0e0, 24.4e0/
  data   y(205),  y(206),  y(207),  y(208),  y(209),  y(210) &
      /    42.0e0, 63.5e0, 53.8e0, 62.0e0, 48.5e0, 43.9e0/
  data   y(211),  y(212),  y(213),  y(214),  y(215),  y(216) &
      /    18.6e0,  5.7e0,  3.6e0,  1.4e0,  9.6e0, 47.4e0/
  data   y(217),  y(218),  y(219),  y(220),  y(221),  y(222) &
      /    57.1e0,103.9e0, 80.6e0, 63.6e0, 37.6e0, 26.1e0/
  data   y(223),  y(224),  y(225),  y(226),  y(227),  y(228) &
      /    14.2e0,  5.8e0, 16.7e0, 44.3e0, 63.9e0, 69.0e0/
  data   y(229),  y(230),  y(231),  y(232),  y(233),  y(234) &
      /    77.8e0, 64.9e0, 35.7e0, 21.2e0, 11.1e0,  5.7e0/
  data   y(235),  y(236),  y(237),  y(238),  y(239),  y(240) &
      /     8.7e0, 36.1e0, 79.7e0,114.4e0,109.6e0, 88.8e0/
  data   y(241),  y(242),  y(243),  y(244),  y(245),  y(246) &
      /    67.8e0, 47.5e0, 30.6e0, 16.3e0,  9.6e0, 33.2e0/
  data   y(247),  y(248),  y(249),  y(250),  y(251),  y(252) &
      /    92.6e0,151.6e0,136.3e0,134.7e0, 83.9e0, 69.4e0/
  data   y(253),  y(254),  y(255),  y(256),  y(257),  y(258) &
      /    31.5e0, 13.9e0,  4.4e0, 38.0e0,141.7e0,190.2e0/
  data   y(259),  y(260),  y(261) &
      /   184.8e0,159.0e0,112.3e0/

  itest = 1
  ldstak = lds

  n = 261
  nprt = 1
  fd = 1.0e0/11.0e0
  fc = 1.0e0/22.0e0
  k = 41
!
!  test of demod
!
5 continue

  write ( *, 1016)
  call demod (y, n, fd, fc, k, ldstak)
  write ( *, 1002) ierr
!
!  test of demods
!
  write ( *, 1017)
  call demods (y, n, fd, fc, k, ampl, phas, ndem, nprt, ldstak)
  write ( *, 1002) ierr
!
!  print storage from demods
!
  if (ierr==0) then
    write ( *, 1004) (ampl(i), i = 1, ndem)
    write ( *, 1004) (phas(i), i = 1, ndem)
  end if

  go to (100, 200, 300), itest
!
!  test minimum problem specifications
!
  100 continue

  itest = itest + 1
  n = 17
  k = 15
  nprt = -1
  go to 5
!
!  test error conditions
!
  200 continue

  itest = itest + 1
  n = 0
  fd = 0.5e0
  fc = 0.3e0
  k = 1
  go to 5

  300 return

 1002 format (' ierr is', i5)
 1004 format (10f10.5)
 1016 format ('test of demod')
 1017 format ('test of demods')
end
