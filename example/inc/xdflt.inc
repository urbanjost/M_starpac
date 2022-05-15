subroutine xdflt ( lds ,typ)
use M_starpac, only : arflt,dif,difc,difm,difmc,gfarf,gfarfs,gfslf,gfslfs, hipass,hpcoef,iprint,lopass,lpcoef,maflt,sample,slflt
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XDFLT tests time series digital filtering and complex demodulation routines.
!
!  Discussion:
!
!    series y is the wolf sunspot data from 1700 to 1960 as
!    tabulated by waldmeier
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
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     real fc
!        the cutoff frequency used for the low pass filter.
!     real fmax, fmin
!        the minimum and maximum frequency for which the gain
!        function is to be estimated.
!     real freq(101)
!        the vector of frequencies at which the gain function
!        has been estimated.
!     real gain(101)
!        the vector in which the gain function estimates are
!        stored.
!     real hhp(50)
!        the array in which the -ideal- high pass filter coefficients
!        will be returned.
!     real hlp(50)
!        the array in which the input low pass filter coefficients
!        are stored.
!     integer i
!        an indexing variable.
!     integer iar
!        the number of filter coefficients.
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr == 1, errors have been detected.
!     integer iod(10)
!        the order of each of the difference factors.
!     integer itest
!        the number of the test being performed
!     integer k
!        the number of terms in the symetric linear filter.
!     integer lds, ldstak
!        the length of the array dstak.
!     integer lphi
!        the length of the vector phi.
!     integer n
!        the number of observations in the input series.
!     integer nd(10)
!        the array containing the number of times the difference
!        factors are to be applied.
!     integer nfac
!        the number of difference factors.
!     integer nf
!        the number of frequencies at which the gain function
!        is to be estimated.
!     integer nprt
!        a code used to specify the type of plot, where if
!        nprt = 0 the plot is suppressed
!        nprt = 1 the plot is decibels/linear
!        nprt = 2 the plot is log/linear
!     integer nyf
!        the number of values in the filtered series.
!     integer nys
!        the number of values in the sampled series.
!     real phas(300)
!        the array in which the primary phase estimates are returned.
!     real phi(50)
!        the vector containing the filter coefficients.
!     real y(300)
!        the vector containing the observed time series.
!     real yf(300)
!        the vector in which the filtered series is returned.
!     real yfmiss
!        the missing value code used in the filtered series.
!     real ymiss
!        the missing value code used in the input time series.
!     real ys(300)
!        the array containing the sampled series.
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
     fc,fmax,fmin,yfmiss,ymiss
  integer &
     i,iar,itest,k,ldstak,lphi,n,nf,nfac,nprt,nyf,nys
!
!  local arrays
  real &
     freq(101),gain(101),hhp(50),hlp(50),phas(300),phi(50),y(300), &
     yf(300),ys(300)
  integer &
     iod(10),nd(10)
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
  nprt = 2
  fc = 1.0e0/22.0e0
  nf = 101
  fmin = 0.0e0
  fmax = 0.2e0
  lphi = 50
  nfac = 1
  nd(1) = 1
  iod(1) = 1
  iar = 1
  phi(1) = 0.6e0
  k = 41
  ymiss = 11.0e0
!
!  test of lpcoef
!
   10 write ( *, 1001)
  call lpcoef (fc, k, hlp)
  write ( *, 1002) ierr
!
!  print storage from lpcoef
!
  if (ierr==0) write ( *, 1004) hlp(1:k)
!
!  test of lopass
!
  write ( *, 1007)
  call lopass (y, n, fc, k, hlp, yf, nyf)
  write ( *, 1002) ierr
!
!  print storage from lopass
!
  if (ierr==0) then
    write ( *, 1004) (hlp(i), i = 1, k)
    write ( *, 1004) yf(1:nyf)
  end if
!
!  test of hipass
!
  write ( *, 1008)
  call hipass (y, n, fc, k, hhp, yf, nyf)
  write ( *, 1002) ierr
!
!  print storage from hipass
!
  if (ierr==0) then
    write ( *, 1004) (hhp(i), i = 1, k)
    write ( *, 1004) yf(1:nyf)
  end if
!
!  test of hpcoef
!
   20 write ( *, 1003)
  call hpcoef (hlp, k, hhp)
  write ( *, 1002) ierr
!
!  print storage from hpcoef
!
  if (ierr==0) write ( *, 1004) (hhp(i), i = 1, k)
!
!  test of maflt
!
  write ( *, 1020)
  call maflt (y, n, k, yf, nyf)
  write ( *, 1002) ierr
!
!  print storage from maflt
!
  if (ierr==0) write ( *, 1004) yf(1:nyf)
!
!  test of slflt
!
  write ( *, 1005)
  call slflt (y, n, k, hlp, yf, nyf)
  write ( *, 1002) ierr
!
!  print storage from slflt
!
  if (ierr==0) write ( *, 1004) yf(1:nyf)
!
!  test of sample
!
  write ( *, 1006)
  call sample (yf, n, k, ys, nys)
  write ( *, 1002) ierr
!
!  print storage from sample
!
  if (ierr==0) write ( *, 1004) yf(1:nys)
!
!  test of arflt
!
  write ( *, 1009)
  call arflt (y, n,  iar, phi, yf, nyf)
  write ( *, 1002) ierr
!
!  print storage from arflt
!
  if (ierr==0) write ( *, 1004) yf(1:nyf)
!
!  test of dif
!
  write ( *, 1015)
  call dif (y, n, yf, nyf)
  write ( *, 1002) ierr
!
!  print storage from dif
!
  if (ierr==0) write ( *, 1004) yf(1:nyf)
!
!  test of difm
!
  write ( *, 1018)
  call difm (y, ymiss, n, yf, yfmiss, nyf)
  write ( *, 1002) ierr
!
!  print storage from difm
!
  if (ierr==0) then
     write ( *, 1004) yf(1:nyf)
     write ( *, 1004) yfmiss
  end if
!
!  test of gfslf
!
  write ( *, 1011)
  call gfslf (hlp, k)
  write ( *, 1002) ierr
!
!  test of gfarf
!
  write ( *, 1013)
  call gfarf (phi, iar)
  write ( *, 1002) ierr
!
!  test of difc
!
   30 write ( *, 1010)
  call difc (y, n, nfac, nd, iod, iar, phi, lphi, yf, nyf, ldstak)
  write ( *, 1002) ierr
!
!  print storage from difc
!
  if (ierr==0) then
    write ( *, 1004) (phi(i), i = 1, k)
    write ( *, 1004) yf(1:nyf)
  end if
!
!  test of difmc
!
  write ( *, 1019)
  call difmc (y, ymiss, n, nfac, nd, iod, iar, phi, lphi, yf, &
     yfmiss, nyf, ldstak)
  write ( *, 1002) ierr
!
!  print storage from difmc
!
  if (ierr==0) then
    write ( *, 1004) (phi(i), i = 1, k)
    write ( *, 1004) yf(1:nyf)
    write ( *, 1004) yfmiss
  end if
!
!  test of gfslfs
!
  write ( *, 1012)
  call gfslfs (hlp, k, nf, fmin, fmax, gain, freq, nprt, ldstak)
  write ( *, 1002) ierr
!
!  print storage from gfslfs
!
  if (ierr==0) then
    write ( *, 1004) (gain(i), i = 1, nf)
    write ( *, 1004) (freq(i), i = 1, nf)
  end if
!
!  test of gfarfs
!
  write ( *, 1014)
  call gfarfs (phi, iar, nf, fmin, fmax, gain, phas, freq, nprt, &
     ldstak)
  write ( *, 1002) ierr
!
!  print storage from gfarfs
!
  if (ierr==0) then
    write ( *, 1004) (gain(i), i = 1, nf)
    write ( *, 1004) (phas(i), i = 1, nf)
    write ( *, 1004) (freq(i), i = 1, nf)
  end if

  go to (100, 200, 300, 400), itest
!
!  test special cases
!
  100 itest = itest + 1
!
!  test of gfslfs
!
  fmin = 0.4e0
  fmax = 0.1e0
  nprt = 1
  write ( *, 1012)
  call gfslfs (hlp, k, nf, fmin, fmax, gain, freq, nprt, ldstak)
  write ( *, 1002) ierr
!
!  print storage from gfslfs
!
  if (ierr==0) then
    write ( *, 1004) (gain(i), i = 1, nf)
    write ( *, 1004) (freq(i), i = 1, nf)
  end if
!
!  test of gfarfs
!
  nprt = -1
  write ( *, 1014)
  call gfarfs (phi, iar, nf, fmin, fmax, gain, phas, freq, nprt, &
     ldstak)
  write ( *, 1002) ierr
!
!  print storage from gfarfs
!
  if (ierr==0) then
    write ( *, 1004) (gain(i), i = 1, nf)
    write ( *, 1004) (phas(i), i = 1, nf)
    write ( *, 1004) (freq(i), i = 1, nf)
  end if
!
!  test minimum problem size
!
  n = 3
  k = 1
  nprt = -1
  iar = 1
  nf = 1
  go to 20
!
!  test error conditions
!
  200 itest = itest + 1
  n = -5
  fc = 1.0e0
  nf = 0
  lphi = 0
  nfac = 1
  nd(1) = -1
  iod(1) = -1
  iar = 0
  k = -1
  go to 10
!
!  test ldstak
!
  300 itest = itest + 1
  n = 261
  nprt = 2
  fc = 1.0e0/22.0e0
  nf = 101
  fmin = 0.0e0
  fmax = 0.2e0
  lphi = 50
  nfac = 1
  nd(1) = 1
  iod(1) = 1
  iar = 1
  phi(1) = 0.6e0
  k = 41
  ymiss = 11.0e0
  ldstak = 0
  go to 30
!
  400 return

 1001 format ('test of lpcoef')
 1002 format (/' ierr is ', i5)
 1003 format ('test of hpcoef')
 1004 format (10e10.3)
 1005 format ('test of slflt')
 1006 format ('test of sample')
 1007 format ('test of lopass')
 1008 format ('test of hipass')
 1009 format ('test of arflt')
 1010 format ('test of difc')
 1011 format ('test of gfslf')
 1012 format ('test of gfslfs')
 1013 format ('test of gfarf')
 1014 format ('test of gfarfs')
 1015 format ('test of dif')
 1018 format ('test of difm')
 1019 format ('test of difmc')
 1020 format ('test of maflt')
end