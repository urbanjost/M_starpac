subroutine xpgm ( ldstak ,typ)
use M_starpac, only : &
     center,fftlen,fftr,ipgm,ipgmp,ipgmps,ipgms,iprint,mdflt, &
     pgm,pgms,ppl,scopy,taper
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XPGM tests the time series periodogram subroutines.
!
!  Discussion:
!
!    series y is the first 50 values of the series listed on page
!    318 of jenkins and watts.  the spectrum of this series is shown
!    for various bandwidth on page 270 of jenkins and watts.
!
!    series z is the wolf sunspot numbers from 1700 to 1960 as
!    tabulated by waldmeier.  the raw and smoothed periodograms of
!    tapered series are shown on pages 95 and 176, respectively, of
!    bloomfield.
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
!    Gwilym Jenkins, Donald Watts,
!    Spectral Analysis and its Applications,
!    Holden-Day 1968.
!
!    Max Waldmeier,
!    The Sunspot-Activity in the Years 1610-1960,
!    Shulthess, Zurich, 1961.
!
!  Parameters:
!
!     real ab(600)
!        the vector of the nf real and imaginary components of the
!        fourier coefficients.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     real freq(300)
!        the vector of frequencies at which the spectrum is computed.
!     integer i
!        an index variable
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list
!        if ierr == 0, no errors were detected
!        if ierr == 1, errors have been detected
!     integer iextnd
!        the indicator variable used to designate whether zero
!        (iextnd == 0) or the series mean (iextnd .ne. 0) is to be
!        used to extend the series.
!     integer itest
!        the number of the test set being performed
!     integer k(10)
!        the vector of the modified daniel filter lengths.
!     integer lab
!        the length of the vector ab.
!     integer ldstak
!        the length of the vector dstak in common cstak.
!     integer lfreq
!        the length of the vector freq.
!     integer lper
!        the length of the vector per.
!     integer lperi
!        the length of the vector peri.
!     integer lzfft
!        the length of the vectors yfft and zfft, respectively..
!     integer nf
!        the number of frequencies at which the spectrum is
!        to be computed.
!     integer nfft
!        the extended series length.
!     integer nk
!        the number of modified daniel filters to be applied.
!     integer nprt
!        the variable controling printed output, where
!        for the periodogram routines
!        if nprt <= -2, the output consists of a page plot of the
!                         periodogram on a log-linear scale,
!        if nprt == -1, the output consists of a page plot of the
!                         periodogram in decibels on a linear scale,
!        if nprt ==  0, the output is suppressed,
!        if nprt ==  1, the output consists of a vertical plot of the
!                         periodogram in decibels on a linear scale.
!        if nprt .ge.  2, the output consists of a vertical plot of the
!                         periodogram on a log-linear scale,
!        and for the integrated periodogram routines
!        if nprt ==  0, the output is suppressed,
!        if nprt .ne.  0, the output consists of a page plot of the
!                         integrated periodogram
!     integer ntemp
!        a temporary storage location
!     integer ny
!        the number of observations in the series y.
!     integer nz
!        the number of observations in the series z.
!     real per(300)
!        the series periodogram.
!     real perf(300)
!        the filtered (smoothed) periodogram.
!     real peri(300)
!        the series integrated periodogram.
!     real taperp
!        the percent of the series to be tapered.
!     real y(150)
!         the array containing the time series from jenkins and watts.
!     real yfft(400)
!        the vector of the observed time series to be analyzed using
!        the fft.
!     real z(275)
!        the array of the wolf sunspot numbers.
!     real zfft(600)
!        the vector of the tapered wolf sunspot numbers to be
!        analyzed using the fft.
!     real zt(275)
!        the array of the tapered sunspot numbers.
!

  integer &
     ldstak
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
     taperp
  integer &
     i,iextnd,itest,lab,lfreq,lper,lperi,lzfft,nf,nfft,nk, &
     nprt,ntemp,ny,nz
!
!  local arrays
  real &
     ab(600),freq(300),per(300),perf(300),peri(300),y(150), &
     yfft(400),z(275),zfft(600),zt(275)
  integer &
     k(10)
!
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr

  data   y(  1),  y(  2),  y(  3),  y(  4),  y(  5),  y(  6) &
      / -.88e0,  -.12e0,  -.89e0, -1.38e0,  -.07e0,  1.03e0/
  data   y(  7),  y(  8),  y(  9),  y( 10),  y( 11),  y( 12) &
      / 2.14e0,   .35e0, -1.10e0, -1.78e0, -2.76e0, -1.77e0/
  data   y( 13),  y( 14),  y( 15),  y( 16),  y( 17),  y( 18) &
      /  .98e0,  1.00e0,  -.70e0, -1.01e0, -1.30e0,  -.85e0/
  data   y( 19),  y( 20),  y( 21),  y( 22),  y( 23),  y( 24) &
      / -.46e0,  1.63e0,   .06e0,  -.17e0, -1.01e0, -1.04e0/
  data   y( 25),  y( 26),  y( 27),  y( 28),  y( 29),  y( 30) &
      / -.66e0, -1.12e0,  -.51e0,  -.71e0,  -.20e0,  -.13e0/
  data   y( 31),  y( 32),  y( 33),  y( 34),  y( 35),  y( 36) &
      /  .14e0,  1.59e0,  -.76e0, -1.08e0, -1.77e0, -1.20e0/
  data   y( 37),  y( 38),  y( 39),  y( 40),  y( 41),  y( 42) &
      /  .45e0,  -.07e0,  -.63e0,  -.35e0,  -.87e0,  -.62e0/
  data   y( 43),  y( 44),  y( 45),  y( 46),  y( 47),  y( 48) &
      /  .28e0,  1.90e0,  2.14e0,  1.05e0,   .31e0,  1.07e0/
  data   y( 49),  y( 50) &
      / 2.67e0,  2.44e0/
!
  data   z(  1),  z(  2),  z(  3),  z(  4),  z(  5),  z(  6) &
      /  5.0e0,  11.0e0,  16.0e0,  23.0e0,  36.0e0,  58.0e0/
  data   z(  7),  z(  8),  z(  9),  z( 10),  z( 11),  z( 12) &
      / 29.0e0,  20.0e0,  10.0e0,   8.0e0,   3.0e0,   0.0e0/
  data   z( 13),  z( 14),  z( 15),  z( 16),  z( 17),  z( 18) &
      /  0.0e0,   2.0e0,  11.0e0,  27.0e0,  47.0e0,  63.0e0/
  data   z( 19),  z( 20),  z( 21),  z( 22),  z( 23),  z( 24) &
      / 60.0e0,  39.0e0,  28.0e0,  26.0e0,  22.0e0,  11.0e0/
  data   z( 25),  z( 26),  z( 27),  z( 28),  z( 29),  z( 30) &
      / 21.0e0,  40.0e0,  78.0e0, 122.0e0, 103.0e0,  73.0e0/
  data   z( 31),  z( 32),  z( 33),  z( 34),  z( 35),  z( 36) &
      / 47.0e0,  35.0e0,  11.0e0,   5.0e0,  16.0e0,  34.0e0/
  data   z( 37),  z( 38),  z( 39),  z( 40),  z( 41),  z( 42) &
      / 70.0e0,  81.0e0, 111.0e0, 101.0e0,  73.0e0,  40.0e0/
  data   z( 43),  z( 44),  z( 45),  z( 46),  z( 47),  z( 48) &
      / 20.0e0,  16.0e0,   5.0e0,  11.0e0,  22.0e0,  40.0e0/
  data   z( 49),  z( 50),  z( 51),  z( 52),  z( 53),  z( 54) &
      / 60.0e0,  80.9e0,  83.4e0,  47.7e0,  47.8e0,  30.7e0/
  data   z( 55),  z( 56),  z( 57),  z( 58),  z( 59),  z( 60) &
      / 12.2e0,   9.6e0,  10.2e0,  32.4e0,  47.6e0,  54.0e0/
  data   z( 61),  z( 62),  z( 63),  z( 64),  z( 65),  z( 66) &
      / 62.9e0,  85.9e0,  61.2e0,  45.1e0,  36.4e0,  20.9e0/
  data   z( 67),  z( 68),  z( 69),  z( 70),  z( 71),  z( 72) &
      / 11.4e0,  37.8e0,  69.8e0, 106.1e0, 100.8e0,  81.6e0/
  data   z( 73),  z( 74),  z( 75),  z( 76),  z( 77),  z( 78) &
      / 66.5e0,  34.8e0,  30.6e0,   7.0e0,  19.8e0,  92.5e0/
  data   z( 79),  z( 80),  z( 81),  z( 82),  z( 83),  z( 84) &
      /154.4e0, 125.9e0,  84.8e0,  68.1e0,  38.5e0,  22.8e0/
  data   z( 85),  z( 86),  z( 87),  z( 88),  z( 89),  z( 90) &
      / 10.2e0,  24.1e0,  82.9e0, 132.0e0, 130.9e0, 118.1e0/
  data   z( 91),  z( 92),  z( 93),  z( 94),  z( 95),  z( 96) &
      / 89.9e0,  66.6e0,  60.0e0,  46.9e0,  41.0e0,  21.3e0/
  data   z( 97),  z( 98),  z( 99),  z(100),  z(101),  z(102) &
      / 16.0e0,   6.4e0,   4.1e0,   6.8e0,  14.5e0,  34.0e0/
  data   z(103),  z(104),  z(105),  z(106),  z(107),  z(108) &
      / 45.0e0,  43.1e0,  47.5e0,  42.2e0,  28.1e0,  10.1e0/
  data   z(109),  z(110),  z(111),  z(112),  z(113),  z(114) &
      /  8.1e0,   2.5e0,   0.0e0,   1.4e0,   5.0e0,  12.2e0/
  data   z(115),  z(116),  z(117),  z(118),  z(119),  z(120) &
      / 13.9e0,  35.4e0,  45.8e0,  41.1e0,  30.1e0,  23.9e0/
  data   z(121),  z(122),  z(123),  z(124),  z(125),  z(126) &
      / 15.6e0,   6.6e0,   4.0e0,   1.8e0,   8.5e0,  16.6e0/
  data   z(127),  z(128),  z(129),  z(130),  z(131),  z(132) &
      / 36.3e0,  49.6e0,  64.2e0,  67.0e0,  70.9e0,  47.8e0/
  data   z(133),  z(134),  z(135),  z(136),  z(137),  z(138) &
      / 27.5e0,   8.5e0,  13.2e0,  56.9e0, 121.5e0, 138.3e0/
  data   z(139),  z(140),  z(141),  z(142),  z(143),  z(144) &
      /103.2e0,  85.7e0,  64.6e0,  36.7e0,  24.2e0,  10.7e0/
  data   z(145),  z(146),  z(147),  z(148),  z(149),  z(150) &
      / 15.0e0,  40.1e0,  61.5e0,  98.5e0, 124.7e0,  96.3e0/
  data   z(151),  z(152),  z(153),  z(154),  z(155),  z(156) &
      / 66.6e0,  64.5e0,  54.1e0,  39.0e0,  20.6e0,   6.7e0/
  data   z(157),  z(158),  z(159),  z(160),  z(161),  z(162) &
      /  4.3e0,  22.7e0,  54.8e0,  93.8e0,  95.8e0,  77.2e0/
  data   z(163),  z(164),  z(165),  z(166),  z(167),  z(168) &
      / 59.1e0,  44.0e0,  47.0e0,  30.5e0,  16.3e0,   7.3e0/
  data   z(169),  z(170),  z(171),  z(172),  z(173),  z(174) &
      / 37.6e0,  74.0e0, 139.0e0, 111.2e0, 101.6e0,  66.2e0/
  data   z(175),  z(176),  z(177),  z(178),  z(179),  z(180) &
      / 44.7e0,  17.0e0,  11.3e0,  12.4e0,   3.4e0,   6.0e0/
  data   z(181),  z(182),  z(183),  z(184),  z(185),  z(186) &
      / 32.3e0,  54.3e0,  59.7e0,  63.7e0,  63.5e0,  52.2e0/
  data   z(187),  z(188),  z(189),  z(190),  z(191),  z(192) &
      / 25.4e0,  13.1e0,   6.8e0,   6.3e0,   7.1e0,  35.6e0/
  data   z(193),  z(194),  z(195),  z(196),  z(197),  z(198) &
      / 73.0e0,  85.1e0,  78.0e0,  64.0e0,  41.8e0,  26.2e0/
  data   z(199),  z(200),  z(201),  z(202),  z(203),  z(204) &
      / 26.7e0,  12.1e0,   9.5e0,   2.7e0,   5.0e0,  24.4e0/
  data   z(205),  z(206),  z(207),  z(208),  z(209),  z(210) &
      / 42.0e0,  63.5e0,  53.8e0,  62.0e0,  48.5e0,  43.9e0/
  data   z(211),  z(212),  z(213),  z(214),  z(215),  z(216) &
      / 18.6e0,   5.7e0,   3.6e0,   1.4e0,   9.6e0,  47.4e0/
  data   z(217),  z(218),  z(219),  z(220),  z(221),  z(222) &
      / 57.1e0, 103.9e0,  80.6e0,  63.6e0,  37.6e0,  26.1e0/
  data   z(223),  z(224),  z(225),  z(226),  z(227),  z(228) &
      / 14.2e0,   5.8e0,  16.7e0,  44.3e0,  63.9e0,  69.0e0/
  data   z(229),  z(230),  z(231),  z(232),  z(233),  z(234) &
      / 77.8e0,  64.9e0,  35.7e0,  21.2e0,  11.1e0,   5.7e0/
  data   z(235),  z(236),  z(237),  z(238),  z(239),  z(240) &
      /  8.7e0,  36.1e0,  79.7e0, 114.4e0, 109.6e0,  88.8e0/
  data   z(241),  z(242),  z(243),  z(244),  z(245),  z(246) &
      / 67.8e0,  47.5e0,  30.6e0,  16.3e0,   9.6e0,  33.2e0/
  data   z(247),  z(248),  z(249),  z(250),  z(251),  z(252) &
      / 92.6e0, 151.6e0, 136.3e0, 134.7e0,  83.9e0,  69.4e0/
  data   z(253),  z(254),  z(255),  z(256),  z(257),  z(258) &
      / 31.5e0,  13.9e0,   4.4e0,  38.0e0, 141.7e0, 190.2e0/
  data   z(259),  z(260),  z(261) &
      /184.8e0, 159.0e0, 112.3e0/

  itest = 1
!
!  make calls with valid data
!
  ny = 50
  nz = 261
  nfft = 514
  nprt = 2
  lzfft = 600
  lper = 514
  lperi = 514
  taperp = 0.10e0
  lfreq = 300
  iextnd = 0
  lab = 600
  nk = 3
  k(1) = 8
  k(2) = 8
  k(3) = 8
!
!  test of center
!
    5 write ( *, 1018)
  call center (z, nz, zt)
  write ( *, 1002) ierr
!
!  print returned variables from center
!
  if (ierr==0) write ( *, 1004) (zt(i), i = 1, nz)
!
!  test of taper
!
  write ( *, 1015)
  call taper (z, nz, taperp, zt)
  write ( *, 1002) ierr
!
!  print returned variables from taper
!
  if (ierr==0) write ( *, 1004) (zt(i), i = 1, nz)
!
!  test of pgm
!
  write ( *, 1013)
  call scopy (nz, zt, 1, zfft, 1)
  call pgm (zfft, nz, lzfft, ldstak)
  write ( *, 1002) ierr
!
!  test of fftlen
!
  write ( *, 1026)
  call fftlen(nfft-2, 2, ntemp)
  write ( *, 1002) ierr
!
!  print returned variables from fftlen
!
  if (ierr==0) write ( *, 1027) ntemp
!
!  test of pgms
!
  ntemp = nfft-1
  write ( *, 1025)
  call scopy (nz, zt, 1, zfft, 1)
  call pgms (zfft, nz, ntemp, lzfft, iextnd, nf, per, lper, freq, &
     lfreq, -2)
  write ( *, 1027) ntemp
  write ( *, 1002) ierr
!
!  print returned variables from pgms
!
  if (ierr==0) then
    write ( *, 1004) (freq(i), i = 1, nf)
    write ( *, 1004) (per(i), i = 1, nf)
  end if
!
!  test of mdflt
!
  write ( *, 1016)
  call mdflt (per, nf, nk, k, perf, ldstak)
  write ( *, 1002) ierr
  if (ierr==0) then
!
!  print returned variables from mdflt
!
    write ( *, 1004) (perf(i), i = 1, nf)
!
! display smoothed periodogram on a log plot
!
    write ( *, 1028)
    call ppl (perf, freq, nf, 1)
  end if
!
!  test of ipgmp
!
  write ( *, 1029)
  call ipgmp (per, freq, nf, nz, ldstak)
  write ( *, 1002) ierr
!
!  test of ipgmps
!
  write ( *, 1030)
  call ipgmps (per, freq, nf, nz, ldstak, peri, nprt)
  write ( *, 1002) ierr
  if (ierr==0) then
!
!  print returned variables from ipgmps
!
     write ( *, 1004) (freq(i), i = 1, nf)
     write ( *, 1004) (peri(i), i = 1, nf)
  end if
!
!  test of ipgm
!
  write ( *, 1017)
  call scopy (nz, zt, 1, zfft, 1)
  call ipgm (zfft, nz, lzfft, ldstak)
  write ( *, 1002) ierr
!
!  test of ipgms
!
  write ( *, 1014)
  call scopy (nz, zt, 1, zfft, 1)
  call ipgms (zfft, nz, lzfft, ldstak, nf, peri, lperi, freq, lfreq, &
     nprt)
  write ( *, 1002) ierr
  if (ierr==0) then
!
!  print returned variables from ipgms
!
     write ( *, 1004) (freq(i), i = 1, nf)
     write ( *, 1004) (peri(i), i = 1, nf)
  end if
!
!  test of fftr (centered data - o percent taper)
!
  taperp = -1.0e0
  write ( *, 1031)
  if (ny.ge.1) call taper (y, ny, taperp, yfft)
  call fftr (yfft, ny, nfft, iextnd, nf, ab, lab)
  write ( *, 1002) ierr
!
!  print returned variables from fftr
!
  if (ierr==0) write ( *, 1004) (ab(i), i = 1, nf)
!
  go to (10, 20, 30, 40) itest
!
!  check minimum problem size
!
   10 itest = itest + 1
  nz = 17
  ny = 17
  go to 5
!
!  check various options
!
   20 itest = itest + 1
!
!  test of mdflt (elements of k not even)
!
  k(1) = 7
  write ( *, 1016)
  call mdflt (per, nf, nk, k, perf, ldstak)
  write ( *, 1002) ierr
!
!  print returned variables from mdflt
!
  write ( *, 1004) (perf(i), i = 1, nf)
!
!  test of pgms (uncentered data)
!
  iextnd = 1
  nprt = 1
  write ( *, 1025)
  call scopy (nz, z, 1, zfft, 1)
  call pgms (zfft, nz, nfft, lzfft, iextnd, nf, per, lper, freq, &
     lfreq, nprt)
  write ( *, 1002) ierr
!
!  print returned variables from pgms
!
  write ( *, 1004) (freq(i), i = 1, nf)
  write ( *, 1004) (per(i), i = 1, nf)

  nprt = 2
  write ( *, 1025)
  call scopy (nz, z, 1, zfft, 1)
  call pgms (zfft, nz, nfft, lzfft, iextnd, nf, per, lper, freq, &
     lfreq, nprt)
  write ( *, 1002) ierr
!
!  print returned variables from pgms
!
  write ( *, 1004) (freq(i), i = 1, nf)
  write ( *, 1004) (per(i), i = 1, nf)
!
!  test of fftr (centered data - 100 percent taper)
!
  taperp = 1.1e0
  write ( *, 1031)
  call taper (y, ny, taperp, yfft)
  call fftr (yfft, ny, nfft, iextnd, nf, ab, lab)
  write ( *, 1002) ierr
!
!  print returned variables from fftr
!
  write ( *, 1004) (ab(i), i = 1, nf)
!
!  test of fftr (centered data - o percent taper)
!
  taperp = -1.0e0
  write ( *, 1031)
  call taper (y, ny, taperp, yfft)
  call fftr (yfft, ny, nfft-1, iextnd, nf, ab, nfft-2)
  write ( *, 1002) ierr
!
!  print returned variables from fftr
!
  write ( *, 1004) (ab(i), i = 1, nf)
!
!  perform error checking
!
  ny = -1
  nz = -1
  nfft = 0
  nprt = 2
  lzfft = 0
  lper = 0
  lperi = 0
  taperp = 0.10e0
  lfreq = 0
  iextnd = 0
  lab = 0
  nk = 0
  k(1) = 0
  k(2) = 0
  k(3) = 0
  go to 5
!
!     perform more error checking
!
   30 itest = itest + 1
  ny = 50
  nz = 261
  nprt = 2
  lzfft = 0
  lper = 0
  lperi = 0
  taperp = 0.10e0
  lfreq = 0
  iextnd = 0
  lab = 0
  nk = 3
  k(1) = 0
  k(2) = 0
  k(3) = 0
  go to 5

   40 return

 1002 format ('ierr is ', i5)
 1004 format (3(1x, e16.8))
 1013 format ('test of pgm')
 1014 format ('test of ipgms')
 1015 format ('test of taper')
 1016 format ('test of mdflt')
 1017 format ('test of ipgm')
 1018 format ('test of center')
 1025 format ('test of pgms')
 1026 format ('test of fftlen')
 1027 format (/'nfft is ', i6)
 1028 format ('display of periodogram smoothed with modified', &
     ' daniel filter')
 1029 format ('test of ipgmp')
 1030 format ('test of ipgmps')
 1031 format ('test of fftr')
end
