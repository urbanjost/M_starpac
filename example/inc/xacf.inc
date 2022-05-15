subroutine xacf ( lds ,typ)
use M_starpac
implicit none
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XACF tests the time series correlation subroutines.
!
!  Discussion:
!
!    Series y is listed as series x1 on page 362 in jenkins and watts.
!
!    Series yd is listed as series g on page 531 of box and jenkins.
!
!  Modified:
!
!    01 December 2006
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
!    George Box, Gwilym Jenkins,
!    Time Series Analysis: Forecasting and Control,
!    Holden-Day, 1970,
!    QA280.B67.
!
!    Gwilym Jenkins, Donald Watts,
!    Spectral Analysis and its Applications,
!    Holden-Day 1968.
!
!  Parameters:
!
!     real acov(21)
!        the autocovariance vector.
!     real amiss
!        the missing value code for the returned acvf estimates
!        (vector acov).
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer i
!        an indexing variable.
!     integer iar
!        the order of the autoregressive process chosen.
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list
!        if ierr == 0, no errors were detected
!     integer iod(2)
!        the order of each of the difference factors.
!     integer itest
!        the number of the test being run
!     integer lacov
!        the length of the acvf related vectors.
!     integer lagmax
!        the maximum lag value requested.
!     integer ldstak
!        the length of the array dstak.
!     integer lyfft
!        the length of the arrays used when the computations are
!        performed by the fft.
!     integer n
!        the integer number of observations in each series
!     integer nd(2)
!        the array containing the number of times the difference
!        factors are to be applied.
!     integer nfac
!        the number of difference factors.
!     integer nlppa(21)
!        the array containing the number of lagged product pairs
!        used to compute each acvf estimate.
!     integer nprt
!        the indicator variable used to specify whether or not
!        printed output is to be given, where if the value of
!        nprt is zero, no output is made.
!     integer nyd
!        the number of observations in the series to be differenced.
!     real phi(21)
!        the array of autoregressive coefficients for the selected
!        order.
!     real y(100), yd(150)
!        the vector containing the observed time series
!     real yfft(150)
!        the vectors used for storing the series for the routines
!        using the fft.
!     real ymiss
!        the missing value codes for series y and ym.
!

  real acov(21)
  real amiss
  double precision dstak(12)
  integer i
  integer iar
  integer ierr
  integer iod(2)
  integer itest
  integer lacov
  integer lagmax
  integer lds
  integer ldstak
  integer lyfft
  integer n
  integer nd(2)
  integer nfac
  integer nlppa(21)
  integer nprt
  integer nyd
  real phi(21)
  real y(100)
  real yd(150)
  real yfft(150)
  real ymiss


  common /cstak/ dstak
  common /errchk/ ierr

  data    y(  1),   y(  2),   y(  3),   y(  4),   y(  5),   y(  6) &
      / -2.07e0, -1.15e0,  0.69e0, -0.46e0, -1.49e0, -0.70e0/
  data    y(  7),   y(  8),   y(  9),   y( 10),   y( 11),   y( 12) &
      / -1.07e0, -0.69e0, -0.68e0,  1.27e0, -1.05e0, -0.05e0/
  data    y( 13),   y( 14),   y( 15),   y( 16),   y( 17),   y( 18) &
      / -0.84e0, -0.62e0, -0.49e0, -1.29e0, -0.49e0, -1.06e0/
  data    y( 19),   y( 20),   y( 21),   y( 22),   y( 23),   y( 24) &
      / -0.38e0, -0.52e0, -0.13e0,  1.30e0, -1.51e0, -0.43e0/
  data    y( 25),   y( 26),   y( 27),   y( 28),   y( 29),   y( 30) &
      / -1.33e0, -0.78e0,  0.31e0, -0.95e0, -0.90e0, -0.30e0/
  data    y( 31),   y( 32),   y( 33),   y( 34),   y( 35),   y( 36) &
      / -1.02e0, -0.53e0,  0.15e0,  1.40e0,  1.22e0,  0.59e0/
  data    y( 37),   y( 38),   y( 39),   y( 40),   y( 41),   y( 42) &
      /  0.70e0,  1.70e0,  2.78e0,  1.98e0,  1.39e0,  1.85e0/
  data    y( 43),   y( 44),   y( 45),   y( 46),   y( 47),   y( 48) &
      /  2.60e0,  0.51e0,  2.77e0,  1.16e0,  1.07e0, -0.48e0/
  data    y( 49),   y( 50),   y( 51),   y( 52),   y( 53),   y( 54) &
      / -0.52e0,  0.37e0,  0.00e0, -1.99e0, -1.75e0,  0.70e0/
  data    y( 55),   y( 56),   y( 57),   y( 58),   y( 59),   y( 60) &
      /  0.73e0,  1.16e0,  0.06e0, -0.02e0,  1.10e0, -0.35e0/
  data    y( 61),   y( 62),   y( 63),   y( 64),   y( 65),   y( 66) &
      / -1.67e0, -1.57e0,  1.16e0,  1.84e0,  3.35e0,  0.40e0/
  data    y( 67),   y( 68),   y( 69),   y( 70),   y( 71),   y( 72) &
      /  0.45e0,  1.30e0,  0.93e0,  1.17e0, -1.74e0, -1.28e0/
  data    y( 73),   y( 74),   y( 75),   y( 76),   y( 77),   y( 78) &
      / -0.07e0,  1.50e0,  0.53e0,  0.20e0, -0.42e0,  1.18e0/
  data    y( 79),   y( 80),   y( 81),   y( 82),   y( 83),   y( 84) &
      /  0.82e0,  1.50e0,  2.92e0,  1.18e0,  1.23e0,  3.16e0/
  data    y( 85),   y( 86),   y( 87),   y( 88),   y( 89),   y( 90) &
      /  0.79e0,  0.68e0,  1.14e0,  1.02e0,  1.02e0, -0.71e0/
  data    y( 91),   y( 92),   y( 93),   y( 94),   y( 95),   y( 96) &
      / -0.17e0, -1.50e0, -0.26e0, -0.38e0,  0.93e0, -0.33e0/
  data    y( 97),   y( 98),   y( 99),   y(100) &
      / -1.12e0, -2.95e0, -2.09e0, -1.11e0                    /
!
  data   yd(  1),  yd(  2),  yd(  3),  yd(  4),  yd(  5),  yd(  6) &
      /  112.0e0, 118.0e0, 132.0e0, 129.0e0, 121.0e0, 135.0e0/
  data   yd(  7),  yd(  8),  yd(  9),  yd( 10),  yd( 11),  yd( 12) &
      /  148.0e0, 148.0e0, 136.0e0, 119.0e0, 104.0e0, 118.0e0/
  data   yd( 13),  yd( 14),  yd( 15),  yd( 16),  yd( 17),  yd( 18) &
      /  115.0e0, 126.0e0, 141.0e0, 135.0e0, 125.0e0, 149.0e0/
  data   yd( 19),  yd( 20),  yd( 21),  yd( 22),  yd( 23),  yd( 24) &
      /  170.0e0, 170.0e0, 158.0e0, 133.0e0, 114.0e0, 140.0e0/
  data   yd( 25),  yd( 26),  yd( 27),  yd( 28),  yd( 29),  yd( 30) &
      /  145.0e0, 150.0e0, 178.0e0, 163.0e0, 172.0e0, 178.0e0/
  data   yd( 31),  yd( 32),  yd( 33),  yd( 34),  yd( 35),  yd( 36) &
      /  199.0e0, 199.0e0, 184.0e0, 162.0e0, 146.0e0, 166.0e0/
  data   yd( 37),  yd( 38),  yd( 39),  yd( 40),  yd( 41),  yd( 42) &
      /  171.0e0, 180.0e0, 193.0e0, 181.0e0, 183.0e0, 218.0e0/
  data   yd( 43),  yd( 44),  yd( 45),  yd( 46),  yd( 47),  yd( 48) &
      /  230.0e0, 242.0e0, 209.0e0, 191.0e0, 172.0e0, 194.0e0/
  data   yd( 49),  yd( 50),  yd( 51),  yd( 52),  yd( 53),  yd( 54) &
      /  196.0e0, 196.0e0, 236.0e0, 235.0e0, 229.0e0, 243.0e0/
  data   yd( 55),  yd( 56),  yd( 57),  yd( 58),  yd( 59),  yd( 60) &
      /  264.0e0, 272.0e0, 237.0e0, 211.0e0, 180.0e0, 201.0e0/
  data   yd( 61),  yd( 62),  yd( 63),  yd( 64),  yd( 65),  yd( 66) &
      /  204.0e0, 188.0e0, 235.0e0, 227.0e0, 234.0e0, 264.0e0/
  data   yd( 67),  yd( 68),  yd( 69),  yd( 70),  yd( 71),  yd( 72) &
      /  302.0e0, 293.0e0, 259.0e0, 229.0e0, 203.0e0, 229.0e0/
  data   yd( 73),  yd( 74),  yd( 75),  yd( 76),  yd( 77),  yd( 78) &
      /  242.0e0, 233.0e0, 267.0e0, 269.0e0, 270.0e0, 315.0e0/
  data   yd( 79),  yd( 80),  yd( 81),  yd( 82),  yd( 83),  yd( 84) &
      /  364.0e0, 347.0e0, 312.0e0, 274.0e0, 237.0e0, 278.0e0/
  data   yd( 85),  yd( 86),  yd( 87),  yd( 88),  yd( 89),  yd( 90) &
      /  284.0e0, 277.0e0, 317.0e0, 313.0e0, 318.0e0, 374.0e0/
  data   yd( 91),  yd( 92),  yd( 93),  yd( 94),  yd( 95),  yd( 96) &
      /  413.0e0, 405.0e0, 355.0e0, 306.0e0, 271.0e0, 306.0e0/
  data   yd( 97),  yd( 98),  yd( 99),  yd(100),  yd(101),  yd(102) &
      /  315.0e0, 301.0e0, 356.0e0, 348.0e0, 355.0e0, 422.0e0/
  data   yd(103),  yd(104),  yd(105),  yd(106),  yd(107),  yd(108) &
      /  465.0e0, 467.0e0, 404.0e0, 347.0e0, 305.0e0, 336.0e0/
  data   yd(109),  yd(110),  yd(111),  yd(112),  yd(113),  yd(114) &
      /  340.0e0, 318.0e0, 362.0e0, 348.0e0, 363.0e0, 435.0e0/
  data   yd(115),  yd(116),  yd(117),  yd(118),  yd(119),  yd(120) &
      /  491.0e0, 505.0e0, 404.0e0, 359.0e0, 310.0e0, 337.0e0/
  data   yd(121),  yd(122),  yd(123),  yd(124),  yd(125),  yd(126) &
      /  360.0e0, 342.0e0, 406.0e0, 396.0e0, 420.0e0, 472.0e0/
  data   yd(127),  yd(128),  yd(129),  yd(130),  yd(131),  yd(132) &
      /  548.0e0, 559.0e0, 463.0e0, 407.0e0, 362.0e0, 405.0e0/
  data   yd(133),  yd(134),  yd(135),  yd(136),  yd(137),  yd(138) &
      /  417.0e0, 391.0e0, 419.0e0, 461.0e0, 472.0e0, 535.0e0/
  data   yd(139),  yd(140),  yd(141),  yd(142),  yd(143),  yd(144) &
      /  622.0e0, 606.0e0, 508.0e0, 461.0e0, 390.0e0, 432.0e0/

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'XACF'
  write ( *, '(a)' ) '  Test the time series correlation routines.'
  write ( *, '(a)' ) ' '

  itest = 1
  ldstak = lds

  n = 100
  lagmax = 20
  nprt = 1
  lyfft = 150
  lacov = 21
  nyd = 144
  nfac = 2
  nd(1) = 1
  nd(2) = 1
  iod(1) = 12
  iod(2) = 1
  ymiss = 1.16e0
!
!  Test of acf
!
5 continue

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of ACF:'
  write ( *, '(a)' ) ' '

  call acf ( y, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  IERR = ', ierr
!
!  Test of acfs
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of ACFS:'
  write ( *, '(a)' ) ' '

  call acfs ( y, n, lagmax, lacov, acov, iar, phi, nprt, ldstak )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  IERR = ', ierr
!
!  print storage from acfs
!
  if (ierr==0) then
    write ( *,'(9f10.5)') (acov(i),i=1,lagmax+1)
    write ( *,'(9f10.5)') (phi(i),i=1,iar)
  end if
!
!  Test of acfd
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of ACFD:'
  write ( *, '(a)' ) ' '

  call acfd(yd, nyd, lagmax, nfac, nd, iod, ldstak)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  IERR = ', ierr
!
!  Test of acfm
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of ACFM:'
  write ( *, '(a)' ) ' '

  call acfm ( y, ymiss, n )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  IERR = ', ierr
!
!  Test of acfms
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of ACFMS:'
  write ( *, '(a)' ) ' '

  call acfms(y, ymiss, n, lagmax, lacov, acov, amiss, nlppa, nprt, &
     ldstak)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  IERR = ', ierr
!
!  print storage from acfms
!
  if (ierr==0) then
    write ( *,'(9f10.5)') (acov(i),i=1,lagmax+1)
    write ( *,'(9i10)') (nlppa(i),i=1,lagmax+1)
  end if
!
!  copy data into yfft for acff
!
  yfft(1:n) = y(1:n)
!
!  Test of acff
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of ACFF:'
  write ( *, '(a)' ) ' '

  call acff ( yfft, n, lyfft, ldstak )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  IERR = ', ierr
!
!  copy data into yfft for acffs
!
  yfft(1:n) = y(1:n)
!
!  Test of acffs
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of ACFFS:'
  write ( *, '(a)' ) ' '

  call acffs ( yfft, n, lyfft, ldstak, lagmax, lacov, acov, iar, phi, &
     nprt )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  IERR = ', ierr
!
!  print storage from acffs
!
  if ( ierr == 0 ) then
    write ( *,'(9f10.5)') (acov(i),i=1,lagmax+1)
    write ( *,'(9f10.5)') (phi(i),i=1,iar)
  end if
!
!  Test minimum problem size
!
  if ( itest == 1 ) then

    itest = itest + 1
    n = 13
    lagmax = 1
    nfac = 1
    nd(1) = 1
    iod(1) = 1
!
!  Check error handling
!
  else if ( itest == 2 ) then

    itest = itest + 1
    n = 0
    lagmax = 20
    lyfft = 0
    lacov = 0
    nyd = 0
    nfac = 1
    nd(1) = 0
    iod(1) = 0
    go to 5
!
!  Check error handling
!
  else if ( itest == 3 ) then

    itest = itest + 1
    n = 100
    lagmax = 0
    lyfft = 0
    lacov = 0
    nyd = 144
    nfac = 0
    ldstak = 0
    go to 5

  end if

  return
end
