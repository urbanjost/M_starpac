subroutine xccf ( lds ,typ)
use M_starpac, only : ccf,ccff,ccffs,ccfm,ccfms,ccfs,ccfxp,iprint,scopy,setra, setrv
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XCCF tests the time series correlation subroutines.
!
!  Discussion:
!
!    series y1 and y2 are listed as series x1 and x2 on page of 361 of
!    jenkins and watts.  ccf for series y1 and y2 are plotted on page 3
!    and listed on page 420.
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
!  Parameters:
!
!     real ccov(30,5,5)
!        the cross covariance array.
!     real cmiss
!        the missing value code for the returned ccvf estimates
!        (vector ccov).
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer iccov
!        the first dimension of the array ccov.
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list
!        if ierr == 0, no errors were detected
!     integer inlppc
!        the first dimension of the array nlppc.
!     integer itest
!        the number of the test being performed
!     integer iym, iymfft
!        the first dimension of the arrays ym and ymfft, respectively.
!     integer jccov, jnlppc
!        the second dimensions of the arrays ccov and nlppc,
!        respectively.
!     integer lagmax
!        the maximum lag value requested.
!     integer lds, ldstak
!        the length of the array dstak.
!     integer lyfft
!        the length of the arrays used when the computations are
!        performed by the fft.
!     integer m
!        the number of series in the multivariate time series ym.
!     integer n
!        the integer number of observations in each series
!     integer nlag
!        the number of lags at which the acvf was computed.
!     integer nlppc(30,5,5)
!        the array containing the number of lagged product pairs
!        used to compute each acvf estimate.
!     integer nprt
!        the indicator variable used to specify whether or not
!        printed output is to be given, where if the value of
!        nprt is zero, no output is made.
!     integer nyd
!        the number of observations in the series to be differenced.
!     real yfft1(150), yfft2(150)
!        the vectors used for storing the series for the routines
!        using the fft.
!     real ym(150,5), ymfft(150,5)
!        the arrays used for multivariate time series.
!     real ymiss0, ymmiss(5)
!        the missing value codes for series y and ym.
!     real y1(100), y1miss
!        the first series, and its missing value code.
!     real y2(100), y2miss
!        the second series, and its missing value code.
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
     cmiss,y1miss,y2miss,ymiss0
  integer &
     iccov,inlppc,itest,iym,iymfft,jccov,jnlppc,lagmax, &
     ldstak,lyfft,m,n,nlag,nprt,nyd
!
!  local arrays
  real &
     ccov(30,5,5),y1(100),y2(100),yfft1(150),yfft2(150),ym(150,5), &
     ymfft(150,5),ymmiss(5)
  integer &
     nlppc(30,5,5)
!
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr

  data   y1(  1),  y1(  2),  y1(  3),  y1(  4),  y1(  5),  y1(  6) &
      /-0.88e0, -0.16e0, -1.87e0, -1.12e0,  1.38e0,  2.13e0/
  data   y1(  7),  y1(  8),  y1(  9),  y1( 10),  y1( 11),  y1( 12) &
      / 2.76e0,  0.56e0, -0.69e0, -1.79e0, -3.82e0, -2.38e0/
  data   y1( 13),  y1( 14),  y1( 15),  y1( 16),  y1( 17),  y1( 18) &
      / 1.00e0,  0.70e0, -0.15e0,  0.98e0,  0.11e0, -0.35e0/
  data   y1( 19),  y1( 20),  y1( 21),  y1( 22),  y1( 23),  y1( 24) &
      /-0.73e0,  0.89e0, -1.63e0, -0.44e0, -1.37e0, -1.71e0/
  data   y1( 25),  y1( 26),  y1( 27),  y1( 28),  y1( 29),  y1( 30) &
      /-1.22e0, -2.00e0, -0.22e0,  0.38e0,  1.31e0,  0.71e0/
  data   y1( 31),  y1( 32),  y1( 33),  y1( 34),  y1( 35),  y1( 36) &
      / 0.32e0,  0.48e0, -1.88e0, -0.94e0, -1.54e0, -0.13e0/
  data   y1( 37),  y1( 38),  y1( 39),  y1( 40),  y1( 41),  y1( 42) &
      / 1.02e0,  0.02e0, -0.77e0,  0.11e0, -0.60e0, -0.52e0/
  data   y1( 43),  y1( 44),  y1( 45),  y1( 46),  y1( 47),  y1( 48) &
      /-0.09e0,  1.23e0,  1.46e0,  0.61e0,  0.42e0,  2.16e0/
  data   y1( 49),  y1( 50),  y1( 51),  y1( 52),  y1( 53),  y1( 54) &
      / 3.18e0,  2.10e0,  0.37e0, -0.24e0,  0.57e0, -0.53e0/
  data   y1( 55),  y1( 56),  y1( 57),  y1( 58),  y1( 59),  y1( 60) &
      / 2.44e0,  1.02e0, -0.53e0, -2.49e0, -2.12e0, -1.04e0/
  data   y1( 61),  y1( 62),  y1( 63),  y1( 64),  y1( 65),  y1( 66) &
      /-0.12e0, -1.88e0, -1.50e0,  1.54e0,  3.33e0,  3.08e0/
  data   y1( 67),  y1( 68),  y1( 69),  y1( 70),  y1( 71),  y1( 72) &
      / 1.71e0,  0.79e0,  1.55e0,  0.89e0, -0.89e0, -1.18e0/
  data   y1( 73),  y1( 74),  y1( 75),  y1( 76),  y1( 77),  y1( 78) &
      / 0.89e0,  1.71e0,  3.05e0,  0.15e0, -1.04e0,  0.12e0/
  data   y1( 79),  y1( 80),  y1( 81),  y1( 82),  y1( 83),  y1( 84) &
      / 0.08e0,  0.11e0, -2.62e0, -1.28e0,  1.07e0,  3.20e0/
  data   y1( 85),  y1( 86),  y1( 87),  y1( 88),  y1( 89),  y1( 90) &
      / 1.92e0,  0.53e0, -1.08e0,  0.49e0, -0.58e0,  0.17e0/
  data   y1( 91),  y1( 92),  y1( 93),  y1( 94),  y1( 95),  y1( 96) &
      / 1.15e0, -0.97e0, -1.63e0,  1.14e0, -0.67e0, -0.88e0/
  data   y1( 97),  y1( 98),  y1( 99),  y1(100) &
      /-0.07e0,  0.24e0,  0.55e0, -2.16e0/
  data   y2(  1),  y2(  2),  y2(  3),  y2(  4),  y2(  5),  y2(  6) &
      / 0.79e0,  1.12e0, -1.10e0, -2.39e0, -1.75e0, -0.82e0/
  data   y2(  7),  y2(  8),  y2(  9),  y2( 10),  y2( 11),  y2( 12) &
      /-0.36e0,  1.27e0,  1.75e0,  2.44e0,  0.36e0, -2.10e0/
  data   y2( 13),  y2( 14),  y2( 15),  y2( 16),  y2( 17),  y2( 18) &
      /-1.93e0, -1.30e0, -1.75e0, -0.34e0,  0.74e0,  0.49e0/
  data   y2( 19),  y2( 20),  y2( 21),  y2( 22),  y2( 23),  y2( 24) &
      / 0.70e0,  0.71e0,  0.09e0,  0.59e0,  1.54e0,  0.14e0/
  data   y2( 25),  y2( 26),  y2( 27),  y2( 28),  y2( 29),  y2( 30) &
      / 0.55e0, -1.40e0, -2.55e0, -1.66e0, -0.43e0,  0.58e0/
  data   y2( 31),  y2( 32),  y2( 33),  y2( 34),  y2( 35),  y2( 36) &
      / 2.18e0, -0.24e0,  0.58e0, -0.18e0, -1.55e0, -0.64e0/
  data   y2( 37),  y2( 38),  y2( 39),  y2( 40),  y2( 41),  y2( 42) &
      /-1.09e0,  0.90e0, -0.66e0, -0.35e0,  0.48e0,  0.50e0/
  data   y2( 43),  y2( 44),  y2( 45),  y2( 46),  y2( 47),  y2( 48) &
      / 0.05e0, -0.68e0,  0.24e0,  0.58e0, -1.26e0, -0.25e0/
  data   y2( 49),  y2( 50),  y2( 51),  y2( 52),  y2( 53),  y2( 54) &
      / 0.25e0,  2.18e0,  2.96e0,  1.56e0, -0.36e0, -0.59e0/
  data   y2( 55),  y2( 56),  y2( 57),  y2( 58),  y2( 59),  y2( 60) &
      /-0.12e0,  3.03e0,  2.11e0,  0.78e0,  0.89e0, -1.45e0/
  data   y2( 61),  y2( 62),  y2( 63),  y2( 64),  y2( 65),  y2( 66) &
      /-0.36e0, -0.37e0, -1.39e0, -4.19e0, -0.73e0, -0.98e0/
  data   y2( 67),  y2( 68),  y2( 69),  y2( 70),  y2( 71),  y2( 72) &
      / 0.36e0,  0.06e0, -1.94e0, -0.08e0,  0.17e0,  1.00e0/
  data   y2( 73),  y2( 74),  y2( 75),  y2( 76),  y2( 77),  y2( 78) &
      /-0.05e0,  0.43e0,  0.15e0,  2.69e0,  0.57e0,  0.29e0/
  data   y2( 79),  y2( 80),  y2( 81),  y2( 82),  y2( 83),  y2( 84) &
      / 1.10e0,  0.48e0, -1.06e0, -2.28e0, -2.03e0, -0.75e0/
  data   y2( 85),  y2( 86),  y2( 87),  y2( 88),  y2( 89),  y2( 90) &
      / 1.00e0,  1.71e0,  0.58e0,  1.97e0,  0.99e0,  1.94e0/
  data   y2( 91),  y2( 92),  y2( 93),  y2( 94),  y2( 95),  y2( 96) &
      / 2.18e0,  3.14e0,  0.60e0,  0.51e0,  1.35e0,  0.56e0/
  data   y2( 97),  y2( 98),  y2( 99),  y2(100) &
      / 0.11e0,  0.00e0,  2.34e0,  1.88e0/

  itest = 1
  ldstak = lds

  n = 100
  lagmax = 20
  nlag = 30
  nprt = 1
  lyfft = 150
  iccov = 30
  jccov = 5
  iym = 150
  m = 4
  iymfft = 150
  inlppc = 30
  jnlppc = 5
  nyd = 144
  ymiss0 = 1.16e0
  y1miss = 0.89e0
  y2miss = 0.89e0
!
!  copy data into ym for ccfs and ccfms
!
  call scopy(n, y1, 1, ym(1,1), 1)
  call scopy(n, y2, 1, ym(1,2), 1)
  call scopy(n, y1, 1, ym(1,3), 1)
  call scopy(n, y2, 1, ym(1,4), 1)
  call setrv(ymmiss, 4, ymiss0)
!
!  test of ccf
!
  write ( *,1060)
  call ccf(y1, y2, n)
!
!  print returned results
!
  call ccfxp (.false., lagmax, m, ccov, iccov, jccov, .false., &
     nlppc,  inlppc, jnlppc, cmiss)
!
!  test of ccfs
!
  write ( *,1080)
  call ccfs(ym, n, m, iym, lagmax, ccov, iccov, jccov, nprt, &
     ldstak)
!
!  print returned results
!
  call ccfxp (.true., lagmax, m, ccov, iccov, jccov, .false., &
     nlppc,  inlppc, jnlppc, cmiss)
!
!  test of ccfm without missing values
!
  write ( *,1070)
  write ( *, 1050)
  call ccfm(y1, ymiss0, y2, ymiss0, n)
!
!  print returned results
!
  call ccfxp ( .false., lagmax, m, ccov, iccov, jccov, .true., &
     nlppc,  inlppc, jnlppc, cmiss )
!
!  test of ccfms without missing values
!
  write ( *,1140)
  write ( *, 1050)
  call ccfms ( ym, ymmiss, n, m, iym, lagmax, ccov, cmiss, &
     iccov, jccov, nlppc, inlppc, jnlppc, nprt, ldstak )
!
!  print returned results
!
  call ccfxp ( .true., lagmax, m, ccov, iccov, jccov, .true., &
     nlppc,  inlppc, jnlppc, cmiss )
!
!  copy data into yfft1, yfft2 and ymfft for ccff and ccffs
!
  call scopy(n, y1, 1, yfft1, 1)
  call scopy(n, y2, 1, yfft2, 1)
  call scopy(n, y1, 1, ymfft(1,1), 1)
  call scopy(n, y2, 1, ymfft(1,2), 1)
  call scopy(n, y1, 1, ymfft(1,3), 1)
  call scopy(n, y2, 1, ymfft(1,4), 1)
!
!  test of ccff
!
  write ( *,1100)
  call ccff(yfft1, yfft2, n, lyfft, ldstak)
!
!  print returned results
!
  call ccfxp (.false., lagmax, m, ccov, iccov, jccov, .false., &
     nlppc,  inlppc, jnlppc, cmiss)
!
!  test of ccffs
!
  write ( *,1150)
  call ccffs(ymfft, n, m, iymfft, lagmax, ccov, &
     iccov, jccov, nprt, ldstak)
!
!  print returned results
!
  call ccfxp (.true., lagmax, m, ccov, iccov, jccov, .false., &
     nlppc,  inlppc, jnlppc, cmiss)
!
!  reset ymmiss
!
  ymmiss(1) = y1miss
  ymmiss(2) = y2miss
  ymmiss(3) = y1miss
  ymmiss(4) = y2miss
!
!  test of ccfm with missing values
!
  write ( *,1070)
  write ( *, 1040)
  call ccfm(y1, y1miss, y2, y2miss, n)
!
!  print returned results
!
  call ccfxp (.false., lagmax, m, ccov, iccov, jccov, .true., &
     nlppc,  inlppc, jnlppc, cmiss)
!
!  test of ccfms with missing values
!
  write ( *,1140)
  write ( *, 1040)
  call ccfms(ym, ymmiss, n, m, iym, lagmax, ccov, cmiss, &
     iccov, jccov, nlppc, inlppc, jnlppc, nprt, ldstak)
!
!     print returned results
!
  call ccfxp (.true., lagmax, m, ccov, iccov, jccov, .true., &
     nlppc,  inlppc, jnlppc, cmiss)
!
!  test print control
!
  nprt = 0
!
!  test of ccfs
!
  write ( *,1080)
  write ( *, 1020)
  call ccfs(ym, n, m, lagmax, iym, ccov, iccov, jccov, nprt, &
     ldstak)
!
!  print returned results
!
  call ccfxp (.true., lagmax, m, ccov, iccov, jccov, .false., &
     nlppc,  inlppc, jnlppc, cmiss)
!
!  test of ccfms with missing values
!
  write ( *,1140)
  write ( *, 1040)
  write ( *, 1020)
  call ccfms(ym, ymmiss, n, m, iym, lagmax, ccov, cmiss, &
     iccov, jccov, nlppc, inlppc, jnlppc, nprt, ldstak)
!
!  print returned results
!
  call ccfxp (.true., lagmax, m, ccov, iccov, jccov, .true., &
     nlppc,  inlppc, jnlppc, cmiss)
!
!  copy data into ymfft for ccffs
!
  call scopy(n, y1, 1, ymfft(1,1), 1)
  call scopy(n, y2, 1, ymfft(1,2), 1)
  call scopy(n, y1, 1, ymfft(1,3), 1)
  call scopy(n, y2, 1, ymfft(1,4), 1)
!
!  test of ccffs
!
  write ( *,1150)
  write ( *, 1020)
  call ccffs(ymfft, n, m, iymfft, lagmax, ccov, &
     iccov, jccov, nprt, ldstak)
!
!  print returned results
!
  call ccfxp (.true., lagmax, m, ccov, iccov, jccov, .false., &
     nlppc,  inlppc, jnlppc, cmiss)
!
!  test lead/lag message
!
  nprt = 1

  call setra(ymfft, iymfft, m, n, 0.0e0)
  ymfft(5,1) = 1.0e0
  ymfft(15,2) = 1.0e0
  ymfft(5,3) = ymfft(5,1)
  ymfft(15,4) = ymfft(15,2)
!
!  test of ccffs
!
  write ( *,1150)
  write ( *, 1020)
  call ccffs(ymfft, n, m, iymfft, lagmax, ccov, &
     iccov, jccov, nprt, ldstak)
!
!  print returned results
!
  call ccfxp (.true., lagmax, m, ccov, iccov, jccov, .false., &
     nlppc,  inlppc, jnlppc, cmiss)
!
  go to (100, 200, 300, 400), itest
!
!     test minimum problem size
!
  100 itest = itest + 1
  n = 3
  lagmax = 1
  lyfft = 150
  iccov = 30
  jccov = 5
  iym = 150
  m = 1
  iymfft = 150
  inlppc = 30
  jnlppc = 5
  nyd = 144
  ymiss0 = 1.16e0
  y1miss = 0.89e0
  y2miss = 0.89e0
!
!  test error handling
!
  200 itest = itest + 1
  n = 0
  lagmax = 1
  lyfft = 0
  iccov = 0
  jccov = 0
  iym = 0
  m = 0
  iymfft = 0
  inlppc = 0
  jnlppc = 0
  nyd = 0
!
!  test error handling
!
  300 itest = itest + 1
  n = 100
  lagmax = 100
  lyfft = 0
  iccov = 0
  jccov = 0
  iym = 0
  m = 0
  iymfft = 0
  inlppc = 0
  jnlppc = 0
  nyd = 144
  ldstak = 0

  400 return

 1020 format (' output suppressed')
 1040 format (' with missing values')
 1050 format (' without missing values')
 1060 format ('test of ccf')
 1070 format ('test of ccfm')
 1080 format ('test of ccfs')
 1100 format ('test of ccff')
 1140 format ('test of ccfms')
 1150 format ('test of ccffs')
end
