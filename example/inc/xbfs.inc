subroutine xbfs ( lds ,typ)
use M_starpac, only : bfs,bfsf,bfsfs,bfsm,bfsms,bfsmv,bfsmvs,bfss,bfsv,bfsvs, ccfms,ccfs,iprint,nrand,scopy,setiv,setrv
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XBFS tests the Fourier spectrum analysis routines.
!
!  Discussion:
!
!    series y1 and y2 are listed as series x1 and x2 on page of 361 of
!    jenkins and watts.
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
!     real cmiss
!         the missing value code for the returned ccvf estimates.
!     real ccov(101,2,2)
!        the covariances.
!     real cspc2(300,2)
!        the squared coherency component of the spectrum.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     real fmax, fmin
!        the maximum and minimum frequencies at which the
!        spectrum is to be computed.
!     real freq(300)
!        the vector of frequencies at which the spectrum is computed.
!     integer i
!        an index variable
!     integer iccov
!        the first dimension of the array ccov.
!     integer icspc2
!        the first dimension of the array cspc2.
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list
!        if ierr == 0, no errors were detected
!        if ierr == 1, errors have been detected
!     integer iphas
!        the first dimension of the array phas.
!     integer ispcf
!         the actual dimension for the spectrum arrays.
!     integer j
!        index variable.
!     integer lacov
!        the length of the vector acov.
!     integer lagmax
!        the indexing variable indicating the lag value of the
!        autocovariance being computed and the maximum lag to be used,
!        respectively.
!     integer lags(4)
!        the array used to store the lag window truccation
!        points used for each set of spectrum values.
!     integer lds, ldstak
!        the length of the vector dstak in common cstak.
!     integer lyfft
!        the length of the vectors yfft and zfft, respectively..
!     integer nf
!        the number of frequencies at which the spectrum is
!        to be computed.
!     integer nlppc(101,2,2)
!        the numbers of lagged product pairs used for each acvf.
!     integer nprt
!        a code used to specify the type of plot.
!        if nprt < 0 the plot is decibles/linear
!        if nprt = 0 the plot is suppressed.
!        if nprt > 0 the plot is log/linear
!     integer nw
!        the number of different lag window truncation points specified,
!        and therefore, the number of plots.
!     integer n
!        the number of observations in the series y.
!     real phas(300,2)
!        the phase component of the spectrum.
!     real yfft1(400), yfft2(400)
!        the vectors of the observed time series to be analyzed using
!        the fft.
!     real ymiss, ymiss1, ymiss2, ymmiss(4)
!        the user supplied code which is used to determine whether or
!        not an observation in the series is missing.  if y(i) = ymiss,
!        the value is assumed missing, otherwise it is not.
!     real ym(150,2)
!        the multivariate representation of the data
!     real y1(150), y2(150)
!         the vectors containing the time series from jenkins and watts.
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
     cmiss,fmax,fmin,ymiss,ymiss1,ymiss2
  integer &
     i,iccov,icspc2,index1,index2,inlppc,iphas,ispcf,j, &
     jccov,jnlppc,lacov,lagmax,ldstak,lyfft,n,nf,nprt,nw
!
!  local arrays
  real &
     ccov(101,2,2),cspc2(300,2),freq(300),phas(300,2),y1(150), &
     y2(150),yfft1(400),yfft2(400),ym(150,2),ymmiss(4)
  integer &
     lags(4),nlppc(101,2,2)
!
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr
!
!  equivalences
  equivalence (ym(1,1),y1(1))
  equivalence (ym(1,2),y2(1))

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

  call setrv(ymmiss, 4, 0.89e0)
!
!  check error handling
!
!  test 1  -  miscellaneous error checking
!
  write ( *, 2000)
  lagmax = -1
  n = -10
  index1 = 0
  index2 = 0
  iccov = 0
  jccov = 0
  inlppc = 0
  jnlppc = 0
  icspc2 = -10
  iphas = -10
  lacov = -11
  lyfft = -11
  nw = -1
  nf = -5
  fmin = 0.5e0
  fmax = 0.0e0
  nprt = -1
  ispcf = -20
  ldstak = 0
  ymiss1 = 0.89e0
  ymiss2 = 0.89e0

  write ( *, 1001)
  call bfs (y1, y2, n)
  write ( *, 1002) ierr

  write ( *, 1003)
  call bfss(y1, y2, n, nw, lags, nf, fmin, fmax, nprt, &
     cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr

  write ( *, 1019)
  call bfsf (yfft1, yfft2, n, lyfft, ldstak)
  write ( *, 1002) ierr

  write ( *, 1020)
  call bfsfs(yfft1, yfft2, n, lyfft, ldstak, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq)
  write ( *, 1002) ierr

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of BFSM'
  write ( *, '(a)' ) ' '

  call bfsm (y1, ymiss1, y2, ymiss2, n)
  write ( *, 1002) ierr

  write ( *, 1006)
  call bfsms(y1, ymiss1, y2, ymiss2, n, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr

  write ( *, 1007)
  call bfsv(ccov, index1, index2, n, lagmax, iccov, jccov)
  write ( *, 1002) ierr

  write ( *, 1008)
  call bfsvs (ccov, index1, index2, n, iccov, jccov, &
     nw, lags, nf, fmin, fmax, nprt, cspc2, icspc2, phas, iphas, &
     freq, ldstak)
  write ( *, 1002) ierr

  write ( *, 1021)
  call bfsmv(ccov, nlppc, index1, index2, n, lagmax, iccov, jccov, &
     inlppc, jnlppc)
  write ( *, 1002) ierr

  write ( *, 1022)
  call bfsmvs (ccov, nlppc, index1, index2, n, iccov, jccov, &
     inlppc, jnlppc, nw, lags, nf, fmin, fmax, nprt, cspc2, icspc2, &
     phas, iphas, freq, ldstak)
  write ( *, 1002) ierr
!
!  test 2  -  miscellaneous error checking (continued)
!
  write ( *, 2010)
  n = 100
  lagmax = 40
  index1 = 0
  index2 = 0
  iccov = 0
  jccov = 0
  inlppc = 0
  jnlppc = 0
  icspc2 = 300
  iphas = 300
  lacov = 101
  lyfft = -11
  nw = 2
  lags(1) = 0
  lags(2) = 100
  nf = 202
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 3
  ispcf = 101
  ldstak = 0

  write ( *, 1003)
  call bfss(y1, y2, n, nw, lags, nf, fmin, fmax, nprt, &
     cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr

  write ( *, 1019)
  call bfsf (yfft1, yfft2, n, lyfft, ldstak)
  write ( *, 1002) ierr

  write ( *, 1020)
  call bfsfs(yfft1, yfft2, n, lyfft, ldstak, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq)
  write ( *, 1002) ierr

  write ( *, 1006)
  call bfsms(y1, ymiss1, y2, ymiss2, n, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr

  write ( *, 1008)
  call bfsvs (ccov, index1, index2, n, iccov, jccov, &
     nw, lags, nf, fmin, fmax, nprt, cspc2, icspc2, phas, iphas, &
     freq, ldstak)
  write ( *, 1002) ierr

  write ( *, 1022)
  call bfsmvs (ccov, nlppc, index1, index2, n, iccov, jccov, &
     inlppc, jnlppc, nw, lags, nf, fmin, fmax, nprt, cspc2, icspc2, &
     phas, iphas, freq, ldstak)
  write ( *, 1002) ierr
!
!  test 3  -  ldstak too small
!
  write ( *, 2030)
  n = 100
  index1 = 2
  index2 = 1
  iccov = 101
  jccov = 2
  inlppc = 101
  jnlppc = 2
  icspc2 = 300
  iphas = 300
  lagmax = 99
  lacov = 101
  lyfft = 400
  nw = 2
  lags(1) = 8
  lags(2) = 16
  nf = 26
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 1
  ispcf = 101
  ldstak = 0

  write ( *, 1003)
  call bfss(y1, y2, n, nw, lags, nf, fmin, fmax, nprt, &
     cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr

  write ( *, 1019)
  call bfsf (yfft1, yfft2, n, lyfft, ldstak)
  write ( *, 1002) ierr

  write ( *, 1020)
  call bfsfs(yfft1, yfft2, n, lyfft, ldstak, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq)
  write ( *, 1002) ierr

  write ( *, 1006)
  call bfsms(y1, ymiss1, y2, ymiss2, n, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr

  write ( *, 1008)
  call bfsvs (ccov, index1, index2, n, iccov, jccov, &
     nw, lags, nf, fmin, fmax, nprt, cspc2, icspc2, phas, iphas, &
     freq, ldstak)
  write ( *, 1002) ierr

  write ( *, 1022)
  call bfsmvs (ccov, nlppc, index1, index2, n, iccov, jccov, &
     inlppc, jnlppc, nw, lags, nf, fmin, fmax, nprt, cspc2, icspc2, &
     phas, iphas, freq, ldstak)
  write ( *, 1002) ierr
!
!  test 4  -  all data and covariances missing
!
  write ( *, 2040)
  n = 100
  lagmax = 99
  icspc2 = 300
  iphas = 300
  lacov = 101
  lyfft = 400
  nw = 2
  lags(1) = 8
  lags(2) = 16
  nf = 26
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 1
  ispcf = 101
  ldstak = lds
  call setrv(yfft1, n, ymiss1)
  call setrv(yfft2, n, ymiss2)
  call setrv(ccov, 404, 0.0e0)
  call setiv(nlppc, 404, 0)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of BFSM'
  write ( *, '(a)' ) ' '

  call bfsm (yfft1, ymiss1, yfft2, ymiss2, n)
  write ( *, 1002) ierr

  write ( *, 1006)
  call bfsms(yfft1, ymiss1, yfft2, ymiss2, n, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr

  write ( *, 1021)
  call bfsmv(ccov, nlppc, index1, index2, n, lagmax, iccov, jccov, &
     inlppc, jnlppc)
  write ( *, 1002) ierr

  write ( *, 1022)
  call bfsmvs (ccov, nlppc, index1, index2, n, iccov, jccov, &
     inlppc, jnlppc, nw, lags, nf, fmin, fmax, nprt, cspc2, icspc2, &
     phas, iphas, freq, ldstak)
  write ( *, 1002) ierr
!
!  test 5  -  every other value missing
!
  write ( *, '(a)' ) '  Every other data value missing.'
  n = 100
  lagmax = 99
  icspc2 = 300
  iphas = 300
  lacov = 101
  lyfft = 400
  nw = 2
  lags(1) = 8
  lags(2) = 16
  nf = 26
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 1
  ispcf = 101
  ldstak = lds
  call setrv(yfft1, n, ymiss1)
  call setrv(yfft2, n, ymiss2)
  yfft1(1:n:2) = y1(1:n:2)
  yfft2(1:n:2) = y2(1:n:2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of BFSM'
  write ( *, '(a)' ) ' '

  call bfsm ( yfft1, ymiss1, yfft2, ymiss2, n )
  write ( *, 1002) ierr

  write ( *, 1006)
  call bfsms(yfft1, ymiss1, yfft2, ymiss2, n, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr
!
!  check results from valid call
!
  write ( *, 2020)
  ymiss = 1.16e0
  n = 100
  lagmax = 99
  icspc2 = 300
  iphas = 300
  lacov = 101
  lyfft = 400
  nw = 2
  lags(1) = 8
  lags(2) = 16
  nf = 26
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 1
  ispcf = 101
  ldstak = lds
!
!  test of bfs
!
  write ( *, 1001)
  call bfs (y1, y2, n)
  write ( *, 1002) ierr
!
!  test of bfss
!
  write ( *, 2020)
  write ( *, 1003)
  call bfss(y1, y2, n, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr
!
!  print returned values from bfss
!
  write ( *, 1004) (freq(i), (cspc2(i,j), phas(i,j),j=1,nw), &
     i=1,nf)
!
!  test of bfsf
!
  write ( *, 2020)
  write ( *, 1019)
  call scopy(n, y1, 1, yfft1, 1)
  call scopy(n, y2, 1, yfft2, 1)
  call bfsf (yfft1, yfft2, n, lyfft, ldstak)
  write ( *, 1002) ierr
!
!  test of bfsfs
!
  write ( *, 2020)
  write ( *, 1020)
  call scopy(n, y1, 1, yfft1, 1)
  call scopy(n, y2, 1, yfft2, 1)
  call bfsfs(yfft1, yfft2, n, lyfft, ldstak, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq)
  write ( *, 1002) ierr
!
!  print returned values from bfsfs
!
  write ( *, 1004) (freq(i), (cspc2(i,j), phas(i,j),j=1,nw), &
     i=1,nf)
!
!  test of bfsm
!
  write ( *, 2020)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test of BFSM'
  write ( *, '(a)' ) ' '

  call bfsm (y1, ymiss1, y2, ymiss2, n)
  write ( *, 1002) ierr
!
!  test of bfsms
!
  write ( *, 2020)
  write ( *, 1006)
  call bfsms(y1, ymiss1, y2, ymiss2, n, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr
!
!  print returned values from bfsms
!
  write ( *, 1004) (freq(i), (cspc2(i,j), phas(i,j),j=1,nw), &
     i=1,nf)
!
!  test of bfsv
!
  write ( *, 2020)
  call ccfs (ym, n, 2, 150, lagmax, ccov, iccov, jccov, 0, &
     ldstak)
  write ( *, 1007)
  call bfsv(ccov, index1, index2, n, lagmax, iccov, jccov)
  write ( *, 1002) ierr
!
!  test of bfsvs
!
  write ( *, 2020)
  write ( *, 1008)
  call bfsvs(ccov, index1, index2, n, iccov, jccov, nw, lags, nf, &
     fmin, fmax, nprt, cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr
!
!  print returned values from bfsvs
!
  write ( *, 1004) (freq(i), (cspc2(i,j), phas(i,j),j=1,nw), &
     i=1,nf)
!
!  test of bfsmv
!
  write ( *, 2020)
  call ccfms (ym, ymmiss, n, 2, 150, lagmax, ccov, cmiss, iccov, &
    jccov, nlppc, inlppc, jnlppc, 0, ldstak)
  write ( *, 1021)
  call bfsmv(ccov, nlppc, index1, index2, n, lagmax, iccov, &
    jccov, inlppc, jnlppc)
  write ( *, 1002) ierr
!
!  test of bfsmvs
!
  write ( *, 2020)
  write ( *, 1022)
  call bfsmvs(ccov, nlppc, index1, index2, n, iccov, &
    jccov, inlppc, jnlppc, nw, lags, nf, fmin, fmax, nprt, &
    cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr
!
!  print returned values from bfsmvs
!
  write ( *, 1004) (freq(i), (cspc2(i,j), phas(i,j),j=1,nw), &
     i=1,nf)
!
!  minimum problem size
!
  ymiss = 1.16e0
  n = 17
  lagmax = 1
  icspc2 = 1
  iphas = 1
  lacov = 101
  lyfft = 400
  nw = 1
  lags(1) = 1
  lags(2) = 16
  nf = 1
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 1
  ispcf = 101
  ldstak = lds
!
!  test of bfs
!
  write ( *, 2060)
  write ( *, 1001)
  call bfs(y1, y2, n)
  write ( *, 1002) ierr
!
!  test of bfss
!
  write ( *, 2060)
  write ( *, 1003)
  call bfss(y1, y2, n, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr
!
!  print returned values from bfss
!
  write ( *, 1004) (freq(i), (cspc2(i,j), phas(i,j),j=1,nw), &
     i=1,nf)
!
!  check handling of fmin and fmax
!
  ymiss = 1.16e0
  n = 100
  lagmax = 99
  icspc2 = 300
  iphas = 300
  lacov = 101
  lyfft = 400
  nw = 2
  lags(1) = 8
  lags(2) = 16
  nf = 26
  fmin = 0.45e0
  fmax = 0.5e0
  nprt = 1
  ispcf = 101
  ldstak = lds
!
!  test of bfss
!
  write ( *, 2070)
  write ( *, 1003)
  call bfss(y1, y2, n, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq, ldstak)
  write ( *, 1002) ierr
!
!  print returned values from bfss
!
  write ( *, 1004) (freq(i), (cspc2(i,j), phas(i,j),j=1,nw), &
     i=1,nf)
!
!  check results for white noise spectrum
!
  ymiss = 1.16e0
  call nrand(yfft1, n, 12343)
  call nrand(yfft2, n, 34523)
  n = 100
  lagmax = 99
  icspc2 = 300
  iphas = 300
  lacov = 101
  lyfft = 400
  nw = 2
  lags(1) = 8
  lags(2) = 16
  nf = 26
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 1
  ispcf = 101
  ldstak = lds
!
!  test of bfsfs
!
  write ( *, 2080)
  write ( *, 1003)
  call bfsfs(yfft1, yfft2, n, lyfft, ldstak, nw, lags, nf, fmin, &
     fmax, nprt, cspc2, icspc2, phas, iphas, freq)
  write ( *, 1002) ierr
!
!  print returned values from bfss
!
  write ( *, 1004) (freq(i), (cspc2(i,j), phas(i,j),j=1,nw), &
     i=1,nf)

  return

 1001 format (' test of bfs')
 1002 format (' ierr is', i5/)
 1003 format (' test of bfss')
 1004 format (5(1x, e15.7))
 1006 format (' test of bfsms')
 1007 format (' test of bfsv')
 1008 format (' test of bfsvs')
 1019 format (' test of bfsf')
 1020 format (' test of bfsfs')
 1021 format (' test of bfsmv')
 1022 format (' test of bfsmvs')
 2000 format ('1check error handling  -  test 1')
 2010 format ('1check error handling  -  test 2')
 2020 format ('1valid problem')
 2030 format ('1lds too small')
 2040 format ('1all data and covariances missing')
 2060 format ('1minimum problem size')
 2070 format ('1check handling of fmin and fmax')
 2080 format ('1white noise spectrum')
end
