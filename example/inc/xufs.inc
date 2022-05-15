subroutine xufs ( ldstak ,typ)
use M_starpac, only : acfms,acfs,iprint,nrand,scopy,setrv,ufs,ufsf,ufsfs,ufsm, &
     ufsms,ufsmv,ufsmvs,ufss,ufsv,ufsvs
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XUFS tests the Fourier spectrum analysis routines.
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
!     real acov(101)
!        the autocovariance vector.
!     real amiss
!         the missing value code for the returned acvf estimates.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     real fmax, fmin
!        the maximum and minimum frequencies at which the
!        spectrum is to be computed.
!     real freq(300)
!        the vector of frequencies at which the spectrum is computed.
!     integer i
!        an index variable
!     integer iar
!        the order of the autoregressive model to be used.
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list
!        if ierr == 0, no errors were detected
!        if ierr == 1, errors have been detected
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
!     integer lds
!        the length of the vector dstak in common cstak.
!     integer lyfft
!        the length of the vector yfft.
!     integer nf
!        the number of frequencies at which the spectrum is
!        to be computed.
!     integer nlppa(101)
!        the numbers of lagged product pairs used for each acvf.
!     integer nprt
!        a code used to specify the type of plot, where if
!        nprt = 0 the plot is suppressed, if
!        nprt = 2 the plot is decibels/linear, if
!        nprt = 2 the plot is log/linear, if
!        nprt = 3 the plot is decibels/log, and if
!        nprt = 4 the plot is log/log.
!     integer nw
!        the number of different lag window truncation points specified,
!        and therefore, the number of plots.
!     integer ny
!        the number of observations in the series y.
!     real phi(100)
!        the vector of the order iar autoregressive model coefficients.
!     real spcf(101, 4)
!        the arrays in which the fourier spectrum is stored
!        for each lag window.
!     real y(150)
!         the array containing the time series from jenkins and watts.
!     real yfft(400)
!        the vector of the observed time series to be analyzed using
!        the fft.
!     real ymiss
!        the user supplied code which is used to determine whether or
!        not an observation in the series is missing.  if y(i) = ymiss,
!        the value is assumed missing, otherwise it is not.
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
     amiss,fmax,fmin,ymiss
  integer &
     i,iar,ispcf,j,lacov,lagmax,lds,lyfft,nf,nprt,nw,ny
!
!  local arrays
  real &
     acov(101),freq(300),phi(100),spcf(101,4),y(150),yfft(400)
  integer &
     lags(4),nlppa(101)
!
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr

  data   y(  1), y(  2), y(  3), y(  4), y(  5), y(  6) &
      /-0.88e0, -0.12e0, -0.89e0, -1.38e0, -0.07e0,  1.03e0/
  data   y(  7), y(  8), y(  9), y( 10), y( 11), y( 12) &
      / 2.14e0,  0.35e0, -1.10e0, -1.78e0, -2.76e0, -1.77e0/
  data   y( 13), y( 14), y( 15), y( 16), y( 17), y( 18) &
      / 0.98e0,  1.00e0, -0.70e0, -1.01e0, -1.30e0, -0.85e0/
  data   y( 19), y( 20), y( 21), y( 22), y( 23), y( 24) &
      /-0.46e0,  1.63e0,  0.06e0, -0.17e0, -1.01e0, -1.04e0/
  data   y( 25), y( 26), y( 27), y( 28), y( 29), y( 30) &
      /-0.66e0, -1.12e0, -0.51e0, -0.71e0, -0.20e0, -0.13e0/
  data   y( 31), y( 32), y( 33), y( 34), y( 35), y( 36) &
      / 0.14e0,  1.59e0, -0.76e0, -1.08e0, -1.77e0, -1.20e0/
  data   y( 37), y( 38), y( 39), y( 40), y( 41), y( 42) &
      / 0.45e0, -0.07e0, -0.63e0, -0.35e0, -0.87e0, -0.62e0/
  data   y( 43), y( 44), y( 45), y( 46), y( 47), y( 48) &
      / 0.28e0,  1.90e0,  2.14e0,  1.05e0,  0.31e0,  1.07e0/
  data   y( 49), y( 50) &
      / 2.67e0,  2.44e0/

!     check error handling
!
!        test 1  -  miscellaneous error checking
!
  write ( *, 2000)
  ymiss = 1.16e0
  lagmax = -1
  ny = -10
  lacov = 101
  lyfft = -11
  nw = -1
  nf = -5
  fmin = 0.5e0
  fmax = 0.0e0
  nprt = -1
  ispcf = -20
  lds = 0
  write ( *, 1001)
  call ufs (y, ny)
  write ( *, 1002) ierr
  write ( *, 1003)
  call ufss(y, ny, nw, lags, nf, fmin, fmax, nprt, spcf, &
     ispcf, freq, lds)
  write ( *, 1002) ierr
  write ( *, 1019)
  call ufsf (yfft, ny, lyfft, lds)
  write ( *, 1002) ierr
  write ( *, 1020)
  call ufsfs(yfft, ny, lyfft, lds, nw, lags, nf, fmin, fmax, nprt, &
     spcf, ispcf, freq)
  write ( *, 1002) ierr
  write ( *, 1005)
  call ufsm (y, ymiss, ny)
  write ( *, 1002) ierr
  write ( *, 1006)
  call ufsms(y, ymiss, ny, nw, lags, nf, fmin, fmax, nprt, spcf, &
     ispcf, freq, lds)
  write ( *, 1002) ierr
  write ( *, 1007)
  call ufsv(acov, lagmax, ny)
  write ( *, 1002) ierr
  write ( *, 1008)
  call ufsvs (acov, lagmax, ny, nw, lags, nf, &
     fmin, fmax, nprt, spcf, ispcf, freq, lds)
  write ( *, 1002) ierr
  write ( *, 1021)
  call ufsmv(acov, nlppa, lagmax, ny)
  write ( *, 1002) ierr
  write ( *, 1022)
  call ufsmvs (acov, nlppa, lagmax, ny, nw, lags, nf, &
     fmin, fmax, nprt, spcf, ispcf, freq, lds)
  write ( *, 1002) ierr
!
!  test 2  -  miscellaneous error checking (continued)
!
  write ( *, 2010)
  ymiss = 1.16e0
  ny = 50
  lagmax = 55
  lyfft = -11
  nw = 2
  lags(1) = 0
  lags(2) = 50
  nf = 51
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 3
  ispcf = 20
  lds = 0
  write ( *, 1003)
  call ufss(y, ny, nw, lags, nf, fmin, fmax, nprt, spcf, &
     ispcf, freq, lds)
  write ( *, 1002) ierr
  write ( *, 1019)
  call ufsf (yfft, ny, lyfft, lds)
  write ( *, 1002) ierr
  write ( *, 1020)
  call ufsfs(yfft, ny, lyfft, lds, nw, lags, nf, fmin, fmax, nprt, &
     spcf, ispcf, freq)
  write ( *, 1002) ierr
  write ( *, 1006)
  call ufsms(y, ymiss, ny, nw, lags, nf, fmin, fmax, nprt, spcf, &
     ispcf, freq, lds)
  write ( *, 1002) ierr
  write ( *, 1008)
  call ufsvs (acov, lagmax, ny, nw, lags, nf, &
     fmin, fmax, nprt, spcf, ispcf, freq, lds)
  write ( *, 1002) ierr
  write ( *, 1022)
  call ufsmvs (acov, nlppa, lagmax, ny, nw, lags, nf, &
     fmin, fmax, nprt, spcf, ispcf, freq, lds)
  write ( *, 1002) ierr
!
!  test 3  -  lds too small
!
  write ( *, 2030)
  ymiss = 1.16e0
  ny = 50
  lagmax = 49
  lyfft = 400
  nw = 2
  lags(1) = 0
  lags(2) = 50
  nf = 51
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 3
  ispcf = 101
  lds = 0
  write ( *, 1003)
  call ufss(y, ny, nw, lags, nf, fmin, fmax, nprt, spcf, &
     ispcf, freq, lds)
  write ( *, 1002) ierr
  write ( *, 1019)
  call ufsf (yfft, ny, lyfft, lds)
  write ( *, 1002) ierr
  write ( *, 1020)
  call ufsfs(yfft, ny, lyfft, lds, nw, lags, nf, fmin, fmax, nprt, &
     spcf, ispcf, freq)
  write ( *, 1002) ierr
  write ( *, 1006)
  call ufsms(y, ymiss, ny, nw, lags, nf, fmin, fmax, nprt, spcf, &
     ispcf, freq, lds)
  write ( *, 1002) ierr
  write ( *, 1008)
  call ufsvs (acov, lagmax, ny, nw, lags, nf, &
     fmin, fmax, nprt, spcf, ispcf, freq, lds)
  write ( *, 1002) ierr
  write ( *, 1022)
  call ufsmvs (acov, nlppa, lagmax, ny, nw, lags, nf, &
     fmin, fmax, nprt, spcf, ispcf, freq, lds)
  write ( *, 1002) ierr
!
!  test 4  -  all data and covariances missing
!
  write ( *, 2040)
  ymiss = 1.16e0
  ny = 50
  lagmax = 49
  lyfft = 400
  nw = 2
  lags(1) = 8
  lags(2) = 16
  nf = 51
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 3
  ispcf = 101
  lds = 700
  call setrv(yfft, ny, ymiss)
  call setrv(acov, lagmax, 0.0e0)
  nlppa(1:lagmax) = 0
  write ( *, 1005)
  call ufsm(yfft, ymiss, ny)
  write ( *, 1002) ierr
  write ( *, 1006)
  call ufsms(yfft, ymiss, ny, nw, lags, nf, fmin, fmax, nprt, spcf, &
     ispcf, freq, lds)
  write ( *, 1002) ierr
  write ( *, 1021)
  call ufsmv (acov, nlppa, lagmax, ny)
  write ( *, 1002) ierr
  write ( *, 1022)
  call ufsmvs (acov, nlppa, lagmax, ny, nw, lags, nf, &
     fmin, fmax, nprt, spcf, ispcf, freq, lds)
  write ( *, 1002) ierr
!
!  test 5  -  every other value missing
!
  write ( *, 2050)
  ymiss = 1.16e0
  ny = 50
  lagmax = 49
  lyfft = 400
  nw = 2
  lags(1) = 8
  lags(2) = 16
  nf = 51
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 3
  ispcf = 101
  lds = 700
  call setrv(yfft, ny, ymiss)
  yfft(1:ny:2) = y(1:ny:2)
  write ( *, 1005)
  call ufsm(yfft, ymiss, ny)
  write ( *, 1002) ierr
  write ( *, 1006)
  call ufsms(yfft, ymiss, ny, nw, lags, nf, fmin, fmax, nprt, spcf, &
     ispcf, freq, lds)
  write ( *, 1002) ierr
!
!  check results from valid call
!
  ymiss = 1.16e0
  ny = 50
  lagmax = 49
  lyfft = 400
  nw = 2
  lags(1) = 8
  lags(2) = 16
  nf = 51
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 2
  ispcf = 101
  lds = ldstak
!
!  test of ufs
!
  write ( *, 2020)
  write ( *, 1001)
  call ufs (y, ny)
  write ( *, 1002) ierr
!
!  test of ufss
!
  write ( *, 2020)
  write ( *, 1003)
  call ufss(y, ny, nw, lags, nf, fmin, fmax, nprt, spcf, &
     ispcf, freq, lds)
  write ( *, 1002) ierr
!
!  print returned values from ufss
!
  write ( *, 1004) (freq(i), (spcf(i,j),j=1,nw), i=1,nf)
!
!  test of ufsf
!
  write ( *, 2020)
  write ( *, 1019)
  call scopy(ny, y, 1, yfft, 1)
  call ufsf (yfft, ny, lyfft, lds)
  write ( *, 1002) ierr
!
!  test of ufsfs
!
  write ( *, 2020)
  write ( *, 1020)
  call scopy(ny, y, 1, yfft, 1)
  call ufsfs(yfft, ny, lyfft, lds, nw, lags, nf, fmin, fmax, nprt, &
     spcf, ispcf, freq)
  write ( *, 1002) ierr
!
!  print returned values from ufsfs
!
  write ( *, 1004) (freq(i), (spcf(i,j),j=1,nw), i=1,nf)
!
!  test of ufsm
!
  write ( *, 2020)
  write ( *, 1005)
  call ufsm (y, ymiss, ny)
  write ( *, 1002) ierr
!
!     test of ufsms
!
  write ( *, 2020)
  write ( *, 1006)
  call ufsms(y, ymiss, ny, nw, lags, nf, fmin, fmax, nprt, spcf, &
     ispcf, freq, lds)
  write ( *, 1002) ierr
!
!     print returned values from ufsms
!
  write ( *, 1004) (freq(i), (spcf(i,j),j=1,nw), i=1,nf)
!
!     test of ufsv
!
  write ( *, 2020)
  call acfs (y, ny, lagmax, lacov, acov, iar, phi, 0, lds)
  write ( *, 1007)
  call ufsv(acov, lagmax, ny)
  write ( *, 1002) ierr
!
!     test of ufsvs
!
  write ( *, 2020)
  write ( *, 1008)
  call ufsvs (acov, lagmax, ny, nw, lags, nf, &
     fmin, fmax, nprt, spcf, ispcf, freq, lds)
  write ( *, 1002) ierr
!
!     print returned values from ufsvs
!
  write ( *, 1004) (freq(i), (spcf(i,j),j=1,nw), i=1,nf)
!
!     test of ufsmv
!
  write ( *, 2020)
  call acfms (y, ymiss, ny, lagmax, lacov, acov, amiss, nlppa, &
     0, lds)
  write ( *, 1021)
  call ufsmv(acov, nlppa, lagmax, ny)
  write ( *, 1002) ierr
!
!     test of ufsmvs
!
  write ( *, 2020)
  write ( *, 1022)
  call ufsmvs (acov, nlppa, lagmax, ny, nw, lags, nf, &
     fmin, fmax, nprt, spcf, ispcf, freq, lds)
  write ( *, 1002) ierr
!
!  print returned values from ufsmvs
!
  write ( *, 1004) (freq(i), (spcf(i,j),j=1,nw), i=1,nf)
!
!  minimum problem size
!
  ymiss = 1.16e0
  ny = 17
  lagmax = 1
  lyfft = 400
  nw = 2
  lags(1) = 1
  lags(2) = 16
  nf = 1
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 2
  ispcf = 101
  lds = 700
!
!  test of ufs
!
  write ( *, 2060)
  write ( *, 1001)
  call ufs (y, ny)
  write ( *, 1002) ierr
!
!  test of ufss
!
  write ( *, 2060)
  write ( *, 1003)
  call ufss(y, ny, nw, lags, nf, fmin, fmax, nprt, spcf, &
     ispcf, freq, lds)
  write ( *, 1002) ierr
!
!  print returned values from ufss
!
  write ( *, 1004) (freq(i), (spcf(i,j),j=1,nw), i=1,nf)
!
!  check handling of fmin and fmax
!
  ymiss = 1.16e0
  ny = 50
  lagmax = 49
  lyfft = 400
  nw = 2
  lags(1) = 8
  lags(2) = 16
  nf = 51
  fmin = 0.45e0
  fmax = 0.5e0
  nprt = 2
  ispcf = 101
  lds = 700
!
!  test of ufss
!
  write ( *, 2070)
  write ( *, 1003)
  call ufss(y, ny, nw, lags, nf, fmin, fmax, nprt, spcf, &
     ispcf, freq, lds)
  write ( *, 1002) ierr
!
!     print returned values from ufss
!
  write ( *, 1004) (freq(i), (spcf(i,j),j=1,nw), i=1,nf)
!
!  white noise spectrum
!
  ymiss = 1.16e0
  call nrand(yfft, ny, 12345)
  ny = 50
  lagmax = 49
  lyfft = 400
  nw = 2
  lags(1) = 8
  lags(2) = 16
  nf = 51
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 2
  ispcf = 101
  lds = 700
!
!  test of ufss
!
  write ( *, 2080)
  write ( *, 1003)
  call ufss(yfft, ny, nw, lags, nf, fmin, fmax, nprt, spcf, &
     ispcf, freq, lds)
  write ( *, 1002) ierr
!
!  print returned values from ufss
!
  write ( *, 1004) (freq(i), (spcf(i,j),j=1,nw), i=1,nf)

  return

 1001 format (' test of ufs')
 1002 format (/' ierr is', i5/)
 1003 format (' test of ufss')
 1004 format (3(1x, e16.8))
 1005 format (' test of ufsm')
 1006 format (' test of ufsms')
 1007 format (' test of ufsv')
 1008 format (' test of ufsvs')
 1019 format (' test of ufsf')
 1020 format (' test of ufsfs')
 1021 format (' test of ufsmv')
 1022 format (' test of ufsmvs')
 2000 format ('1check error handling  -  test 1')
 2010 format ('1check error handling  -  test 2')
 2020 format ('1valid problem')
 2030 format ('1lds too small')
 2040 format ('1all data and covariances missing')
 2050 format ('1every other data value missing')
 2060 format ('1minimum problem size')
 2070 format ('1check handling of fmin and fmax')
 2080 format ('1white noise spectrum')
end
