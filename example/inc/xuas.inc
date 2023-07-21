subroutine xuas ( ldstak ,typ)
use M_starpac, only : acfs,iprint,scopy,setrv,uas,uasf,uasfs,uass,uasv,uasvs
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XUAS tests the autoregressive spectrum analysis routines.
!
!  Discussion:
!
!    Series y is the first 50 values of the series listed on page
!    318 of jenkins and watts.  the spectrum of this series is shown
!    for various bandwidth on page 270 of jenkins and watts.
!
!  Modified:
!
!    02 December 2006
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
!     real acov(101)
!        the autocovariance vector.
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
!     integer lacov
!        the length of the vector acov.
!     integer lag, lagmax
!        the indexing variable indicating the lag value of the
!        autocovariance being computed and the maximum lag to be used,
!        respectively.
!     integer lds
!        the length of the vector dstak in common cstak.
!     integer lyfft
!        the length of the vector yfft.
!     integer nf
!        the number of frequencies at which the spectrum is
!        to be computed.
!     integer nprt
!        a code used to specify the type of plot, where if
!        nprt < 0 the plot is decibles/linear
!        nprt = 0 the plot is suppressed
!        nprt > 0 the plot is log/linear
!     integer ny
!        the number of observations in the series y.
!     real phi(101)
!        the vector of the order iar autoregressive model coefficients.
!     real spca(101)
!        the arrays in which the autoregressive spectrum is stored
!        for each lag window.
!     real spcf(101)
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

  real acov(101)
  double precision dstak(12)
  real fmax
  real fmin
  real freq(300)
  integer i
  integer iar
  integer ierr
  integer lacov
  integer lag
  integer lagmax
  integer lds
  integer ldstak
  integer lyfft
  integer nf
  integer nprt
  integer ny
  real phi(101)
  real spca(101)
  real spcf(101)
  real y(150)
  real yfft(400)
  real ymiss


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

!
!  Check error handling
!
!  Test 1  -  miscellaneous error checking
!
  write ( *, 2000)
  ymiss = 1.16e0
  lagmax = -1
  ny = -10
  lacov = 101
  lag = -2
  iar = -2
  lyfft = -11
  nf = -5
  fmin = 0.5e0
  fmax = 0.0e0
  nprt = -1
  lds = 0
  write ( *, 1001)
  call uas(y, ny)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test 1 of UASS:'
  write ( *, '(a)' ) ' '

  call uass(y, ny, iar, phi, lagmax, lag, nf, &
     fmin, fmax, nprt,spca, spcf, freq, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr

  write ( *, 1019)
  call uasf (yfft, ny, lyfft, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr

  write ( *, 1020)
  call uasfs(yfft, ny, lyfft, lds, iar, phi, lagmax, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr

  write ( *, 1007)
  call uasv(acov, lagmax, ny)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr

  write ( *, 1008)
  call uasvs(acov, lagmax, y, ny, iar, phi, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Test 2  -  miscellaneous error checking (continued)
!
  write ( *, 2010)
  ymiss = 1.16e0
  ny = 50
  lagmax = 50
  lag = 101
  iar = 101
  call setrv(phi, iar, 2.0e0)
  call setrv(acov, lagmax+1, 2.0e0)
  acov(1) = 1.0e0
  lyfft = -11
  nf = 51
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 5
  lds = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test 2 of UASS:'
  write ( *, '(a)' ) ' '

  call uass(y, ny, iar, phi, lagmax, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr

  write ( *, 1019)
  call uasf (yfft, ny, lyfft, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr

  write ( *, 1020)
  call uasfs(yfft, ny, lyfft, lds, iar, phi, lagmax, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr

  write ( *, 1007)
  call uasv(acov, lagmax, ny)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr

  write ( *, 1008)
  call uasvs(acov, lagmax, y, ny, iar, phi, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Test 3  -  lds too small
!
  write ( *, 2030)
  ymiss = 1.16e0
  ny = 50
  lagmax = 49
  lyfft = 400
  lag = 16
  iar = 2
  phi(1) = 1.0e0
  phi(2) = -0.5e0
  call acfs (y, ny, lagmax, lacov, acov, iar, phi, 0, 700)
  nf = 51
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 3
  lds = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test 3 of UASS:'
  write ( *, '(a)' ) ' '

  call uass(y, ny, iar, phi, lagmax, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr

  write ( *, 1019)
  call uasf (yfft, ny, lyfft, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr

  write ( *, 1020)
  call uasfs(yfft, ny, lyfft, lds, iar, phi, lagmax, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr

  write ( *, 1008)
  call uasvs(acov, lagmax, y, ny, iar, phi, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Check results from valid call
!
  ymiss = 1.16e0
  ny = 50
  lagmax = 49
  lyfft = 400
  lag = 16
  iar = 2
  phi(1) = 1.0e0
  phi(2) = -0.5e0
  nf = 51
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 2
  lds = 700
!
!  Test of UAS.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  VALID PROBLEM:'
  write ( *, '(a)' ) ' '

  write ( *, 1001)
  call uas ( y, ny )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Test of UASS.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  VALID PROBLEM:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test 4 of UASS:'
  write ( *, '(a)' ) ' '

  call uass ( y, ny, iar, phi, lagmax, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq, lds )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Print returned values from UASS.
!
  write ( *, 1004) (freq(i), spca(i), spcf(i), i=1,nf)
  write ( *, 1005) iar, lag
  write ( *, 1006) (phi(i), i=1,abs(iar))
!
!  Test of uasf
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  VALID PROBLEM:'
  write ( *, '(a)' ) ' '

  write ( *, 1019)
  call scopy(ny, y, 1, yfft, 1)
  call uasf (yfft, ny, lyfft, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Test of uasfs
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  VALID PROBLEM:'
  write ( *, '(a)' ) ' '

  write ( *, 1020)
  call scopy(ny, y, 1, yfft, 1)
  call uasfs(yfft, ny, lyfft, lds, iar, phi, lagmax, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Print returned values from uasfs
!
  write ( *, 1004) (freq(i), spca(i), spcf(i), i=1,nf)
  write ( *, 1005) iar, lag
  write ( *, 1006) (phi(i), i=1,abs(iar))
!
!  Test of uasv
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  VALID PROBLEM:'
  write ( *, '(a)' ) ' '

  write ( *, 1007)
  call uasv(acov, lagmax, ny)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Test of uasvs
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  VALID PROBLEM:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test 5 of UASS:'
  write ( *, '(a)' ) ' '

  call uasvs(acov, lagmax, y, ny, iar, phi, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Print returned values from uasvs
!
  write ( *, 1004) (freq(i), spca(i), spcf(i), i=1,nf)
  write ( *, 1005) iar, lag
  write ( *, 1006) (phi(i), i=1,abs(iar))
!
!  Minimum problem size
!
  ymiss = 1.16e0
  ny = 17
  lagmax = 1
  lyfft = 400
  lag = 1
  iar = -1
  phi(1) = 1.0e0
  phi(2) = -0.5e0
  nf = 1
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 2
  lds = ldstak
!
!  Test of uas
!
  write ( *, 2060)
  write ( *, 1001)
  call uas(y, ny)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Test of UASS.
!
  write ( *, 2060)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test 6 of UASS:'
  write ( *, '(a)' ) ' '

  call uass(y, ny, iar, phi, lagmax, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Print returned values from UASS.
!
  write ( *, 1004) (freq(i), spca(i), spcf(i), i=1,nf)
  write ( *, 1005) iar, lag
  write ( *, 1006) (phi(i), i=1,abs(iar))
!
!  Check handling of fmin and fmax, and lag==0 and iar==0
!
  ny = 50
  lagmax = 49
  lyfft = 400
  lag = 0
  iar = 0
  phi(1) = 1.0e0
  phi(2) = -0.5e0
  nf = 51
  fmin = 0.45e0
  fmax = 0.5e0
  nprt = 2
!
!  Test of UASS.
!
  write ( *, 2070)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test 7 of UASS:'
  write ( *, '(a)' ) ' '

  call uass(y, ny, iar, phi, lagmax, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Print returned values from UASS.
!
  write ( *, 1004) (freq(i), spca(i), spcf(i), i=1,nf)
  write ( *, 1005) iar, lag
  write ( *, 1006) (phi(i), i=1,abs(iar))
!
!  White noise spectrum
!
  ymiss = 1.16e0
  call setrv(yfft, ny, 0.0e0)
  ny = 50
  lagmax = 49
  lyfft = 400
  lag = 16
  iar = 2
  phi(1) = 1.0e0
  phi(2) = -0.5e0
  nf = 51
  fmin = 0.0e0
  fmax = 0.5e0
  nprt = 2
!
!  Test of uass
!
  write ( *, 2080)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test 8 of UASS:'
  write ( *, '(a)' ) ' '

  call uass(yfft, ny, iar, phi, lagmax, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Print returned values from UASS.
!
  write ( *, 1004) (freq(i), spca(i), spcf(i), i=1,nf)
  write ( *, 1005) iar, lag
  write ( *, 1006) (phi(i), i=1,abs(iar))
!
!  Suppress output and check handling of lag .lt.0 and iar .lt. 0
!
  ny = 50
  lagmax = 49
  lyfft = 400
  lag = 0
  iar = 0
  phi(1) = 1.0e0
  phi(2) = -0.5e0
  nf = 51
  fmin = 0.45e0
  fmax = 0.5e0
  nprt = 0
!
!  Test of UASS.
!
  write ( *, 2090)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Test 9 of UASS:'
  write ( *, '(a)' ) ' '

  call uass(y, ny, iar, phi, lagmax, lag, nf, &
     fmin, fmax, nprt, spca, spcf, freq, lds)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Return value of IERR is ', ierr
!
!  Print returned values from UASS.
!
  write ( *, 1004) (freq(i), spca(i), spcf(i), i=1,nf)
  write ( *, 1005) iar, lag
  write ( *, 1006) (phi(i), i=1,abs(iar))

  return

 1001 format ('test of uas')
 1004 format (3(1x, e16.8))
 1005 format (/' iar = ', i5/' lag = ', i5)
 1006 format (/' phi = ', (1x, 5e21.8))
 1007 format ('test of uasv')
 1008 format ('test of uasvs')
 1019 format ('test of uasf')
 1020 format ('test of uasfs')
 2000 format ('check error handling  -  test 1')
 2010 format ('check error handling  -  test 2')
 2030 format ('lds too small')
 2060 format ('minimum problem size')
 2070 format ( &
     'check handling of fmin and fmax, lag and iar equal to zero')
 2080 format ('white noise spectrum')
 2090 format ('suppress output, lag and iar less than zero')
end
