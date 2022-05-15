subroutine xaimt ( ldstak ,typ)
use M_starpac, only : aime,aimec,aimes,aimf,aimfs,iprint
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XAIMT tests the time series model estimation routines.
!
!  Discussion:
!
!    Series y is the airline data listed on page 531 of box and jenkins
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
!    George Box, Gwilym Jenkins,
!    Time Series Analysis: Forecasting and Control,
!    Holden-Day, 1970,
!    QA280.B67.
!
!  Parameters:
!
!     real delta
!        the maximum change allowed in the model parameters at the
!        first iteration.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     real fcst(50,5)
!        the forecasts.
!     real fcstsd(50)
!        the standard deviations of the forecasts.
!     integer i
!        *
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr .ge. 1, errors were detected.
!     integer ifcst
!        the first dimension of the array fcst.
!     integer ifixed(50)
!        the indicator values used to designate whether the
!        parameters are to be optimized or are to be held fixed.  if
!        ifixed(i).ne.0, then par(i) will be optimized.  if
!        ifixed(i)==0, then par(i) will be held fixed.
!        ifixed(i).lt.0, then all par(i),i=1,npar, will be optimized..
!     integer ivaprx
!        an indicator value used to designate which option is to be used
!        to compute the variance covariance matrix (vcv), where
!        ivaprx le 0 indicates the the default option will be used
!        ivaprx eq 1 indicates the vcv is to be computed by
!                       inverse(hessian)*transpose(jacobian)*jacobian
!                          *inverse(hessian)
!                    using both the model subroutine the user supplied
!                    derivative subroutine when it is available
!        ivaprx eq 2 indicates the vcv is to be computed by
!                       inverse(hessian)
!                    using both the model subroutine the user supplied
!                    derivative subroutine when it is available
!        ivaprx eq 3 indicates the vcv is to be computed by
!                       inverse(transpose(jacobian)*jacobian)
!                    using both the model subroutine the user supplied
!                    derivative subroutine when it is available
!        ivaprx eq 4 indicates the vcv is to be computed by
!                       inverse(hessian)*transpose(jacobian)*jacobian
!                          *inverse(hessian)
!                    using only the model subroutine
!        ivaprx eq 5 indicates the vcv is to be computed by
!                       inverse(hessian)
!                    using only the model subroutine
!        ivaprx eq 6 indicates the vcv is to be computed by
!                       inverse(transpose(jacobian)*jacobian)
!                    using only the model subroutine
!        ivaprx ge 7 indicates the default option will be used
!     integer ivcv
!        the first dimension of the variance covariance matrix vcv.
!     integer ldstak
!        the length of the array dstak.
!     integer mit
!        the maximum number of iterations allowed.
!     integer mspec(4,50)
!        the array containing the values of p, d, q, and s for each
!        factor.
!     integer nfac
!        the number of factors in the model
!     integer npar
!        the number of unknown parameters in the model.
!     integer npare
!        the number of parameters estimated by the routine.
!     integer nprt
!        the parameter used to indicate how much printed output is
!        to be provided.
!     integer ny
!        the number of observations.
!     real par(50)
!        the array in which the current estimates of the unknown
!        parameters are stored.
!     real pv(200)
!        the predicted value based on the current parameter estimates
!     real res(200)
!        the residuals from the fit.
!     real rsd
!        the value of the residual standard deviation at the solution.
!     real scale(50)
!        a value to indicate use of the default values of
!        the typical size of the unknown parameters.
!     real sdpv(200)
!        the standard deviation of the predicted value.
!     real sdres(200)
!        the standard deviations of the residuals.
!     real stopp
!        the stopping criterion for the test based on the maximum scaled
!        relative change in the elements of the model parameter vector
!     real stopss
!        the stopping criterion for the test based on the ratio of the
!        predicted decrease in the residual sum of squares (computed
!        by starpac) to the current residual sum of squares estimate.
!     real stp(50)
!        the rcstep size array.
!     real vcv(10,10)
!        the covariance matrix.
!     real y(200),ylog(200),yt(200)
!        the array of the dependent variable.
!

  double precision dstak(12)
  integer ierr
  integer ldstak
!
!  local scalars
  real &
     delta,rsd,stopp,stopss
  integer ifcst,ivaprx,ivcv,mit,nfac,npar,npare,nprt,ny
!
!  local arrays
  real &
     fcst(50,5),fcstsd(50),par(50),pv(200),res(200),scale(50), &
     sdpv(200),sdres(200),stp(50),vcv(10,10),y(200),ylog(200), &
     yt(200)
  integer &
     ifixed(50),mspec(4,50)
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr

  data    y(  1),   y(  2),   y(  3),   y(  4),   y(  5),   y(  6) &
      / 112.0e0, 118.0e0, 132.0e0, 129.0e0, 121.0e0, 135.0e0/
  data    y(  7),   y(  8),   y(  9),   y( 10),   y( 11),   y( 12) &
      / 148.0e0, 148.0e0, 136.0e0, 119.0e0, 104.0e0, 118.0e0/
  data    y( 13),   y( 14),   y( 15),   y( 16),   y( 17),   y( 18) &
      / 115.0e0, 126.0e0, 141.0e0, 135.0e0, 125.0e0, 149.0e0/
  data    y( 19),   y( 20),   y( 21),   y( 22),   y( 23),   y( 24) &
      / 170.0e0, 170.0e0, 158.0e0, 133.0e0, 114.0e0, 140.0e0/
  data    y( 25),   y( 26),   y( 27),   y( 28),   y( 29),   y( 30) &
      / 145.0e0, 150.0e0, 178.0e0, 163.0e0, 172.0e0, 178.0e0/
  data    y( 31),   y( 32),   y( 33),   y( 34),   y( 35),   y( 36) &
      / 199.0e0, 199.0e0, 184.0e0, 162.0e0, 146.0e0, 166.0e0/
  data    y( 37),   y( 38),   y( 39),   y( 40),   y( 41),   y( 42) &
      / 171.0e0, 180.0e0, 193.0e0, 181.0e0, 183.0e0, 218.0e0/
  data    y( 43),   y( 44),   y( 45),   y( 46),   y( 47),   y( 48) &
      / 230.0e0, 242.0e0, 209.0e0, 191.0e0, 172.0e0, 194.0e0/
  data    y( 49),   y( 50),   y( 51),   y( 52),   y( 53),   y( 54) &
      / 196.0e0, 196.0e0, 236.0e0, 235.0e0, 229.0e0, 243.0e0/
  data    y( 55),   y( 56),   y( 57),   y( 58),   y( 59),   y( 60) &
      / 264.0e0, 272.0e0, 237.0e0, 211.0e0, 180.0e0, 201.0e0/
  data    y( 61),   y( 62),   y( 63),   y( 64),   y( 65),   y( 66) &
      / 204.0e0, 188.0e0, 235.0e0, 227.0e0, 234.0e0, 264.0e0/
  data    y( 67),   y( 68),   y( 69),   y( 70),   y( 71),   y( 72) &
      / 302.0e0, 293.0e0, 259.0e0, 229.0e0, 203.0e0, 229.0e0/
  data    y( 73),   y( 74),   y( 75),   y( 76),   y( 77),   y( 78) &
      / 242.0e0, 233.0e0, 267.0e0, 269.0e0, 270.0e0, 315.0e0/
  data    y( 79),   y( 80),   y( 81),   y( 82),   y( 83),   y( 84) &
      / 364.0e0, 347.0e0, 312.0e0, 274.0e0, 237.0e0, 278.0e0/
  data    y( 85),   y( 86),   y( 87),   y( 88),   y( 89),   y( 90) &
      / 284.0e0, 277.0e0, 317.0e0, 313.0e0, 318.0e0, 374.0e0/
  data    y( 91),   y( 92),   y( 93),   y( 94),   y( 95),   y( 96) &
      / 413.0e0, 405.0e0, 355.0e0, 306.0e0, 271.0e0, 306.0e0/
  data    y( 97),   y( 98),   y( 99),   y(100),   y(101),   y(102) &
      / 315.0e0, 301.0e0, 356.0e0, 348.0e0, 355.0e0, 422.0e0/
  data    y(103),   y(104),   y(105),   y(106),   y(107),   y(108) &
      / 465.0e0, 467.0e0, 404.0e0, 347.0e0, 305.0e0, 336.0e0/
  data    y(109),   y(110),   y(111),   y(112),   y(113),   y(114) &
      / 340.0e0, 318.0e0, 362.0e0, 348.0e0, 363.0e0, 435.0e0/
  data    y(115),   y(116),   y(117),   y(118),   y(119),   y(120) &
      / 491.0e0, 505.0e0, 404.0e0, 359.0e0, 310.0e0, 337.0e0/
  data    y(121),   y(122),   y(123),   y(124),   y(125),   y(126) &
      / 360.0e0, 342.0e0, 406.0e0, 396.0e0, 420.0e0, 472.0e0/
  data    y(127),   y(128),   y(129),   y(130),   y(131),   y(132) &
      / 548.0e0, 559.0e0, 463.0e0, 407.0e0, 362.0e0, 405.0e0/
  data    y(133),   y(134),   y(135),   y(136),   y(137),   y(138) &
      / 417.0e0, 391.0e0, 419.0e0, 461.0e0, 472.0e0, 535.0e0/
  data    y(139),   y(140),   y(141),   y(142),   y(143),   y(144) &
      / 622.0e0, 606.0e0, 508.0e0, 461.0e0, 390.0e0, 432.0e0/
!
!  test against published results
!
  ny = 144
  ylog(1:144) = log(y(1:144))

  nfac = 2
  mspec(1,1) = 0
  mspec(2,1) = 1
  mspec(3,1) = 1
  mspec(4,1) = 1

  mspec(1,2) = 0
  mspec(2,2) = 1
  mspec(3,2) = 1
  mspec(4,2) = 12

  npar = 3
  par(1) = 0.0e0
  par(2) = 0.40e0
  par(3) = 0.60e0

  ifixed(1) = 1
  ifixed(2) = 0
  ifixed(3) = 0

  stopss = -1.0e0
  stopp = -1.0e0
  scale(1) = -1.0e0
  scale(2) = 1.0e-7
  scale(3) = 1.0e-7
  stp(1) = -1.0e0
  stp(2) = 1.0e-7
  stp(3) = 1.0e-7
  mit = 0
  nprt = -1
  delta = -1.0e0
  ivaprx = -1

  write ( *, 1000)
  call aimec ( ylog, ny, mspec, nfac, &
     par, npar, res,ldstak, ifixed,stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt )

  write ( *, 1005)
  par(1) = 0.0e0
  par(2) = 0.395e0
  par(3) = 0.615e0
  call aimfs ( ylog, ny, mspec, nfac, &
     par, npar, ldstak, ny/10+1, 1, ny, nprt, fcst, 50, fcstsd )

  scale(1) = 1.0e-7
  scale(2) = 1.0e-7
  scale(3) = 1.0e-7

  nfac = 2
  mspec(1,1) = 0
  mspec(2,1) = 1
  mspec(3,1) = 1
  mspec(4,1) = 1

  mspec(1,2) = 0
  mspec(2,2) = 0
  mspec(3,2) = 1
  mspec(4,2) = 12

  write ( *, 1000)
  call aimec ( ylog, ny, mspec, nfac, &
     par, npar, res,ldstak, ifixed,stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt )
  ny = 20
  write ( *, 1000)
  call aimec ( ylog, ny, mspec, nfac, &
     par, npar, res,ldstak, ifixed,stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt )

  nfac = 2
  mspec(1,1) = 0
  mspec(2,1) = 0
  mspec(3,1) = 1
  mspec(4,1) = 1

  mspec(1,2) = 0
  mspec(2,2) = 0
  mspec(3,2) = 1
  mspec(4,2) = 12

  ny = 144
  write ( *, 1000)
  call aimec ( ylog, ny, mspec, nfac, &
     par, npar, res,ldstak, ifixed,stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt )
!
!  example from page 212 of box and jenkins (1970)
!  add print statements to mdlts2 to check computations
!  at first call against those listed on page 214.
!
  write ( *, 1000)
  ny = 10
  yt(1) = 460.0e0
  yt(2) = 457.0e0
  yt(3) = 452.0e0
  yt(4) = 459.0e0
  yt(5) = 462.0e0
  yt(6) = 459.0e0
  yt(7) = 463.0e0
  yt(8) = 479.0e0
  yt(9) = 493.0e0
  yt(10) = 490.0e0

  nfac = 1
  mspec(1,1) = 0
  mspec(2,1) = 1
  mspec(3,1) = 1
  mspec(4,1) = 1

  npar = 2
  par(1) = 0.0e0
  par(2) = 0.5e0

  ifixed(1) = 1
  ifixed(2) = 0

  call aimec ( yt, ny, mspec, nfac, &
     par, npar, res,ldstak, ifixed,stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt )
!
!  example from page 216 of box and jenkins (1970)
!  add print statements to mdlts2 to check computations
!  at first call against those listed on page 218.
!
  write ( *, 1000)
  ny = 12
  yt(1) = 2.0e0
  yt(2) = 0.8e0
  yt(3) = -0.3e0
  yt(4) = -0.3e0
  yt(5) = -1.9e0
  yt(6) = 0.3e0
  yt(7) = 3.2e0
  yt(8) = 1.6e0
  yt(9) = -0.7e0
  yt(10) = 3.0e0
  yt(11) = 4.3e0
  yt(12) = 1.1e0

  nfac = 1
  mspec(1,1) = 1
  mspec(2,1) = 0
  mspec(3,1) = 1
  mspec(4,1) = 1

  npar = 3
  par(1) = 0.3e0
  par(2) = 0.0e0
  par(3) = 0.7e0

  ifixed(1) = 0
  ifixed(2) = 1
  ifixed(3) = 0

  call aimec ( yt, ny, mspec, nfac, &
     par, npar, res, ldstak, ifixed, stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt )
!
!  test error messages
!
  write ( *, 1010)
  ny = 0
  nfac = 0
  call aime ( yt, ny, mspec, nfac, &
     par, npar, res, ldstak )
  call aimec ( yt, ny, mspec, nfac, &
     par, npar, res, ldstak, ifixed, stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt )
  call aimes ( yt, ny, mspec, nfac, &
     par, npar, res, ldstak, ifixed, stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt, npare, rsd, pv, sdpv, sdres, vcv, &
     ivcv )
  call aimf ( y, ny, mspec, nfac, par, npar, ldstak )
  call aimfs ( y, ny, mspec, nfac, &
     par, npar, ldstak, ny/10+1, 1, ny, nprt, fcst, 50, fcstsd )

  ny = 144
  nfac = 2
  mspec(1,1) = -1
  call aime ( yt, ny, mspec, nfac, &
     par, npar, res, ldstak )
  call aimec ( yt, ny, mspec, nfac, &
     par, npar, res, ldstak, ifixed, stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt )
  call aimes ( yt, ny, mspec, nfac, &
     par, npar, res, ldstak, ifixed, stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt, npare, rsd, pv, sdpv, sdres, vcv, &
     ivcv )
  call aimf ( y, ny, mspec, nfac, par, npar, ldstak )
  call aimfs ( y, ny, mspec, nfac, &
     par, npar, ldstak, ny/10+1, 1, ny, nprt, fcst, 50, fcstsd )
  ny = 144
  nfac = 2
  mspec(1,1) = 0
  npar = 1
  call aime ( yt, ny, mspec, nfac, &
     par, npar, res, ldstak )
  call aimec ( yt, ny, mspec, nfac, &
     par, npar, res, ldstak, ifixed, stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt )
  call aimes ( yt, ny, mspec, nfac, &
     par, npar, res, ldstak, ifixed, stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt, npare, rsd, pv, sdpv, sdres, vcv, &
     ivcv )
  call aimf ( y, ny, mspec, nfac, par, npar, ldstak )
  call aimfs ( y, ny, mspec, nfac, &
     par, npar, ldstak, ny/10+1, 1, ny, nprt, fcst, 50, fcstsd )
  ny = 144
  nfac = 2
  mspec(1,1) = 0
  npar = 3
  ifixed(1:npar) = 1
  ivcv = 0
  ifcst = 0
  call aimec ( yt, ny, mspec, nfac, &
     par, npar, res, ldstak, ifixed, stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt )
  call aimes ( yt, ny, mspec, nfac, &
     par, npar, res, ldstak, ifixed, stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt, npare, rsd, pv, sdpv, sdres, vcv, &
     ivcv )
  call aimfs ( y, ny, mspec, nfac, &
     par, npar, ldstak, ny/10+1, 1, ny, nprt, fcst, ifcst, fcstsd )
  ifixed(1:npar) = 1
  ivcv = 0
  stp(2) = -1.0e0
  scale(2) = -1.0e0
  call aimec ( yt, ny, mspec, nfac, &
     par, npar, res, ldstak, ifixed, stp, mit, stopss, stopp, &
     scale, delta, ivaprx, nprt )
  return

 1000 format ('1test of arima estimation routines')
 1005 format ('1test of arima forecasting routines')
 1010 format ('1test of error checking facilities')
end
