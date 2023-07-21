subroutine xaimd ( ldstak ,typ)
use M_starpac, only : aime,aimec,aimes,aimf,aimfs,aimx1,fitxsp,iprint
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XAIMD demonstrates the user callable routines in the ARIMA family.
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
!  Parameters:
!
!     real delta
!        the maximum change allowed in the model parameters at the
!        first iteration.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     real fcst(200, 10)
!        the forecasts.
!     real fcstsd(200, 10)
!        the standard deviations of the forecasts.
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr .ge. 1, errors were detected.
!     integer ifcst
!        the first dimension of the array fcst.
!     integer ifcsto(10)
!        the indices of the origins for the forecasts.
!     integer ifixed(10)
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
!     integer mspec(4,10)
!        the array containing the values of p, d, q, and s for each
!        factor.
!     integer n
!        the number of observations.
!     integer nfac
!        the number of factors in the model
!     integer nfcst
!        the number of forecasts.
!     integer nfcsto
!        the number of the origins.
!     integer npar
!        the number of unknown parameters in the model.
!     integer npare
!        the number of parameters estimated by the routine.
!     integer nprt
!        the parameter used to indicate how much printed output is
!        to be provided.
!     integer ntest
!        the number of the current test.
!     real par(10)
!        the array in which the current estimates of the unknown
!        parameters are stored.
!     real pv(200)
!        the predicted value based on the current parameter estimates
!     real res(200)
!        the residuals from the fit.
!     real rsd
!        the value of the residual standard deviation at the solution.
!     real scale(10)
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
!     real stp(10)
!        the rcstep size array.
!     real vcv(10,10)
!        the covariance matrix.
!     real y(200)
!        the array of the dependent variable.
!

  real delta
  double precision dstak(12)
  integer ierr
  integer ldstak
  integer mspec(4,10)
  integer n
  integer nfac
  integer npar
  real par(10)
  real res(200)
  real rsd
  real stopp
  real stopss
  real y(200)

  integer &
     ifcst,ivaprx,ivcv,mit,mxfac,mxfc,mxfco,mxn,mxpar, &
     nfcst,nfcsto,npare,nprt,ntest
!
!  local arrays
  real &
     fcst(200,10),fcstsd(200,10),pv(200), &
     scale(10),sdpv(200),sdres(200),stp(10),vcv(10,10)
  integer &
     ifcsto(10),ifixed(10)
!
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr
!
!  define constants
!
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

  y(1:144) = log(y(1:144))
!
!  set dimensions
!
  mxn = 200
  mxpar = 10
  mxfc = 200
  mxfco = 10
  mxfac = 10

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'XAIMD'
  write ( *, '(a)' ) '  Demonstrate codes in the ARIMA family.'

  ntest = 0
!
!  Test on normal statement.
!
  ntest = ntest + 1

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  ARIMA test number ', ntest
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) '  Normal problem'
  write ( *, '(a)' ) '  Test of AIM'
!
!  Set the starting parameters for AIMX.
!
  call aimx1 ( mxn, mxpar, mxfc, mxfco, mxfac, &
     1, n, mspec, nfac, par, npar, res, &
     ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, nprt, &
     npare, rsd, pv, sdpv, sdres, vcv, ivcv, &
     nfcst, nfcsto, ifcsto, fcst, ifcst, fcstsd )

  call aime ( y, n, mspec, nfac, par, npar, res, ldstak )

  write ( *,1120) ierr

  call fitxsp ( par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, &
     n, npare, rsd )

  ntest = ntest + 1

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  ARIMA test number ', ntest
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Normal problem'
  write ( *, '(a)' ) '  Test of AIMC'
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt

  call aimx1 ( mxn, mxpar, mxfc, mxfco, mxfac, &
     1, n, mspec, nfac, par, npar, res, &
     ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, nprt, &
     npare, rsd, pv, sdpv, sdres, vcv, ivcv, &
     nfcst, nfcsto, ifcsto, fcst, ifcst, fcstsd )

  call aimec ( y, n, mspec, nfac, par, npar, res, ldstak, &
     ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, nprt )

  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr

  call fitxsp ( par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, &
     n, npare, rsd )

  ntest = ntest + 1

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  ARIMA test number ', ntest
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) '  Normal problem'
  write ( *, '(a)' ) '  Test of AIMS'
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt

  call aimx1 ( mxn, mxpar, mxfc, mxfco, mxfac, &
     1, n, mspec, nfac, par, npar, res, &
     ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, nprt, &
     npare, rsd, pv, sdpv, sdres, vcv, ivcv, &
     nfcst, nfcsto, ifcsto, fcst, ifcst, fcstsd)

  call aimes ( y, n, mspec, nfac, par, npar, res, ldstak, &
     ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, nprt, &
     npare, rsd, pv, sdpv, sdres, vcv, ivcv )
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp ( par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, &
     n, npare, rsd )

  ntest = ntest + 1

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  ARIMA test number ', ntest
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) '  Normal problem'
  write ( *, '(a)' ) '  Test of AIMF'

  call aimx1 ( mxn, mxpar, mxfc, mxfco, mxfac, &
     1, n, mspec, nfac, par, npar, res, &
     ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, nprt, &
     npare, rsd, pv, sdpv, sdres, vcv, ivcv, &
     nfcst, nfcsto, ifcsto, fcst, ifcst, fcstsd )

  call aimf ( y, n, mspec, nfac, par, npar, ldstak )
  write ( *,1120) ierr

  ntest = ntest + 1

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  ARIMA test number ', ntest
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) '  Normal problem'
  write ( *, '(a)' ) '  Test of AIMFS'
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt

  call aimx1 ( mxn, mxpar, mxfc, mxfco, mxfac, &
     1, n, mspec, nfac, par, npar, res, &
     ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, nprt, &
     npare, rsd, pv, sdpv, sdres, vcv, ivcv, &
     nfcst, nfcsto, ifcsto, fcst, ifcst, fcstsd)

  call aimfs ( y, n, mspec, nfac, par, npar, ldstak, &
     nfcst, nfcsto, ifcsto, nprt, fcst, ifcst, fcstsd )
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr

  call fitxsp ( par, fcst(1,1), fcst(1,2), fcst(1,3), fcstsd, vcv, &
    n, npar, ivcv, n, npare, rsd )

  return

!1010 format ('test of aimc'  )
!1020 format ('test of aims'  )
!1030 format ('test of aimf' )
!1040 format ('test of aimfs' )

 1120 format (/' ***** returned results *****', 5x, '(-1 indicates ', &
     'value not changed by called subroutine)'//' ierr is ', i3)
 1340 format (//24h input   -  ifixed(1) = , i6, 9x, ', stp(1) = ', &
     g15.8, ',    mit = ',i5, ', stopss = ', g15.8, 10h, stopp = , &
     g15.8/13x, 'scale(1) = ', g15.8, ',  delta = ', g15.8, &
     ', ivaprx = ', i5, ',   nprt = ', i5//)
 1350 format (//24h output  -  ifixed(1) = , i6, 9x, ', stp(1) = ', &
     g15.8, ',    mit = ',i5, ', stopss = ', g15.8, 10h, stopp = , &
     g15.8/13x, 'scale(1) = ', g15.8, ',  delta = ', g15.8, &
     ', ivaprx = ', i5, ',   nprt = ', i5//)
end
