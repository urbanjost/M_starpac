subroutine xnlst ( ldstak ,typ)
use M_starpac, only : rmdcon, &
  drv1a,drv2,drv3,fitxsp,iprint,ldscmp,mdl1,mdl2,mdl3,nlsc, &
     nlsdc,nlsds,nlss,nlswc,nlswdc,nlswds,nlsws,nlsx1,nlsx2, &
     setrv
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ
!*****************************************************************************80
!
!! XNLST demonstrates the nonlinear least squares family.
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
!     external drv1a, drv2, drv3
!        the name of the user supplied subroutine which computes the
!        derivative (jacobian) matrix of the model.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer i
!        an index value.
!     integer idrvck
!        the variable used to indicate whether the derivatives are
!        to be checked (idrvck = 1) or not (idrvck = 0).
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr .ge. 1, errors were detected.
!     integer ifixed(10)
!        the indicator values used to designate whether the
!        parameters are to be optimized or are to be held fixed.  if
!        ifixed(i).ne.0, then par(i) will be optimized.  if
!        ifixed(i)==0, then par(i) will be held fixed.
!        ifixed(i).lt.0, then all par(i),i=1,npar, will be optimized..
!     integer ivctst(9)
!        variance-covariance code test values.
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
!     integer ixm1, ixm2, ixm3
!        the first dimension of the independent variable array.
!     integer ldsa1, ldsn1b
!        the minimum length allowed for the array dstak
!        for the routines with analytic derivatives and
!        numerical derivatives, respectively.
!     integer ldstak
!        the length of the array dstak.
!     integer m1, m2, m3
!        the number of independent variables.
!     integer mit
!        the maximum number of iterations allowed.
!     external mdl1, mdl2, mdl3
!        the name of the user supplied subroutine which computes the
!        predicted values based on the current parameter estimate.
!     integer n1, n2, n3
!        the number of observations.
!     integer npare
!        the number of parameters estimated by the routine.
!     integer npar1, npar2, npar3
!        the number of unknown parameters in the model.
!     integer nnzw
!        the number of non zero weights.
!     integer nprt
!        the parameter used to indicate how much printed output is
!        to be provided.
!     integer ntest
!        the number of the current test.
!     real par1(10), par2(10), par3(10)
!        the array in which the current estimates of the unknown
!        parameters are stored.
!     real pv(100)
!        the starting location in rstak/dstak of
!        the predicted value based on the current parameter estimates
!     real res(100)
!        the residuals from the fit.
!     real rsd
!        the value of the residual standard deviation at the solution.
!     real scale(10)
!        a value to indicate use of the default values of
!        the typical size of the unknown parameters.
!     real sdpv(100)
!        the starting location in rstak/dstak of
!        the standard deviation of the predicted value.
!     real sdres(100)
!        the starting location in rstak/dstak of the
!        the standard deviations of the residuals.
!     real stop(8)
!        stopping criteria test variable.
!     real stopp
!        the stopping criterion for the test based on the maximum scaled
!        relative change in the elements of the model parameter vector
!     real stopss
!        the stopping criterion for the test based on the ratio of the
!        predicted decrease in the residual sum of squares (computed
!        by starpac) to the current residual sum of squares estimate.
!     real stp(10)
!        the rcstep size array.
!     real vcv(6,6)
!        the covariance matrix.
!     real wt(100)
!        the user supplied weights.
!     real xm1(10,2), xm2(10,3), xm3(101,5)
!        the array in which one row of the independent variable array
!        is stored.
!     real y1(10), y2(10), y3(100)
!        the array of the dependent variable.
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
     delta,rsd,stopp,stopss
  integer &
     i,idrvck,ivaprx,ivcv,ixm1,ixm2,ixm3,ldsa1,ldsmin, &
     ldsn1b,m1,m2,m3,mit,n1,n2,n3,nnzw,npar1,npar2,npar3,npare, &
     nprt,ntest
!
!  local arrays
  real &
     par1(10),par2(10),par3(10),pv(100),res(100),scale(10), &
     sdpv(100),sdres(100),stop(8),stp(10),vcv(6,6),wt(100), &
     xm1(10,2),xm2(10,3),xm3(101,5),y1(10),y2(10),y3(100)
  integer &
     ifixed(10),ivctst(9)
!

!  common blocks
  common /cstak/dstak
  common /errchk/ierr
!
!     define constants
!
  data y1(1), y1(2), y1(3), y1(4), y1(5), y1(6) &
     /2.138e0, 3.421e0, 3.597e0, 4.340e0, 4.882e0, 5.660e0/

  data xm1(1,1), xm1(2,1), xm1(3,1), xm1(4,1), xm1(5,1), xm1(6,1) &
     /1.309e0, 1.471e0, 1.490e0, 1.565e0, 1.611e0, 1.680e0/

  data n2 /10/, m2 /3/, ixm2 /10/, npar2 /3/

  data n3 /50/, m3 /5/, ixm3 /101/, npar3 /5/

  call nlsx2(n1, m1, ixm1, npar1, ifixed, stp, idrvck, mit, stopss, &
     stopp, scale, delta, ivaprx, nprt, ivcv)

  call setrv(wt, n3, 1.0e0)

  call ldscmp(6, 0, 60+2*npar1, 0, 0, 0, &
     's', 94+n1*(3+npar1)+npar1*(3*npar1+35)/2, ldsa1)
  call ldscmp(14, 0, 60+2*npar1, 0, 0, 0, &
     's', 94+n1*(3+npar1)+npar1*(3*npar1+37)/2, ldsn1b)
  call ldscmp(14, 0, max(2*(n3+npar3),60+2*npar3), 0, 0, 0, &
     's', max(10*n3,94+n3*(3+npar3)+npar3*(3*npar3+37)/2), ldsmin)

  if (ldsmin.le.ldstak) go to 5

  write ( *, 1000) ldsmin
  return

    5 continue

  do i=1,n2
     y2(i) = 0.0e0
     xm2(i,1) = i
     xm2(i,2) = i + 0.125e0
     xm2(i,3) = i + 0.25e0
  end do

  do i=1,n3
     xm3(i,1) = 1.0e0
     xm3(i,2) = i
     xm3(i,3) = xm3(i,2)*xm3(i,2)
     xm3(i,4) = xm3(i,3)*xm3(i,2)
     xm3(i,5) = xm3(i,4)*xm3(i,2)
     y3(i) = xm3(i,1) + xm3(i,2) + xm3(i,3) + xm3(i,4) + xm3(i,5)
  end do

  ntest = 0
!
!  test checking of control criteria
!
  write ( *,1240)

  stop(1) = rmdcon(3,typ=0.0_wp)
  stop(2) = 0.1e0
  stop(3) = 0.9e0*rmdcon(3,typ=0.0_wp)
  stop(4) = 0.11e0
  stop(5) = 0.0e0
  stop(6) = 1.0e0
  stop(7) = -1.0e0
  stop(8) = 1.1e0

  nprt = 11000
  mit = 0
  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1250) mit
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  write ( *,1100)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, scale(1), &
     delta, ivaprx, nprt
  call nlswdc(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, &
     res, ldstak, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, scale(1), &
     delta, ivaprx, nprt
  write ( *,1120) ierr

  mit = 1
  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1250) mit
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  write ( *,1100)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, scale(1), &
     delta, ivaprx, nprt
  call nlswdc(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, &
     res, ldstak, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, scale(1), &
     delta, ivaprx, nprt
  write ( *,1120) ierr

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1250) mit
  call nlsx1(2, par2, npar2, pv, sdpv, res, sdres, vcv, n2, ivcv, &
     nnzw, npare, rsd)
  write ( *,1050)
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsws(y2, wt, xm2, n2, m2, ixm2, mdl2, par2, npar2, res, &
     ldstak, ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, &
     nprt, nnzw, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par2, pv, sdpv, res, sdres, vcv, n2, npar2, ivcv, &
     nnzw, npare, rsd)

  call nlsx2(n1, m1, ixm1, npar1, ifixed, stp, idrvck, mit, stopss, &
     stopp, scale, delta, ivaprx, nprt, ivcv)

  do i=1,4
     ntest = ntest + 1
     write ( *,1330) ntest
     write ( *,1260) stop(i)
     call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, &
        ivcv, nnzw, npare, rsd)
     write ( *,1100)
     write ( *,1360) ifixed(1), idrvck, mit, stop(i), stopp, &
        scale(1), delta, ivaprx, nprt
     call nlswdc(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, &
        npar1, res, ldstak, ifixed, idrvck, mit, stop(i), stopp, &
        scale, delta, ivaprx, nprt)
     write ( *,1370) ifixed(1), idrvck, mit, stop(i), stopp, &
        scale(1), delta, ivaprx, nprt
     write ( *,1120) ierr
  end do

  do i=5,8
     ntest = ntest + 1
     write ( *,1330) ntest
     write ( *,1270) stop(i)
     call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, &
        ivcv, nnzw, npare, rsd)
     write ( *,1100)
     write ( *,1360) ifixed(1), idrvck, mit, stopss, stop(i), &
        scale(1), delta, ivaprx, nprt
     call nlswdc(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, &
        npar1, res, ldstak, ifixed, idrvck, mit, stopss, stop(i), &
        scale, delta, ivaprx, nprt)
     write ( *,1370) ifixed(1), idrvck, mit, stopss, stop(i), &
        scale(1), delta, ivaprx, nprt
     write ( *,1120) ierr
  end do

  nprt = 100000

  do i=1,6
     ntest = ntest + 1
     write ( *,1330) ntest
     write ( *,1280) nprt
     call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, &
        ivcv, nnzw, npare, rsd)
     write ( *,1110)
     write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
        scale(1), delta, ivaprx, nprt
     call nlswds(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, &
        npar1, res, ldsa1, ifixed, idrvck, mit, stopss, stopp, &
        scale, delta, ivaprx, nprt, nnzw, npare, rsd, pv, sdpv, &
        sdres, vcv, ivcv)
     write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
        scale(1), delta, ivaprx, nprt
     write ( *,1120) ierr
     call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
        nnzw, npare, rsd)
     nprt = nprt/10
  end do

  nprt = 11000

  do i=1,2
     ntest = ntest + 1
     write ( *,1330) ntest
     write ( *,1280) nprt
     call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, &
        ivcv, nnzw, npare, rsd)
     write ( *,1110)
     write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
        scale(1), delta, ivaprx, nprt
     call nlswds(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, &
        npar1, res, ldsa1, ifixed, idrvck, mit, stopss, stopp, &
        scale, delta, ivaprx, nprt, nnzw, npare, rsd, pv, sdpv, &
        sdres, vcv, ivcv)
     write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
        scale(1), delta, ivaprx, nprt
     write ( *,1120) ierr
     call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
        nnzw, npare, rsd)
     nprt = 11001
  end do

  nprt = 0
  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1280) nprt

  write ( *,1010)
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlsc(y1, xm1, n1, m1, ixm1, mdl1, par1, npar1, res, ldsn1b, &
     ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, nprt)
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1280) nprt
  write ( *,1020)
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlss(y1, xm1, n1, m1, ixm1, mdl1, par1, npar1, res, ldsn1b, &
     ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, nprt, &
     npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1280) nprt
  write ( *,1040)
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlswc(y1, wt, xm1, n1, m1, ixm1, mdl1, par1, npar1, res, &
     ldsn1b, ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, &
     nprt)
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1280) nprt
  write ( *,1050)
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlsws(y1, wt, xm1, n1, m1, ixm1, mdl1, par1, npar1, res, &
     ldsn1b, ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, &
     nprt, nnzw, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1280) nprt
  write ( *,1070)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlsdc(y1, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, res, &
     ldsa1, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1280) nprt
  write ( *,1080)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlsds(y1, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, res, &
     ldsa1, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1280) nprt
  write ( *,1100)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlswdc(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, &
     res, ldsa1, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1280) nprt
  write ( *,1110)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlswds(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, &
     res, ldsa1, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt, nnzw, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  nprt = -1

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1280) nprt
  write ( *,1110)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlswds(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, &
     res, ldsa1, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt, nnzw, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)
!
!  test parameter handling**
!
  write ( *,1190)
!
!  all zero
!
  call nlsx2(n1, m1, ixm1, npar1, ifixed, stp, idrvck, mit, stopss, &
     stopp, scale, delta, ivaprx, nprt, ivcv)
  stp(1) = -1.0e0
  scale(1) = -1.0e0
  delta = -1.0e0
  nprt = 11000

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1200)
  call nlsx1(4, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  write ( *,1010)
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsc(y1, xm1, n1, m1, ixm1, mdl1, par1, npar1, res, ldstak, &
     ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, nprt)
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  call nlsx1(4, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1200)
  write ( *,1100)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlswdc(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, &
     res, ldstak, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)
!
!  test with constant y
!
!  constant y=0
!
  nprt = 21222
  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1210)
  call nlsx1(2, par2, npar2, pv, sdpv, res, sdres, vcv, n2, ivcv, &
     nnzw, npare, rsd)
  write ( *,1050)
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsws(y2, wt, xm2, n2, m2, ixm2, mdl2, par2, npar2, res, &
     ldstak, ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, &
     nprt, nnzw, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par2, pv, sdpv, res, sdres, vcv, n2, npar2, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1210)
  write ( *,1110)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(2, par2, npar2, pv, sdpv, res, sdres, vcv, n2, ivcv, &
     nnzw, npare, rsd)
  call nlswds(y2, wt, xm2, n2, m2, ixm2, mdl2, drv2, par2, npar2, &
     res, ldstak, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt, nnzw, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par2, pv, sdpv, res, sdres, vcv, n2, npar2, ivcv, &
     nnzw, npare, rsd)
!
!  test with linear model
!
  nprt = 11212
  ivaprx = 1

  ifixed(1:npar3) = 0
  ifixed(1) = 1

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1220)
  write ( *,1010)
  write ( *,1340) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(3, par3, npar3, pv, sdpv, res, sdres, vcv, n3, ivcv, &
     nnzw, npare, rsd)
  call nlsc(y3, xm3, n3, m3, ixm3, mdl3, par3, npar3, res, ldstak, &
     ifixed, stp, mit, stopss, stopp, scale, delta, ivaprx, nprt)
  write ( *,1350) ifixed(1), stp(1), mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1220)
  nprt = 11111
  write ( *,1070)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(3, par3, npar3, pv, sdpv, res, sdres, vcv, n3, ivcv, &
     nnzw, npare, rsd)
  call nlsdc(y3, xm3, n3, m3, ixm3, mdl3, drv3, par3, npar3, res, &
     ldstak, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
!
!  test xm
!
!  first column zero
!
  call setrv(y2, n2, 2.0e0)
  call setrv(xm2(1,1), n2, 0.0e0)

  nprt = 11000

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1310)
  call nlsx1(2, par2, npar2, pv, sdpv, res, sdres, vcv, n2, ivcv, &
     nnzw, npare, rsd)
  write ( *,1110)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlswds(y2, wt, xm2, n2, m2, ixm2, mdl2, drv2, par2, npar2, &
     res, ldstak, ifixed, idrvck, mit, stopss, stopp, scale, &
     delta, ivaprx, nprt, nnzw, npare, rsd, pv, sdpv, sdres, vcv, &
     ivcv)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par2, pv, sdpv, res, sdres, vcv, n2, npar2, ivcv, &
     nnzw, npare, rsd)
!
!  2 columns zero
!
  call setrv(xm2(1,2), n2, 0.0e0)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1300)
  write ( *,1110)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(2, par2, npar2, pv, sdpv, res, sdres, vcv, n2, ivcv, &
     nnzw, npare, rsd)
  call nlswds(y2, wt, xm2, n2, m2, ixm2, mdl2, drv2, par2, npar2, &
     res, ldstak, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt, nnzw, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par2, pv, sdpv, res, sdres, vcv, n2, npar2, ivcv, &
     nnzw, npare, rsd)
!
!  test variance-covariance matrix computations
!
  ivctst(1) = -1
  ivctst(2) = 0
  ivctst(3) = 1
  ivctst(4) = 2
  ivctst(5) = 3
  ivctst(6) = 4
  ivctst(7) = 5
  ivctst(8) = 6
  ivctst(9) = 7
  nprt = 2
  do i=1,9
     ntest = ntest + 1
     write ( *,1330) ntest
     write ( *,1380)
     call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, &
        ivcv, nnzw, npare, rsd)
     write ( *,1110)
     write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
        scale(1), delta, ivctst(i), nprt
     call nlswds(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, &
        npar1, res, ldstak, ifixed, idrvck, mit, stopss, stopp, &
        scale, delta, ivctst(i), nprt, nnzw, npare, rsd, pv, sdpv, &
        sdres, vcv, ivcv)
     write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
        scale(1), delta, ivctst(i), nprt
     write ( *,1120) ierr
     call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
        nnzw, npare, rsd)
  end do
!
!  test with 2 zero weights
!
  nprt = 22222
  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1290)
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  wt(3) = 0.0e0
  wt(5) = 0.0e0
  write ( *,1110)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlswds(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, &
     res, ldstak, ifixed, idrvck, mit, stopss, stopp, scale, delta, &
     ivaprx, nprt, nnzw, npare, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  return

 1000 format ('ldstak must be greater than or equal to ', i6)
 1010 format ('test of nlsc'  )
 1020 format ('test of nlss'  )
 1040 format ('test of nlswc' )
 1050 format ('test of nlsws' )
 1070 format ('test of nlsdc' )
 1080 format ('test of nlsds' )
 1100 format ('test of nlswdc')
 1110 format ('test of nlswds')
 1120 format (/'returned results (-1 indicates ', &
     'value not changed by called subroutine)'//'ierr is ', i3)
 1190 format ('test parameter handling')
 1200 format ('all parameters zero')
 1210 format ('test with constant zero y')
 1220 format ('test linear model')
 1240 format ('test control criteria')
 1250 format ('maximum number of iterations = ', i5)
 1260 format (12h --stopss = , g14.8)
 1270 format (11h --stopp = , g14.8)
 1280 format (10h --nprt = , i6)
 1290 format (29h **test with 2 zero weights**)
 1300 format (19h **2 columns zero**)
 1310 format (18h **1 column zero**)
 1330 format (54h1nonlinear least squares estimation subroutine test nu, &
     4hmber, i5/)
 1340 format (24h input   -  ifixed(1) = , i6, 9x, 11h, stp(1) = , &
     g15.8, 11h,    mit = ,i5, 11h, stopss = , g15.8, 10h, stopp = , &
     g15.8/13x, 11hscale(1) = , g15.8, 11h,  delta = , g15.8, &
     11h, ivaprx = , i5, 11h,   nprt = , i6)
 1350 format (//24h output  -  ifixed(1) = , i6, 9x, 11h, stp(1) = , &
     g15.8, 11h,    mit = ,i5, 11h, stopss = , g15.8, 10h, stopp = , &
     g15.8/13x, 11hscale(1) = , g15.8, 11h,  delta = , g15.8, &
     11h, ivaprx = , i5, 11h,   nprt = , i6//)
 1360 format (24h input   -  ifixed(1) = , i6, 9x, 11h, idrvck = , &
     i5, 10x, 11h,    mit = ,i5, 11h, stopss = , g15.8, &
     10h, stopp = , g15.8/ &
     13x, 11hscale(1) = , g15.8, 11h,  delta = , g15.8, &
     11h, ivaprx = , i5, 11h,   nprt = , i6)
 1370 format (//24h output  -  ifixed(1) = , i6, 9x, 11h, idrvck = , &
     i5, 10x, 11h,    mit = ,i5, 11h, stopss = , g15.8, &
     10h, stopp = , g15.8/ &
     13x, 11hscale(1) = , g15.8, 11h,  delta = , g15.8, &
     11h, ivaprx = , i5, 11h,   nprt = , i6//)
 1380 format (54h test handling of variance-covariance computation code, &
     's')
end
