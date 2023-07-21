subroutine xnlsd ( ldstak ,typ)
use M_starpac, only : drv1a,fitxsp,iprint,ldscmp,mdl1,nl2x,nls,nlsc,nlsd,nlsdc, &
        & nlsds,nlss,nlsw,nlswc,nlswd,nlswdc,nlswds,nlsws,nlsx1, nlsx2,setrv
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XNLSD demonstrates the nonlinear least squares family.
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
!     external drv1a
!        the name of the user supplied subroutine which computes the
!        derivative (jacobian) matrix of the model.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
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
!     integer ixm1
!        the first dimension of the independent variable array.
!     integer ldsa1, ldsn1a, ldsn1b
!        the minimum length allowed for the array dstak
!        for the routines with analytic derivatives and
!        numerical derivatives, respectively.
!     integer ldsmin
!        the minimum length of the array dstak allowed.
!     integer ldstak
!        the length of the array dstak.
!     integer m1
!        the number of independent variables.
!     integer mit
!        the maximum number of iterations allowed.
!     external mdl1
!        the name of the user supplied subroutine which computes the
!        predicted values based on the current parameter estimate.
!     integer n1
!        the number of observations.
!     integer npare
!        the number of parameters estimated by the routine.
!     integer npar1
!        the number of unknown parameters in the model.
!     integer nnzw
!        the number of non zero weights.
!     integer nprt
!        the parameter used to indicate how much printed output is
!        to be provided.
!     integer ntest
!        the number of the current test.
!     real par1(10)
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
!     real xm1(10,2)
!        the array in which one row of the independent variable array
!        is stored.
!     real y1(10)
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
     idrvck,ivaprx,ivcv,ixm1,ldsa1,ldsmin,ldsn1a,ldsn1b, &
     m1,mit,n1,nnzw,npar1,npare,nprt,ntest
!
!  local arrays
  real &
     par1(10),pv(100),res(100),scale(10),sdpv(100),sdres(100), &
     stp(10),vcv(6,6),wt(100),xm1(10,2),y1(10)
  integer &
     ifixed(10)

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

  call nlsx2(n1, m1, ixm1, npar1, ifixed, stp, idrvck, mit, stopss, &
     stopp, scale, delta, ivaprx, nprt, ivcv)

  call setrv(wt, n1, 1.0e0)

  call ldscmp(6, 0, 60+2*npar1, 0, 0, 0, &
     's', 94+n1*(3+npar1)+npar1*(3*npar1+35)/2, ldsa1)
  call ldscmp(14, 0, max(2*(n1+npar1),60+2*npar1), 0, 0, 0, &
     's', max(10*n1,94+n1*(3+npar1)+npar1*(3*npar1+37)/2), ldsn1a)
  call ldscmp(14, 0, 60+2*npar1, 0, 0, 0, &
     's', 94+n1*(3+npar1)+npar1*(3*npar1+37)/2, ldsn1b)

  ldsmin = max(ldsa1, ldsn1a, ldsn1b)

  if (ldsmin.le.ldstak) go to 5

  write ( *, 1140) ldsmin
  return

    5 continue

  ntest = 0
!
!  test on normal statement
!
  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1130)
  write ( *,1000)
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nls(y1, xm1, n1, m1, ixm1, mdl1, par1, npar1, res, ldsn1a)
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1130)
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
  write ( *,1130)
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
  write ( *,1130)
  write ( *,1030)
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlsw(y1, wt, xm1, n1, m1, ixm1, mdl1, par1, npar1, res, &
     ldsn1a)
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1130)
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
  write ( *,1130)
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
  write ( *,1130)
  write ( *,1060)
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlsd(y1, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, res, &
     ldsa1)
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1130)
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
  write ( *,1130)
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
  write ( *,1130)
  write ( *,1090)
  write ( *,1360) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  call nlsx1(1, par1, npar1, pv, sdpv, res, sdres, vcv, n1, ivcv, &
     nnzw, npare, rsd)
  call nlswd(y1, wt, xm1, n1, m1, ixm1, mdl1, drv1a, par1, npar1, &
     res, ldsa1)
  write ( *,1370) ifixed(1), idrvck, mit, stopss, stopp, &
     scale(1), delta, ivaprx, nprt
  write ( *,1120) ierr
  call fitxsp(par1, pv, sdpv, res, sdres, vcv, n1, npar1, ivcv, &
     nnzw, npare, rsd)

  ntest = ntest + 1
  write ( *,1330) ntest
  write ( *,1130)
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
  write ( *,1130)
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
!     test nl2sol and nl2sno directly
!
  write ( *,1320)
  call nl2x(0.0_wp)

  return

 1000 format ('test of nls'  )
 1010 format ('test of nlsc'  )
 1020 format ('test of nlss'  )
 1030 format ('test of nlsw' )
 1040 format ('test of nlswc' )
 1050 format ('test of nlsws' )
 1060 format ('test of nlsd' )
 1070 format ('test of nlsdc' )
 1080 format ('test of nlsds' )
 1090 format ('test of nlswd')
 1100 format ('test of nlswdc')
 1110 format ('test of nlswds')
 1120 format (/'returned results (-1 indicates ', &
     'value not changed by called subroutine)'//' ierr is ', i3)
 1130 format ('normal problem')
 1140 format ('ldstak must be greater than or equal to ', i6)
 1320 format ('test of nl2sol and nl2sno called directly')
 1330 format ('nonlinear least squares estimation subroutine test number', &
     i5/)
 1340 format (' input   -  ifixed(1) = ', i6, 9x, ', stp(1) = ', &
     g15.8, ',    mit = ',i5, ', stopss = ', g15.8, ', stopp = ', &
     g15.8/13x, 'scale(1) = ', g15.8, ',  delta = ', g15.8, &
     ', ivaprx = ', i5, ',   nprt = ', i5)
 1350 format (//24h output  -  ifixed(1) = , i6, 9x, 11h, stp(1) = , &
     g15.8, 11h,    mit = ,i5, 11h, stopss = , g15.8, 10h, stopp = , &
     g15.8/13x, 11hscale(1) = , g15.8, 11h,  delta = , g15.8, &
     11h, ivaprx = , i5, 11h,   nprt = , i5//)
 1360 format (24h input   -  ifixed(1) = , i6, 9x, 11h, idrvck = , &
     i5, 10x, 11h,    mit = ,i5, 11h, stopss = , g15.8, &
     10h, stopp = , g15.8/ &
     13x, 11hscale(1) = , g15.8, 11h,  delta = , g15.8, &
     11h, ivaprx = , i5, 11h,   nprt = , i6)
 1370 format (//24h output  -  ifixed(1) = , i6, 9x, 11h, idrvck = , &
     i5, 10x, 11h,    mit = ,i5, 11h, stopss = , g15.8, &
     10h, stopp = , g15.8/ &
     13x, 11hscale(1) = , g15.8, 11h,  delta = , g15.8, &
     11h, ivaprx = , i5, 11h,   nprt = , i5//)
end
