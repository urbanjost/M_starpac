subroutine xlls ( lds ,typ)
use M_starpac, only : fitsxp,genr,iprint,ldscmp,lls,llsp,llsps,llspw,llspws, llss,llsw,llsws,nrand,setrv
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XLLS tests the linear least squares subroutines.
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
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer i
!        an index.
!     integer ierr
!        the integer value designating whether any errors were
!        detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr == 1, errors were detected.
!     integer ivcv
!        the first dimension of the matrix vcv.
!     integer ixm
!        the first dimension of the matrix x.
!     integer j
!        an index.
!     integer lds
!       ..
!     integer ldsmin
!        the minimum length allowed for the array dstak.
!     integer ldstak
!        the length of the vector dstak in common cstak.
!     integer lpar
!        the actual length of the parameter array.
!     integer n
!        the number of observations.
!     integer ndeg
!        the degree of the polynomial model to be fit.
!     integer npar
!        the number of parameters.
!     integer nprt
!        the indicator variable used to designate the amount of
!        printed output.
!     real par(10)
!        the parameters  to be estimated.
!     real pv(50)
!        the predicted values.
!     real rand(1)
!        *
!     real res(50)
!        the residuals.
!     real rsd
!        the residual standard deviation.
!     real sdpv(50)
!        the standard deviations of the predicted values.
!     real sdres(50)
!        the standardized residuals.
!     real vcv(10,10)
!        the variance covariance matrix.
!     real wt(50)
!        the weights (a dummy vector in the unweighted case).
!     real x(50,9)
!        the independent variable.
!     real xm(50,10)
!        the independent variable.
!     real xm1(50,10)
!        the independent variable.
!     real y(50)
!        the dependent variable.
!     real y1(50)
!        the dependent variable.
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
     rsd,sum,term
  integer &
     i,ivcv,ixm,j,ldsmin,ldstak,lpar,n,ndeg,npar,nprt
!
!  local arrays
  real &
     par(10),pv(50),rand(1),res(50),sdpv(50),sdres(50),vcv(10,10), &
     wt(50),x(50,9),xm(50,10),xm1(50,10),y(50),y1(50)
!
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr
!
!  equivalences
  equivalence (xm(1,2),x(1,1))
!
  data      xm(1,1),  xm(1,2),  xm(1,3),  xm(1,4) &
      /      1.0e0, 42.2e0, 11.2e0, 31.9e0/
  data      xm(2,1),  xm(2,2),  xm(2,3),  xm(2,4) &
      /      1.0e0, 48.6e0, 10.6e0, 13.2e0/
  data      xm(3,1),  xm(3,2),  xm(3,3),  xm(3,4) &
      /      1.0e0, 42.6e0, 10.6e0, 28.7e0/
  data      xm(4,1),  xm(4,2),  xm(4,3),  xm(4,4) &
      /      1.0e0, 39.0e0, 10.4e0, 26.1e0/
  data      xm(5,1),  xm(5,2),  xm(5,3),  xm(5,4) &
      /      1.0e0, 34.7e0,  9.3e0, 30.1e0/
  data      xm(6,1),  xm(6,2),  xm(6,3),  xm(6,4) &
      /      1.0e0, 44.5e0, 10.8e0,  8.5e0/
  data      xm(7,1),  xm(7,2),  xm(7,3),  xm(7,4) &
      /      1.0e0, 39.1e0, 10.7e0, 24.3e0/
  data      xm(8,1),  xm(8,2),  xm(8,3),  xm(8,4) &
      /      1.0e0, 40.1e0, 10.0e0, 18.6e0/
  data      xm(9,1),  xm(9,2),  xm(9,3),  xm(9,4) &
      /      1.0e0, 45.9e0, 12.0e0, 20.4e0/
  data         y(1),     y(2),     y(3) &
      /    167.1e0,174.4e0,160.8e0/
  data         y(4),     y(5),     y(6) &
      /    162.0e0,140.8e0,174.6e0/
  data         y(7),     y(8),     y(9) &
      /    163.7e0,174.5e0,185.7e0/
!
!  set parameters necessary for the computations
!
  n = 9
  npar = 4
  ndeg = 3
  nprt = 2
  lpar = 10
  ivcv = 10
  ixm = 50
  ldstak = lds

  call setrv(wt, n, 1.0e0)
!
!  check error handling
!
!  error 1  -  non positive number of observations and parameter
!  number of parameters greater than n
!  ixm less than number of observations
!  ivcv less than number of parameters
!  lpar too small
!
  n = -5
  npar = 0
  ndeg = -1
  ixm = -10
  lpar = -1
  ivcv = -10
  nprt = -1
  write ( *,1200)
  write ( *,1000)
  call lls(y, xm, n, ixm, npar, res, ldstak)
  write ( *,1500) ierr
  write ( *,1010)
  call llss(y, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  write ( *,1020)
  call llsw(y, wt, xm, n, ixm, npar, res, ldstak)
  write ( *,1500) ierr
  write ( *,1030)
  call llsws(y, wt, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv,sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  write ( *,1040)
  call llsp(y, x, n, ndeg, res, ldstak)
  write ( *,1500) ierr
  write ( *,1050)
  call llsps(y, x, n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv,sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  write ( *,1060)
  call llspw(y, wt, x, n, ndeg, res, ldstak)
  write ( *,1500) ierr
  write ( *,1070)
  call llspws(y, wt, x, n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  n = 9
  npar = 4
  ndeg = 3
  ixm = 50
  lpar = -10
  ivcv = 10
!
!  error 2  -  lds too small, lpar too small
!
  ldstak = 0
  write ( *,1220)
  write ( *,1000)
  call lls(y, xm, n, ixm, npar, res, ldstak)
  write ( *,1500) ierr
  write ( *,1010)
  call llss(y, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  write ( *,1020)
  call llsw(y, wt, xm, n, ixm, npar, res, ldstak)
  write ( *,1500) ierr
  write ( *,1030)
  call llsws(y, wt, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv,sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  write ( *,1040)
  call llsp(y, x, n, ndeg, res, ldstak)
  write ( *,1500) ierr
  write ( *,1050)
  call llsps(y, x, n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv,sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  write ( *,1060)
  call llspw(y, wt, x, n, ndeg, res, ldstak)
  write ( *,1500) ierr
  write ( *,1070)
  call llspws(y, wt, x, n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  ldstak = lds
  nprt = 2
  lpar = 10
!
!  error 3  -  negative weights
!
  wt(1) = -1.0e0
  write ( *,1240)
  write ( *,1020)
  call llsw(y, wt, xm, n, ixm, npar, res, ldstak)
  write ( *,1500) ierr
  write ( *,1030)
  call llsws(y, wt, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv,sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  write ( *,1060)
  call llspw(y, wt, x, n, ndeg, res, ldstak)
  write ( *,1500) ierr
  write ( *,1070)
  call llspws(y, wt, x, n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  wt(1) = 1.0e0
!
!  error 4  -  too few positive weights
!
  call setrv(wt(2), n-1, 0.0e0)
  write ( *,1250)
  write ( *,1020)
  call llsw(y, wt, xm, n, ixm, npar, res, ldstak)
  write ( *,1500) ierr
  write ( *,1030)
  call llsws(y, wt, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv,sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  write ( *,1060)
  call llspw(y, wt, x, n, ndeg, res, ldstak)
  write ( *,1500) ierr
  write ( *,1070)
  call llspws(y, wt, x, n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call setrv(wt(2), n-1, 1.0e0)
!
!  check results from valid call
!
  write ( *,1260)
  write ( *,1000)
  call lls(y, xm, n, ixm, npar, res, ldstak)
  write ( *,1500) ierr
  write ( *,1430) res(1:n)
  write ( *,1260)
  write ( *,1010)
  call llss(y, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
  write ( *,1260)
  write ( *,1020)
  call llsw(y, wt, xm, n, ixm, npar, res, ldstak)
  write ( *,1500) ierr
  write ( *,1430) res(1:n)
  write ( *,1260)
  write ( *,1030)
  call llsws(y, wt, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv,sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
  write ( *,1260)
  write ( *,1040)
  call llsp(y, x, n, ndeg, res, ldstak)
  write ( *,1500) ierr
  write ( *,1430) res(1:n)
  write ( *,1260)
  write ( *,1050)
  call llsps(y, x, n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv,sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
  write ( *,1260)
  write ( *,1060)
  call llspw(y, wt, x, n, ndeg, res, ldstak)
  write ( *,1500) ierr
  write ( *,1430) res(1:n)
  write ( *,1260)
  write ( *,1070)
  call llspws(y, wt, x, n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
!
!  check results from exact fit
!
  n = npar
  ndeg = npar-1

  write ( *,1270)
  write ( *,1000)
  call lls(y, xm, n, ixm, npar, res, ldstak)
  write ( *,1500) ierr
  write ( *,1430) res(1:n)
  write ( *,1270)
  write ( *,1010)
  call llss(y, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
  write ( *,1270)
  write ( *,1040)
  call llsp(y, x, n, ndeg, res, ldstak)
  write ( *,1500) ierr
  write ( *,1430) res(1:n)
  write ( *,1270)
  write ( *,1050)
  call llsps(y, x, n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv,sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
!
  n = 9
!
  call setrv(wt(npar+1), n-npar, 0.0e0)
!
  write ( *,1270)
  write ( *,1020)
  call llsw(y, wt, xm, n, ixm, npar, res, ldstak)
  write ( *,1500) ierr
  write ( *,1430) res(1:n)
  write ( *,1270)
  write ( *,1030)
  call llsws(y, wt, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv,sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
  write ( *,1270)
  write ( *,1060)
  call llspw(y, wt, x, n, ndeg, res, ldstak)
  write ( *,1500) ierr
  write ( *,1430) res(1:n)
  write ( *,1270)
  write ( *,1070)
  call llspws(y, wt, x, n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
!
  call setrv(wt(npar+1), n-npar, 1.0e0)
!
!     check results from rank deficient fit
!
  xm(1:n,5) = xm(1:n,4)
  write ( *,1280)
  write ( *,1000)
  call lls(y, xm, n, ixm, npar+1, res, ldstak)
  write ( *,1500) ierr
!
!     check results from a poorly scaled problem.
!
  do i = 1, n
     y1(i) = y(i) * 1.0e-8
     xm1(i,1:4) = xm(i,1:4)
     xm1(i,3) = xm1(i,3) * 1.0e+8
  end do

  write ( *,1290)
  write ( *,1000)
  call lls(y1, xm, n, ixm, npar, res, ldstak)
  write ( *,1500) ierr
  write ( *,1290)
  write ( *,1000)
  call lls(y, xm1, n, ixm, npar, res, ldstak)
  write ( *,1500) ierr
  write ( *,1290)
  write ( *,1000)
  call lls(y1, xm1, n, ixm, npar, res, ldstak)
  write ( *,1500) ierr
!
!  minimum amount of work area.
!
  call ldscmp(15, 0, 0, 0, 0, 0, 's', &
              6*n + npar*(n+2*npar+5) + 1, ldsmin)

  write ( *,1300)
  write ( *,1000)
  call lls(y, xm, n, ixm, npar, res, ldsmin)
  write ( *,1500) ierr
  write ( *,1430) res(1:n)
!
!  check results for weighted analysis
!
  nprt = 1111
  call setrv(wt, n, 100.0e0)
  write ( *,1310)
  write ( *,1030)
  call llsws(y, wt, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv,sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)

  wt(1) = 0.0e0
  wt(5) = 0.0e0
  wt(9) = 0.0e0

  write ( *,1310)
  write ( *,1030)
  call llsws(y, wt, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv,sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
  call setrv(wt, n, 100.0e0)

  call genr(wt, n, 1.0e0, 1.0e0)
  write ( *,1310)
  write ( *,1030)
  call llsws(y, wt, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv,sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
  call setrv(wt, n, 100.0e0)
!
!  check print control
!
  nprt = 1000
  write ( *,1320) nprt
  write ( *,1010)
  call llss(y, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)

  nprt = 2000
  write ( *,1320) nprt
  write ( *,1010)
  call llss(y, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)

  nprt = 200
  write ( *,1320) nprt
  write ( *,1010)
  call llss(y, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)

  nprt = 20
  write ( *,1320) nprt
  write ( *,1010)
  call llss(y, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)

  nprt = 2
  write ( *,1320) nprt
  write ( *,1010)
  call llss(y, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)

  nprt = 0
  write ( *,1320) nprt
  write ( *,1010)
  call llss(y, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
!
!  check results for n = 2, npar = id+1 = 1
!
  nprt = 2222
  n = 2
  npar = 1
  ndeg = 0
  write ( *,1330)
  write ( *,1010)
  call llss(y, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
  write ( *,1330)
  write ( *,1070)
  call llspws(y, wt, x, n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
!
!  check results for n = 1, npar = id+1 = 1
!
  nprt = 2222
  n = 1
  npar = 1
  ndeg = 0
  write ( *,1330)
  write ( *,1010)
  call llss(y, xm, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
  write ( *,1330)
  write ( *,1070)
  call llspws(y, wt, x, n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
  n = 9
  npar = 4
  ndeg = 3
!
!  ill-conditioned
!
  do i = 1, 50
     term = 1.0e0
     sum = 0.0e0
     do j = 1, 6
        xm1(i,j) = term
        sum = sum + term
        term = (i-1)*term
     end do
     y1(i) = sum
  end do

  n = 21
  npar = 6
  ndeg = 5
  write ( *,1340)
  write ( *,1010)
  call llss(y1, xm1, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
  write ( *,1340)
  write ( *,1050)
  call llsps(y1, xm1(1,2), n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)

  n = 50
  npar = 6
  ndeg = 5
  call nrand(rand, 1, 223)
  do i = 1, n
     call nrand(rand, 1, 0)
     y1(i) = y1(i) + rand(1)
  end do
  write ( *,1340)
  write ( *,1010)
  call llss(y1, xm1, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
  write ( *,1340)
  write ( *,1050)
  call llsps(y1, xm1(1,2), n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
!
  nprt = 1000
  write ( *,1340)
  write ( *,1010)
  call llss(y1, xm1, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)
  write ( *,1340)
  write ( *,1050)
  call llsps(y1, xm1(1,2), n, ndeg, res, ldstak, &
     nprt, lpar, par, npar, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)

  n = 45
  call setrv(wt, n, 1.0e0)
  write ( *,1340)
  write ( *,1010)
  call llsws(y1, wt, xm1, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)

  n = 44
  write ( *,1340)
  write ( *,1010)
  call llsws(y1, wt, xm1, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)

  n = 41
  write ( *,1340)
  write ( *,1010)
  call llsws(y1, wt, xm1, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)

  n = 40
  write ( *,1340)
  write ( *,1010)
  call llsws(y1, wt, xm1, n, ixm, npar, res, ldstak, &
     nprt, par, rsd, pv, sdpv, sdres, vcv, ivcv)
  write ( *,1500) ierr
  call fitsxp(par, pv, sdpv, res, sdres, vcv, n, npar, ivcv, rsd)

  return

 1000 format (' call to lls   ')
 1010 format (' call to llss  ')
 1020 format (' call to llsw  ')
 1030 format (' call to llsws ')
 1040 format (' call to llsp  ')
 1050 format (' call to llsps ')
 1060 format (' call to llspw ')
 1070 format (' call to llspws')
 1200 format ('1miscellaneous errors  -  test 1')
 1220 format ('1miscellaneous errors  -  test 2')
 1240 format ('1negative weights')
 1250 format ('1too few positive weights')
 1260 format ('1valid problem')
 1270 format ('1zero residual problem')
 1280 format ('1rank deficient problem')
 1290 format ('1poorly scaled problem')
 1300 format ('1minimum work area size')
 1310 format ('1weighted analysis')
 1320 format ('1check print control  -  nprt = ', i5)
 1330 format ('1check minimum problem size')
 1340 format ('1ill-conditioned problem')
 1430 format (//4h res/ (1x, e22.14))
 1500 format (/' ierr = ', i5)
end
