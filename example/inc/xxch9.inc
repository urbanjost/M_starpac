subroutine xxch9 ( ldstak ,typ)
use M_starpac, only : dckls,drv1a,drv1b,iprint,mdl1,nls,nlsd,stpls
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XXCH9 tests the nonlinear least squares family of routines.
!
!  Discussion:
!
!    Data is from daniel and wood [1980], pages 428-441.
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
!     external drv1a, drv1b
!        the name of the ''user supplied'' derivative routines.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer ierr
!        the integer value designating whether any errors were
!        detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr == 1, errors were detected.
!     integer ixm
!        the first dimension of the matrix x.
!     integer ldstak
!        the length of dstak in common /cstak/.
!     integer m
!        the number of independent variables.
!     external mdl1
!        the name of the ''user supplied'' model routines.
!     integer n
!        the number of observations in each problem.
!     integer npar
!        the number of parameters to be estimated.
!     real par(5)
!        the array in which the current estimates of the unknown
!        parameters are stored.
!     real res(10)
!        the residuals.
!     real stp(5)
!        the step sizes selected for generating finite difference
!        derivatives.
!     real xm(10,2)
!        the independent variable.
!     real y(10)
!        the dependent variable.
!
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
      integer ixm,m,n,npar
!
!  local arrays
      real &
         par(5),res(10),stp(5),xm(10,2),y(10)
!
!
!  common blocks
      common /cstak/dstak
      common /errchk/ierr

      data y(1), y(2), y(3), y(4), y(5), y(6) &
         /2.138e0, 3.421e0, 3.597e0, 4.340e0, 4.882e0, 5.660e0/

      data xm(1,1), xm(2,1), xm(3,1), xm(4,1), xm(5,1), xm(6,1) &
         /1.309e0, 1.471e0, 1.490e0, 1.565e0, 1.611e0, 1.680e0/
!
!  set parameters necessary for the computations
!
      ixm = 10
      n = 6
      m = 1
      npar = 2
!
!  print header
!
      write ( *,1000)
!
!  run simple example of nls
!
      write ( *,1100)
      par(1) = 0.725
      par(2) = 4.000
      call nls(y, xm, n, m, ixm, mdl1, par, npar, res, ldstak)
      write ( *,2000) ierr
!
!  run simple example of nlsd
!
      write ( *,1200)
      par(1) = 0.725
      par(2) = 4.000
      call nlsd(y, xm, n, m, ixm, mdl1, drv1a, par, npar, res, ldstak)
      write ( *,2000) ierr
!
!  run simple example of stpls
!
      write ( *,1300)
      par(1) = 0.725
      par(2) = 4.000
      call stpls(xm, n, m, ixm, mdl1, par, npar, ldstak, stp)
      write ( *,2000) ierr
!
!  run simple example of dckls
!
      write ( *,1400)
      par(1) = 0.000
      par(2) = 4.000
      call dckls(xm, n, m, ixm, mdl1, drv1b, par, npar, ldstak)
      write ( *,2000) ierr

      return

 1000 format ('1*ch9')
 1100 format (' simple test of nls')
 1200 format ('1simple test of nlsd')
 1300 format ('1simple test of stpls')
 1400 format ('1simple test of dckls')
 2000 format (/' the value of ierr is ', i4)
end
