subroutine xstpld ( ldstak,typ )
use M_starpac, only : iprint,ldscmp,lstvec,mdl4,stpls,stpls1,stpls2,stplsc
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XSTPLD tests the derivative checking routines.
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
!        *
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     real exmpt
!        the proportion of observations for which the computed
!        numerical derivatives wrt a given parameter are exempted
!        from meeting the derivative acceptance criteria.
!     integer i
!        an index variable.
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr == 1, errors were detected.
!     integer ixm
!        the first dimension of the independent variable array xm.
!     integer ldsmin
!        the minimum length of the array dstak allowed.
!     integer ldstak
!        the length of the array dstak.
!     integer m
!        the number of independent variables.
!     external mdl4
!        the name of the user supplied subroutine which computes the
!        predicted values based on the current parameter estimates.
!     integer n
!        the number of observations of data.
!     integer neta
!        the number of reliable digits in the model.
!     integer npar
!        the number of unknown parameters in the model.
!     integer nprt
!        the indicator variable used to specify whether or not
!        printed output is to be provided, where if the value of
!        nprt is zero, no printed output is given.
!     integer ntest
!        the number of the current test.
!     real par(10)
!        the array in which the current estimates of the unknown
!        parameters are stored.
!     real scale(10)
!        a dummy array, indicating use of default values for
!        the typical size of the parameters.
!     real stp(10)
!        the selected step sizes for each parameter.
!     real xm(200,2)
!        the array in which one row of the independent variable array
!        is stored.
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
     delta,exmpt
  integer &
     i,ixm,ldsmin,m,n,neta,npar,nprt,ntest
!
!  local arrays
  real &
     par(10),scale(10),stp(10),xm(200,2)
!
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr
!
!     set parameter values
!
  call stpls1(n, m, ixm, par, npar, neta, exmpt, scale, nprt)
  call stpls2(npar, stp)
  call ldscmp(14, 0, 2*(n+npar), 0, 0, 0, 's', 10*n, ldsmin)

  if ( ldstak < ldsmin ) then
    write ( *, 1020) ldsmin
    return
  end if
!
!  create independent variable
!
  delta = 0.0625e0
  xm(1,1) = 0.0e0
  do i=2,n
     xm(i,1) = xm(i-1,1) + delta
  end do

  ntest = 0
!
!  check results from valid calls
!
!  simple example
!
  call stpls1(n, m, ixm, par, npar, neta, exmpt, scale, nprt)

  ntest = ntest + 1
  write ( *,1090) ntest
  write ( *,1040)
  write ( *,1000)
  call stpls2(npar, stp)
  call stpls(xm, n, m, ixm, mdl4, par, npar, ldsmin, stp)
  write ( *,1050) ierr
  write ( *,1080)
  call lstvec(4, stp)

  ntest = ntest + 1
  write ( *,1090) ntest
  write ( *,1040)
  write ( *,1060) neta, exmpt, scale(1), nprt
  write ( *,1010)
  call stpls2(npar, stp)
  call stplsc(xm, n, m, ixm, mdl4, par, npar, ldsmin, stp, neta, &
     exmpt, scale, nprt)
  write ( *,1100) neta, exmpt, scale(1), nprt
  write ( *,1050) ierr
  write ( *,1080)
  call lstvec(4, stp)

  return

 1000 format ('test of stpls' )
 1010 format ('test of stplsc')
 1020 format ('ldstak must be greater than or equal to ', i6)
 1040 format ('simple example')
 1050 format (/'***** returned results *****', 5x, '(-1 indicates ', &
     'value not changed by called subroutine)'//'ierr is ', i3)
 1060 format (' input   -  neta = ', i5, ', exmpt = ', g15.8, &
     ', scale(1) = ', g15.8, ', nprt = ', i5)
 1080 format (//'returned values of stp')
 1090 format ('derivative step size selection subroutine test number', &
     i5)
 1100 format (//'output  -  neta = ', i5, ', exmpt = ', g15.8, &
     ', scale(1) = ', g15.8, ', nprt = ', i5//)
end
