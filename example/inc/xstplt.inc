subroutine xstplt ( ldstak ,typ)
use M_starpac, only : r1mach, iprint,ldscmp,lstvec,mdl4,stpls1,stpls2,stplsc
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XSTPLT tests the derivative checking routines.
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
!     real exmtst(5)
!        various test values for exmpt.
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
!     integer nettst(6)
!        various test values for neta.
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
     exmtst(5),par(10),scale(10),stp(10),xm(200,2)
  integer &
     nettst(6)
!
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr

!
!  set parameter values
!
  call stpls1(n, m, ixm, par, npar, neta, exmpt, scale, nprt)
  call stpls2(npar, stp)
  call ldscmp(14, 0, 2*(n+npar), 0, 0, 0, 's', 10*n, ldsmin)

  if ( ldstak < ldsmin ) then
    write ( *, 1000) ldsmin
    return
  end if
!
!     create independent variable
!
  delta = 0.0625e0
  xm(1,1) = 0.0e0
  do i=2,n
     xm(i,1) = xm(i-1,1) + delta
  end do

  ntest = 0
!
!  test various values of exmpt
!
  call stpls1(n, m, ixm, par, npar, neta, exmpt, scale, nprt)
  exmtst(1) = -1.0e0
  exmtst(2) = 0.0001e0
  exmtst(3) = 0.5e0
  exmtst(4) = 1.0e0
  exmtst(5) = 1.1e0

  do i=1,5

     ntest = ntest + 1
     write ( *,1090) ntest
     write ( *,1040)
     write ( *,1060) neta, exmtst(i), scale(1), nprt
     write ( *,1010)
     call stpls2(npar, stp)
     call stplsc(xm, n, m, ixm, mdl4, par, npar, ldstak, stp, neta, &
        exmtst(i), scale, nprt)
     write ( *,1100) neta, exmtst(i), scale(1), nprt
     write ( *,1050) ierr
     write ( *,1080)
     call lstvec(4, stp)

  end do
!
!  test various values of neta
!
  nettst(1) = -1
  nettst(2) = 0
  nettst(3) = 1
  nettst(4) = 2

  nettst(5) = -log10(r1mach(4))
  nettst(6) = nettst(5) + 1

  scale(1) = 0.0e0

  do i=1,6

     ntest = ntest + 1
     write ( *,1090) ntest
     write ( *,1040)
     write ( *,1060) nettst(i), exmpt, scale(1), nprt
     write ( *,1010)
     call stpls2(npar, stp)
     call stplsc(xm, n, m, ixm, mdl4, par, npar, ldstak, stp, &
        nettst(i), exmpt, scale, nprt)
     write ( *,1100) nettst(i), exmpt, scale(1), nprt
     write ( *,1050) ierr
     write ( *,1080)
     call lstvec(4, stp)

   end do
!
!  suppress output
!
  call stpls1(n, m, ixm, par, npar, neta, exmpt, scale, nprt)
  nprt = 0

  ntest = ntest + 1
  write ( *,1090) ntest
  write ( *,1040)
  write ( *,1060) neta, exmpt, scale(1), nprt
  write ( *,1010)
  call stpls2(npar, stp)
  call stplsc(xm, n, m, ixm, mdl4, par, npar, ldstak, stp, neta, &
     exmpt, scale, nprt)
  write ( *,1100) neta, exmpt, scale(1), nprt
  write ( *,1050) ierr
  write ( *,1080)
  call lstvec(4, stp)
!
!  large calculation error problem
!
  call stpls1(n, m, ixm, par, npar, neta, exmpt, scale, nprt)
  par(3) = 10.0e0**((nettst(5)-1)/2)
  scale(1) = -1.0e0
!
  ntest = ntest + 1
  write ( *,1090) ntest
  write ( *,1070)
  write ( *,1060) neta, exmpt, scale(1), nprt
  write ( *,1010)
  call stpls2(npar, stp)
  call stplsc(xm, n, m, ixm, mdl4, par, npar, ldstak, stp, neta, &
     exmpt, scale, nprt)
  write ( *,1100) neta, exmpt, scale(1), nprt
  write ( *,1050) ierr
  write ( *,1080)
  call lstvec(4, stp)
!
  exmpt = 0.11e0
  nprt = 0
!
  ntest = ntest + 1
  write ( *,1090) ntest
  write ( *,1070)
  write ( *,1060) neta, exmpt, scale(1), nprt
  write ( *,1010)
  call stpls2(npar, stp)
  call stplsc(xm, n, m, ixm, mdl4, par, npar, ldstak, stp, neta, &
     exmpt, scale, nprt)
  write ( *,1100) neta, exmpt, scale(1), nprt
  write ( *,1050) ierr
  write ( *,1080)
  call lstvec(4, stp)
!
  return

 1000 format ('ldstak must be greater than or equal to ', i6)
 1010 format ('test of stplsc')
 1040 format ('simple example')
 1050 format (/'returned results (-1 indicates ', &
     'value not changed by called subroutine)'//' ierr is ', i3)
 1060 format (19h input   -  neta = , i5, 10h, exmpt = , g15.8, &
     13h, scale(1) = , g15.8, 9h, nprt = , i5)
 1070 format (32h large calculation error problem)
 1080 format (//23h returned values of stp)
 1090 format (54h1derivative step size selection subroutine test number, &
     i5)
 1100 format (//19h output  -  neta = , i5, 10h, exmpt = , g15.8, &
     13h, scale(1) = , g15.8, 9h, nprt = , i5//)
end
