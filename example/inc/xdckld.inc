subroutine xdckld ( ldstak ,typ)
use M_starpac, only : dckls,dckls1,dcklsc,drv4a,iprint,ldscmp
use M_starpac, only : mdl4
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XDCKLD tests the derivative checking routines.
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
!     external drv4a
!        the name of the user supplied subroutine which computes the
!        analytic derivatives (jacobian matrix) of the model.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
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
!     integer nrow
!        the number of the row of the independent variable array at
!        which the derivative is to be checked.
!     integer ntau
!        the number of digits of agreement required between the
!        numerical derivatives and the user supplied derivatives.
!     integer ntest
!        the number of the current test.
!     real par(10)
!        the array in which the current estimates of the unknown
!        parameters are stored.
!     real scale(10)
!        a dummy array, indicating use of default values for
!        the typical size of the parameters.
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
     delta
  integer &
     i,ixm,ldsmin,m,n,neta,npar,nprt,nrow,ntau,ntest
!
!  local arrays
  real &
     par(10),scale(10),xm(200,2)
!
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr
!
!  set parameter values
!
  call dckls1(n, m, ixm, par, npar, neta, ntau, nrow, scale, nprt)
  call ldscmp(5, 0, 2*npar+1, 0, 0, 0, 's', &
              n*npar+npar+n, ldsmin)

  if (ldsmin.gt.ldstak) then
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
!  check result for correctly computed derivative
!
  call dckls1(n, m, ixm, par, npar, neta, ntau, nrow, scale, nprt)

  ntest = ntest + 1
  write ( *,1130) ntest
  write ( *,1100)
  write ( *,1040)
  write ( *,1000)
  ierr = -1
  call dckls(xm, n, m, ixm, mdl4, drv4a, par, npar, ldsmin)
  write ( *,1050) ierr

  ntest = ntest + 1
  write ( *,1130) ntest
  write ( *,1100)
  write ( *,1040)
  write ( *,1060) neta, ntau, scale(1), nrow, nprt
  write ( *,1010)
  ierr = -1
  call dcklsc(xm, n, m, ixm, mdl4, drv4a, par, npar, ldsmin, neta, &
     ntau, scale, nrow, nprt)
  write ( *,1050) ierr
  write ( *,1140) neta, ntau, scale(1), nrow, nprt

  return

 1000 format (' test of dckls' )
 1010 format ('test of dcklsc')
 1020 format ('ldstak must be greater than or equal to ', i6)
 1040 format ('simple example')
 1050 format ('***** returned results *****', 5x, ' (-1 indicates ', &
     'value not changed by called subroutine)'//' ierr is ', i3)
 1060 format (' input   -  neta = ', i5, ', ntau = ', i5, &
     ', scale(1) = ', g15.8, ', nrow = ', i5, ', nprt = ', i5)
 1100 format (' correctly coded derivative')
 1130 format ('derivative checking subroutine test number', i5/)
 1140 format ('output  -  neta = ', i5, ', ntau = ', i5, &
     ', scale(1) = ', g15.8, ', nrow = ', i5, ', nprt = ', i5//)
end
