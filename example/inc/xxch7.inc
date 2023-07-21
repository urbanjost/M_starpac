subroutine xxch7 ( ldstak ,typ)
use M_starpac, only : corr,iprint
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XXCH7 tests the correlation analysis family of routines.
!
!  Discussion:
!
!    Data is from draper and smith [1968], page 216.
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
!    Draper, Smith,
!    Applied Regression Analysis
!
!  Parameters:
!
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer ierr
!        the integer value designating whether any errors were
!        detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr == 1, errors were detected.
!     integer iym
!        the first dimension of the array ym.
!     integer ldstak
!        the length of dstak in common /cstak/.
!     integer m
!        the number of variables measured for each observation.
!     integer n
!        the number of observations.
!     real ym(10,5)
!        the observed multivariate data.
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
      integer iym,m,n
!
!  local arrays
      real &
         ym(10,5)
!
!
!  common blocks
      common /cstak/dstak
      common /errchk/ierr

      data     ym(1,1),   ym(1,2),   ym(1,3),   ym(1,4) &
          /      42.2e0,  11.2e0,  31.9e0, 167.1e0/
      data     ym(2,1),   ym(2,2),   ym(2,3),   ym(2,4) &
          /      48.6e0,  10.6e0,  13.2e0, 174.4e0/
      data     ym(3,1),   ym(3,2),   ym(3,3),   ym(3,4) &
          /      42.6e0,  10.6e0,  28.7e0, 160.8e0/
      data     ym(4,1),   ym(4,2),   ym(4,3),   ym(4,4) &
          /      39.0e0,  10.4e0,  26.1e0, 162.0e0/
      data     ym(5,1),   ym(5,2),   ym(5,3),   ym(5,4) &
          /      34.7e0,   9.3e0,  30.1e0, 140.8e0/
      data     ym(6,1),   ym(6,2),   ym(6,3),   ym(6,4) &
          /      44.5e0,  10.8e0,   8.5e0, 174.6e0/
      data     ym(7,1),   ym(7,2),   ym(7,3),   ym(7,4) &
          /      39.1e0,  10.7e0,  24.3e0, 163.7e0/
      data     ym(8,1),   ym(8,2),   ym(8,3),   ym(8,4) &
          /      40.1e0,  10.0e0,  18.6e0, 174.5e0/
      data     ym(9,1),   ym(9,2),   ym(9,3),   ym(9,4) &
          /      45.9e0,  12.0e0,  20.4e0, 185.7e0/
!
!
!     set parameters necessary for the computations
!
      iym = 10
      n = 9
      m = 4
!
!     print header
!
      write ( *,1000)
!
!     run simple example of corr
!
      write ( *,1100)
      call corr(ym, n, m, iym, ldstak)
      write ( *,2000) ierr

      return
!
!     format statements
!
 1000 format ('1*ch7')
 1100 format (' simple test of corr')
 2000 format (/' the value of ierr is ', i4)

end
