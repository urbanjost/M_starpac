subroutine xxch5 ( ldstak ,typ)
use M_starpac, only : iprint,stat
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XXCH5 tests the statistical analysis family of routines.
!
!  Discussion:
!
!    Data set is from page 39 of mandel [1964]
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson, John Koontz,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Parameters:
!
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer ierr
!        error flag
!     integer ldstak
!        the length of dstak in common /cstak/.
!     integer n
!        the length of the vector y.
!     real y(40)
!        data vector for tests.
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
      integer n
!
!  local arrays
      real &
         y(40)
!
!
!  common blocks
      common /cstak/dstak
      common /errchk/ierr
!
!     data initializations.
!
      data n /39/
!
      data y( 1), y( 2), y( 3), y( 4) &
          / 0.4, 0.6, 1.0, 1.0/
      data y( 5), y( 6), y( 7), y( 8) &
          / 1.0, 0.5, 0.6, 0.7/
      data y( 9), y(10), y(11), y(12) &
          / 1.0, 0.6, 0.2, 1.9/
      data y(13), y(14), y(15), y(16) &
          / 0.2, 0.4, 0.0, -0.4/
      data y(17), y(18), y(19), y(20) &
          / -0.3, 0.0, -0.4, -0.3/
      data y(21), y(22), y(23), y(24) &
          / 0.1, -0.1, 0.2, -0.5/
      data y(25), y(26), y(27), y(28) &
          / 0.3, -0.1, 0.2, -0.2/
      data y(29), y(30), y(31), y(32) &
          / 0.8, 0.5, 0.6, 0.8/
      data y(33), y(34), y(35), y(36) &
          / 0.7, 0.7, 0.2, 0.5/
      data y(37), y(38), y(39) &
          / 0.7, 0.8, 1.1/
!
!     print heading
!
      write ( *,1000)
!
!     perform simple test of stat
!
      write ( *,1100)
      call stat(y, n, ldstak)
      write ( *,2000) ierr
!
      return
!
!     formats
!
 1000 format ('1*ch5')
 1100 format (' simple test of stat')
 2000 format (/' the value of ierr is ', i4)
end
