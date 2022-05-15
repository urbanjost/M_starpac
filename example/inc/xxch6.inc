subroutine xxch6 ( ldstak ,typ)
use M_starpac, only : aov1,iprint
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XXCH6 tests the oneway analysis of variance family of routines.
!
!  Discussion:
!
!    Data set is from pages 314-316 of Brownlee [1965].
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Janet Donaldson and John Koontz,
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
!     real tag(20)
!        the tag values for each observation
!     real y(20)
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
         tag(20),y(20)
!
!
!  common blocks
      common /cstak/dstak
      common /errchk/ierr
!
!     data initializations.
!
      data n /16/

      data y( 1), y( 2), y( 3), y( 4) &
          / 83.0, 81.0, 76.0, 78.0/
      data y( 5), y( 6), y( 7), y( 8) &
          / 79.0, 72.0, 61.0, 61.0/
      data y( 9), y(10), y(11), y(12) &
          / 67.0, 67.0, 64.0, 78.0/
      data y(13), y(14), y(15), y(16) &
          / 71.0, 75.0, 72.0, 74.0/

      data tag( 1), tag( 2), tag( 3), tag( 4) &
          / 1.0, 1.0, 1.0, 1.0/
      data tag( 5), tag( 6), tag( 7), tag( 8) &
          / 1.0, 1.0, 2.0, 2.0/
      data tag( 9), tag(10), tag(11), tag(12) &
          / 2.0, 2.0, 2.0, 3.0/
      data tag(13), tag(14), tag(15), tag(16) &
          / 3.0, 3.0, 3.0, 3.0/
!
!
!  print heading
!
      write ( *,1000)
!
!     perform simple test of aov1
!
      write ( *,1100)
      call aov1(y, tag, n, ldstak)
      write ( *,2000) ierr
!
      return
!
!     formats
!
 1000 format ('1*ch6')
 1100 format (' simple test of aov1')
 2000 format (/' the value of ierr is ', i4)
end
