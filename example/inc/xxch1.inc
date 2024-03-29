subroutine xxch1 ( ldstak ,typ)
use M_starpac, only : iprint,pp,stat
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XXCH1 tests the page plot and statistical analysis families of routines.
!
!  Discussion:
!
!    The data set is 84 relative humidity measurememts from pikes peak.
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
!     integer i
!        an index variable.
!     integer ierr
!        error flag
!     integer ldstak
!        the length of dstak in common /cstak/.
!     integer n
!        the length of the vector y.
!     real x(100)
!        the order indices of the data.
!     real y(100)
!        data vector for tests.
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
      integer &
         i,n
!
!  local arrays
      real &
         x(100),y(100)
!
!
!  common blocks
      common /cstak/dstak
      common /errchk/ierr

      data n /84/
!
      data y( 1), y( 2), y( 3), y( 4) &
          / 0.6067e0, 0.6087e0, 0.6086e0, 0.6134e0/
      data y( 5), y( 6), y( 7) &
          / 0.6108e0, 0.6138e0, 0.6125e0/
      data y( 8), y( 9), y(10), y(11) &
          / 0.6122e0, 0.6110e0, 0.6104e0, 0.7213e0/
      data y(12), y(13), y(14) &
          / 0.7078e0, 0.7021e0, 0.7004e0/
      data y(15), y(16), y(17), y(18) &
          / 0.6981e0, 0.7242e0, 0.7268e0, 0.7418e0/
      data y(19), y(20), y(21) &
          / 0.7407e0, 0.7199e0, 0.6225e0/
      data y(22), y(23), y(24), y(25) &
          / 0.6254e0, 0.6252e0, 0.6267e0, 0.6218e0/
      data y(26), y(27), y(28) &
          / 0.6178e0, 0.6216e0, 0.6192e0/
      data y(29), y(30), y(31), y(32) &
          / 0.6191e0, 0.6250e0, 0.6188e0, 0.6233e0/
      data y(33), y(34), y(35) &
          / 0.6225e0, 0.6204e0, 0.6207e0/
      data y(36), y(37), y(38), y(39) &
          / 0.6168e0, 0.6141e0, 0.6291e0, 0.6231e0/
      data y(40), y(41), y(42) &
          / 0.6222e0, 0.6252e0, 0.6308e0/
      data y(43), y(44), y(45), y(46) &
          / 0.6376e0, 0.6330e0, 0.6303e0, 0.6301e0/
      data y(47), y(48), y(49) &
          / 0.6390e0, 0.6423e0, 0.6300e0/
      data y(50), y(51), y(52), y(53) &
          / 0.6260e0, 0.6292e0, 0.6298e0, 0.6290e0/
      data y(54), y(55), y(56) &
          / 0.6262e0, 0.5952e0, 0.5951e0/
      data y(57), y(58), y(59), y(60) &
          / 0.6314e0, 0.6440e0, 0.6439e0, 0.6326e0/
      data y(61), y(62), y(63) &
          / 0.6392e0, 0.6417e0, 0.6412e0/
      data y(64), y(65), y(66), y(67) &
          / 0.6530e0, 0.6411e0, 0.6355e0, 0.6344e0/
      data y(68), y(69), y(70) &
          / 0.6623e0, 0.6276e0, 0.6307e0/
      data y(71), y(72), y(73), y(74) &
          / 0.6354e0, 0.6197e0, 0.6153e0, 0.6340e0/
      data y(75), y(76), y(77) &
          / 0.6338e0, 0.6284e0, 0.6162e0/
      data y(78), y(79), y(80), y(81) &
          / 0.6252e0, 0.6349e0, 0.6344e0, 0.6361e0/
      data y(82), y(83), y(84) &
          / 0.6373e0, 0.6337e0, 0.6383e0/

      do i=1,n
         x(i) = i
      end do
!
!     print heading
!
      write ( *,1000)
!
!     perform simple test of pp
!
      write ( *,1100)
      call pp(y, x, n)
      write ( *,2000) ierr
!
!     perform simple test of stat
!
      write ( *,1200)
      call stat(y, n, ldstak)
      write ( *,2000) ierr

      return
!
!     formats
!
 1000 format ('1*ch1')
 1100 format (' simple test of pp')
 1200 format ('1simple test of stat')
 2000 format (/' the value of ierr is ', i4)
end
