subroutine xhist ( ldstak ,typ)
use M_starpac, only : hist,histc,iprint,ldscmp,nrand
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XHIST tests the HIST family.
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
!        a loop index.
!     integer ierr
!        flag to indicate presence of error detected by preceding
!        starpac call.  (0 is ok, 1 is error)
!     integer ldsmin
!        the minimum amount of work area needed for a given problem.
!     integer ldstak
!        amount of work area.  size of dstak.
!     integer n
!        the length of the vector y.
!     integer ncell
!        the user supplied value for the number of cells in the
!        histogram.  if ncell is less than or equal to zero, the
!        number of cells to be used (ncells) will be calculated from
!        the number of observations.
!     integer nconst
!        length of the vector yconst.
!     integer nprtof
!        flag for no output (except error messages).
!     integer nprton
!        flag for full printout.
!     real y(84)
!        data vector for tests.
!     real yconst(10)
!        vector of constant data.
!     real ylb
!        the lower bound for selecting data from y for the histogram.
!     real ylong(200)
!        long vector of data
!     real ypath(10)
!        a vector of y values designed to force different paths
!        through the summation routines.
!     real yub
!        the upper bound for selecting data from y for the histogram.
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
     ylb,yub
  integer ldsmin,n,ncell,nconst
!
!  local arrays
  real &
     y(84),yconst(10),ylong(200),ypath(10)
!
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr

!
!     data initializations.
!
  data n /84/
  data nconst /10/
  data ncell/10/
  data ylb/0.60e0/, yub/0.63e0/
!
!     davis-harrison r.h. data, pikes peak.
!
!     this is an arbitrarily chosen data set.
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
!
!  check for sufficient work area length.
!
  if (ldstak.lt.300) then
    write ( *, 1000)
     return
  end if

  yconst(1:nconst) = 1.0e0
!
!  heading.
!
  write ( *,1150)
!
!  test 1.  check all error messages.
!
  write ( *,1160)
!
!  error 1, zero or fewer elements.
!
  write ( *,1180)
  call hist(y, 0, ldstak)
  write ( *, 1350)
  write ( *, 1360) y(1:n)
  write ( *,1170) ierr
  call histc(y, 0, ncell, ylb, yub, ldstak)
  write ( *, 1350)
  write ( *, 1360) y(1:n)
  write ( *,1170) ierr
!
!  error 2, not enough space in cstak.
!
  write ( *,1190)
  call ldscmp(2, 0, n, 0, 0, 0, 's', &
              min(nint(5.5+1.5*anint(log10(real(n)))),25),ldsmin)
  call hist(y, n, ldsmin-1)
  write ( *,1170) ierr
  write ( *,1195)
  call hist(y, n, ldsmin)
  write ( *,1170) ierr
  write ( *,1190)
  call ldscmp(2, 0, n, 0, 0, 0, 's', ncell, ldsmin)
  call histc(y, n, ncell, ylb, yub, ldsmin-1)
  write ( *,1170) ierr
  write ( *,1195)
  call histc(y, n, ncell, ylb, yub, ldsmin)
  write ( *,1170) ierr
!
!  constant y. (not an error)
!
  write ( *,1200)
  call hist(yconst, nconst, ldstak)
  write ( *,1170) ierr
  write ( *,1200)
  call histc(yconst, nconst, ncell, ylb, yub, ldstak)
  write ( *,1170) ierr
!
!  error 4, no data within user supplied limits
!
  write ( *, 1110)
  call histc(y, n, 0, 4.0e0, 10.0e0, ldstak)
  write ( *, 1170) ierr
!
!  test 2.  make a working run of each routine to check the output.
!
  write ( *,1300)
  write ( *,1310)
  call hist(y, n, ldstak)
  write ( *, 1350)
  write ( *, 1360) y(1:n)
  write ( *,1170) ierr

  write ( *,1340)
  call histc(y, n, ncell, ylb, yub, ldstak)
  write ( *, 1350)
  write ( *, 1360) y(1:n)
  write ( *,1170) ierr
!
!  run data set 6.7.
!
  ypath(1:10) = 0.0e0
  ypath(1) = -1.0e0
  ypath(10) = 1.0e0
  write ( *,1130)
  call hist(ypath, nconst, ldstak)
  write ( *,1170) ierr
  write ( *, 1130)
  call histc(ypath, nconst, 0, 0.0e0, 0.0e0, ldstak)
  write ( *, 1130)
  call histc(ypath, nconst, 1, 0.0e0, 0.0e0, ldstak)
  write ( *, 1130)
  call histc(ypath, nconst, 0, -0.5e0, 0.5e0, ldstak)
  write ( *, 1130)
  call histc(ypath, nconst, 0, 1.0e0, 4.0e0, ldstak)
!
!  run data set 6.8
!
  write ( *, 1120)
  call nrand (ylong, 200, 3254767)
  call hist (ylong, 200, ldstak)
  return

 1000 format ('1the dimension of dstak and the value of ldstak needed'/ &
    ' for histx must equal or exceed 300.  change driver'/ &
    ' and recall histx.')
 1110 format ('try no data within user supplied limits.')
 1120 format ('run hist on 200 pseudo-randon numbers')
 1130 format('run hist on -1, 8*0, 1.')
 1150 format ('test runs for the histogram family of routines.')
 1160 format('test 1.  generate one of each of the possible ', &
     'error messages.')
 1170 format('the value of ierr is ', i4)
 1180 format('try zero or fewer elements.')
 1190 format('1test with insufficient work area')
 1195 format(' test with exactly the right amount of work area.')
 1200 format('1try constant y. (not an error)')
 1300 format('test 4.  make working runs of all routines to check', &
     ' the output.')
 1310 format('run hist on the davis-harrison pikes peak data.')
 1340 format('run histc on the davis-harrison pikes peak data.')
 1350 format(/' print the data to insure the original order has', &
     ' been restored.')
 1360 format (7f10.5)
end
