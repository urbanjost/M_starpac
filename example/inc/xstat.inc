subroutine xstat ( ldstak ,typ)
use M_starpac, only : r1mach, iprint,stat,stats,statw,statws
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XSTAT tests the STAT family.
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
!     real fplm
!        the floating point largest magnitude.
!     integer i
!        a loop index.
!     integer ierr
!        flag to indicate presence of error detected by preceding
!        starpac call.  (0 is ok, 1 is error)
!     integer ldstak
!        amount of work area.  size of dstak.
!     integer n
!        the length of the vector y.
!     integer nconst
!        length of the vector yconst.
!     integer nprtof
!        flag for no output (except error messages).
!     integer nprton
!        flag for full printout.
!     real sts(53)
!        vector of statistics.
!     real wt(84)
!        weights vector.
!     real wtall0(10)
!        n vector of 0 weights.
!     real wtall1(84)
!        n vector of 1 weights.
!     real wtemp
!        temporary storage for one of the weights.
!     real y(84)
!        data vector for tests.
!     real yconst(10)
!        vector of constant data.
!     real ypath(10)
!        a vector of y values designed to force different paths
!        through the summation routines.
!     real ytempn, ytemp1
!        temporary storage for the first and last y value.
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
     fplm,wtemp,ytemp1,ytempn
  integer &
     i,n,nconst,nprtof,nprton
!
!  local arrays
  real &
     sts(53),wt(84),wtall0(10),wtall1(84),y(84),yconst(10), &
     ypath(10)
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
  data nprton /1/
  data nprtof /0/
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
  data wt( 1), wt( 2), wt( 3), wt( 4), wt( 5), wt( 6), wt( 7) &
     / 0.5e0, 0.5e0, 0.5e0, 0.5e0, 0.5e0, 0.5e0, 0.5e0/
  data wt( 8), wt( 9), wt(10), wt(11), wt(12), wt(13), wt(14) &
     / 0.5e0, 0.5e0, 0.5e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0/
  data wt(15), wt(16), wt(17), wt(18), wt(19), wt(20), wt(21) &
     / 0.0e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0, 0.0e0, 0.5e0/
  data wt(22), wt(23), wt(24), wt(25), wt(26), wt(27), wt(28) &
     / 0.5e0, 0.5e0, 0.5e0, 0.5e0, 0.5e0, 0.5e0, 0.5e0/
  data wt(29), wt(30), wt(31), wt(32), wt(33), wt(34), wt(35) &
     / 0.5e0, 0.5e0, 0.5e0, 0.5e0, 0.5e0, 0.5e0, 0.5e0/
  data wt(36), wt(37), wt(38), wt(39), wt(40), wt(41), wt(42) &
     / 0.5e0, 0.5e0, 0.5e0, 0.5e0, 0.5e0, 0.5e0, 0.5e0/
  data wt(43), wt(44), wt(45), wt(46), wt(47), wt(48), wt(49) &
     / 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0/
  data wt(50), wt(51), wt(52), wt(53), wt(54), wt(55), wt(56) &
     / 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0, 0.0e0, 0.0e0/
  data wt(57), wt(58), wt(59), wt(60), wt(61), wt(62), wt(63) &
     / 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0/
  data wt(64), wt(65), wt(66), wt(67), wt(68), wt(69), wt(70) &
     / 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0/
  data wt(71), wt(72), wt(73), wt(74), wt(75), wt(76), wt(77) &
     / 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0/
  data wt(78), wt(79), wt(80), wt(81), wt(82), wt(83), wt(84) &
     / 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0, 1.0e0/

  fplm = r1mach(2)
!
!  set up the weights vectors.
!
  wtall1(1:n) = 1.0e0

  yconst(1:nconst) = 1.0e0
  wtall0(1:nconst) = 0.0e0
!
!  heading.
!
  write ( *,1150)
!
!  test 1.  check all error messages.
!
!  error 1, two or fewer elements.
!
  write ( *,1180)
  write ( *,1230)
  write ( *,1240)
  call stat(y, 2, ldstak)
  write ( *,1390) y(1:10)
  write ( *,1170) ierr
  write ( *,1230)
  write ( *,1250)
  call stats(y, 2, ldstak, sts, nprton)
  write ( *,1390) y(1:10)
  write ( *,1170) ierr
  write ( *,1230)
  write ( *,1400)
  call statw(y, wt, 2, ldstak)
  write ( *,1390) y(1:10)
  write ( *,1170) ierr
  write ( *,1230)
  write ( *,1410)
  call statws(y, wt, 2, ldstak, sts, nprton)
  write ( *,1390) y(1:10)
  write ( *,1170) ierr
!
!  error 2, not enough space in cstak.
!
  write ( *,1190)
  write ( *,1230)
  write ( *,1240)
  call stat(y, n, n/4)
  write ( *,1170) ierr
  write ( *,1230)
  write ( *,1250)
  call stats(y, n, n/4, sts, nprton)
  write ( *,1170) ierr
  write ( *,1230)
  write ( *,1400)
  call statw(y, wt, n, n/4)
  write ( *,1170) ierr
  write ( *,1230)
  write ( *,1410)
  call statws(y, wt, n, n/4, sts, nprton)
  write ( *,1170) ierr
!
!  error 4, negative weights.
!
  write ( *,1210)
  wtemp = wt(2)
  wt(2) = -1.0e0
  write ( *,1230)
  write ( *,1400)
  call statw(y, wt, n, ldstak)
  write ( *,1390) y(1:10)
  write ( *,1390) y(1:10)
  write ( *,1170) ierr
  write ( *,1230)
  write ( *,1410)
  call statws(y, wt, n, ldstak, sts, nprton)
  write ( *,1390) y(1:10)
  write ( *,1390) y(1:10)
  write ( *,1170) ierr
  wt(2) = wtemp
!
!     error 5, all weights zero (plus constant y).
!
  write ( *,1220)
  write ( *,1230)
  write ( *,1400)
  call statw(yconst, wtall0, nconst, ldstak)
  write ( *,1170) ierr
  write ( *,1230)
  write ( *,1410)
  call statws(yconst, wtall0, nconst, ldstak, sts, nprton)
  write ( *,1170) ierr
!
!  test 2.  check for reading outside of data array.
!
  write ( *,1160)
  ytemp1 = yconst(1)
  yconst(1) = fplm
  ytempn = yconst(nconst)
  yconst(nconst) = fplm
  write ( *,1440)
  write ( *,1240)
  call stat(yconst(2), nconst-2, ldstak)
  write ( *,1170) ierr
  write ( *,1440)
  write ( *,1250)
  call stats(yconst(2), nconst-2, ldstak, sts, nprton)
  write ( *,1170) ierr
  write ( *,1440)
  write ( *,1400)
  call statw(yconst(2), wt, nconst-2, ldstak)
  write ( *,1170) ierr
  write ( *,1440)
  write ( *,1410)
  call statws(yconst(2), wt, nconst-2, ldstak, sts, nprton)
  write ( *,1170) ierr
  yconst(1) = ytemp1
  yconst(nconst) = ytempn
!
!  test 3.  constant y.
!
  write ( *,1200)
  write ( *,1440)
  write ( *,1240)
  call stat(yconst, nconst, ldstak)
  write ( *,1170) ierr
  write ( *,1440)
  write ( *,1250)
  call stats(yconst, nconst, ldstak, sts, nprton)
  write ( *,1170) ierr
  write ( *,1440)
  write ( *,1400)
  call statw(yconst, wt, nconst, ldstak)
  write ( *,1170) ierr
  write ( *,1440)
  write ( *,1410)
  call statws(yconst, wt, nconst, ldstak, sts, nprton)
  write ( *,1170) ierr
!
!  test 4.  see if turning off the printout works.
!
  write ( *,1260)
  write ( *,1270)
  write ( *,1230)
  write ( *,1250)
  call stats(y, n, ldstak, sts, nprtof)
  write ( *,1390) y(1:10)
  write ( *,1280)
  write ( *,1230)
  write ( *,1410)
  call statws(y, wt, n, ldstak, sts, nprtof)
  write ( *,1390) y(1:10)
!
!     test 5.  make a working run of each routine  first with
!              n=2 (the minimun valid value) and then for the whole
!              data set to check the output.
!
  write ( *,1300)
  write ( *,1310)
  write ( *,1240)
  call stat(y, 3, ldstak)
  write ( *,1310)
  write ( *,1240)
  call stat(y, n, ldstak)
  write ( *,1390) y(1:10)
  write ( *,1170) ierr

  write ( *,1320)
  write ( *,1400)
  call statw(y, wt, 3, ldstak)
  write ( *,1320)
  write ( *,1400)
  call statw(y, wt, n, ldstak)
  write ( *,1390) y(1:10)
  write ( *,1170) ierr

  write ( *,1340)
  write ( *,1250)
  call stats(y, 3, ldstak, sts, nprton)
  write ( *,1340)
  write ( *,1250)
  call stats(y, n, ldstak, sts, nprton)
  write ( *,1390) y(1:10)
  write ( *,1170) ierr

  write ( *,1350)
  write ( *,1410)
  call statws(y, wt, 3, ldstak, sts, nprton)
  write ( *,1350)
  write ( *,1410)
  call statws(y, wt, n, ldstak, sts, nprton)
  write ( *,1390) y(1:10)
  write ( *,1170) ierr
!
!  test 5.  check results of weighting all observations
!  with 1.0e0.  compare with stat execution.
!
  write ( *,1370)
  write ( *,1400)
  call statw(y, wtall1, n, ldstak)
  write ( *,1170) ierr
!
!     test 6.  check results of forcing difference paths through
!              the summation routines, using small, simple data sets.
!
  write ( *,1000)
!
!  run data set 6.1
!
  do i=1,10
     ypath(i) = i
  end do
  write ( *,1010)
  write ( *,1240)
  call stat(ypath, nconst, ldstak)
  write ( *,1170) ierr
  write ( *,1020)
  write ( *,1400)
  call statw(ypath, wtall1, nconst, ldstak)
  write ( *,1170) ierr
!
!  run data set 6.2
!
  do i=1,10
     ypath(i) = -i
  end do
  write ( *,1030)
  write ( *,1240)
  call stat(ypath, nconst, ldstak)
  write ( *,1170) ierr
  write ( *,1040)
  write ( *,1400)
  call statw(ypath, wtall1, nconst, ldstak)
  write ( *,1170) ierr
!
!     run data set 6.3
!
  do i=1,10
     ypath(i) = i-1
  end do
  write ( *,1050)
  write ( *,1240)
  call stat(ypath, nconst, ldstak)
  write ( *,1170) ierr
  write ( *,1060)
  write ( *,1400)
  call statw(ypath, wtall1, nconst, ldstak)
  write ( *,1170) ierr
!
!     run data set 6.4
!
  do i=1,10
     ypath(i) = 1-i
  end do
  write ( *,1070)
  write ( *,1240)
  call stat(ypath, nconst, ldstak)
  write ( *,1170) ierr
  write ( *,1080)
  write ( *,1400)
  call statw(ypath, wtall1, nconst, ldstak)
  write ( *,1170) ierr
!
!     run data set 6.5
!
  do i=1,10
     ypath(i) = i-6
  end do
  write ( *,1090)
  write ( *,1240)
  call stat(ypath, nconst, ldstak)
  write ( *,1170) ierr
  write ( *,1100)
  write ( *,1400)
  call statw(ypath, wtall1, nconst, ldstak)
  write ( *,1170) ierr
!
!     run data set 6.6
!
  do i=1,10
     ypath(i) = i-5
  end do
  write ( *,1110)
  write ( *,1240)
  call stat(ypath, nconst, ldstak)
  write ( *,1170) ierr
  write ( *,1120)
  write ( *,1400)
  call statw(ypath, wtall1, nconst, ldstak)
  write ( *,1170) ierr
!
!     run data set 6.7
!
  ypath(1:10) = 0.0e0
  ypath(1) = -5.0e0
  ypath(10) = 5.0e0
  write ( *,1130)
  write ( *,1240)
  call stat(ypath, nconst, ldstak)
  write ( *,1170) ierr
  write ( *,1140)
  write ( *,1400)
  call statw(ypath, wtall1, nconst, ldstak)
  write ( *,1170) ierr
!
!     run data set 6.8
!
  ypath(1:10) = 0.0e0
  ypath(1) = -5.0e0
  wtall1(1) = 0.0e0
  ypath(10) = 5.0e0
  wtall1(10) = 0.0e0
  write ( *,1380)
  write ( *,1400)
  call statw(ypath, wtall1, nconst, ldstak)
  write ( *,1170) ierr
  return

 1000 format('test 6.  try different paths through the summation code.')
 1010 format('1run stat on 1, ..., 10.')
 1020 format('1run statw on 1, ..., 10.  weights are all 1.')
 1030 format('1run stat on -1, ..., -10.')
 1040 format('1run statw on -1, ..., -10.  weights are all 1.')
 1050 format('1run stat on 0, ..., 9.')
 1060 format('1run statw on 0, ..., 9.  weights are all 1.')
 1070 format('1run stat on 0, ..., -9.')
 1080 format('1run statw on 0, ..., -9.  weights are all 1.')
 1090 format('1stat on -5, ..., 4.')
 1100 format('1run statw on -5, ..., 4.  weights are all 1.')
 1110 format('1run stat on -4, ..., 5.')
 1120 format('1run statw on -4, ..., 5.  weights are all 1.')
 1130 format('1run stat on -1, 8*0, 1.')
 1140 format('1run statw on -1, 8*0, 1.  weights are all 1.')
 1150 format('1test runs for the statistical analysis family routines.')
 1160 format('1test runs to be sure code is not reading outside', &
         ' data array.')
 1170 format(/' the value of ierr is ', i4)
 1180 format('1try two or fewer elements.')
 1190 format('1try insufficient work area.')
 1200 format('1try constant y.')
 1210 format('1try negative weights.')
 1220 format('1try all weights zero (and constant y).')
 1230 format (///)
 1240 format (' call to stat')
 1250 format (' call to stats')
 1260 format(45h1test3.  try turning off the print for those , &
     24hroutines which allow it.)
 1270 format(37h try turning the print off for stats.)
 1280 format(38h try turning the print off for statws.)
 1300 format(52h1test 4.  make working runs of all routines to check, &
     16h the statistics.)
 1310 format('1run stat on the davis-harrison pikes peak data.')
 1320 format('1run statw on the davis-harrison pikes peak data.')
 1340 format('1run stats on the davis-harrison pikes peak data.')
 1350 format('1run statws on the davis-harrison pikes peak data.')
 1370 format('1run statw on the davis-harrison pikes peak data.', &
    '  weights all equal to one.  compare to stat above, not to', &
    ' statw.')
 1380 format('series with nonzero values weighted zero.')
 1390 format(/' data = ', 10f7.4)
 1400 format ('call to statw')
 1410 format (' call to statws')
 1440 format ('1')
end
