subroutine xaov1 ( ldstak ,typ)
use M_starpac, only : aov1,aov1s,aov1xp,iprint,ldscmp,msgx,setrv
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XAOV1 exercises the oneway family routines.
!
!  Modified:
!
!    24 April 2006
!
!  Author:
!
!    Linda Mitchell,
!    Statistical Engineering Division,
!    National Bureau of Standards,
!    Boulder, Colorado
!
!  Parameters:
!
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     real fi
!        float of index i
!     real gstat(10,4)
!        the group statistics.  columns correspond to the tag
!        value, sample size, group mean, and group standard deviation.
!     integer i
!        index variable
!     integer ierr
!        common flag indicating whether or not there were any errors
!     integer igstat
!        the first dimension of gstat.
!     integer ldsmin
!        the smallest acceptable size of the common cstak
!     integer ldstak
!        the size of the common cstak
!     integer n
!        the number of observations
!     integer ng
!        the number of different groups
!     real tag(20)
!        the tag values for each observation
!     real y(20)
!        the vector of observations
!     real z(10)
!        test vector
!     real ztag(10)
!        test tag vector
!

  integer ierr
  integer ldstak
!
!  arrays in common
  double precision dstak(12)
!
!  local scalars
  real &
     fi
  integer &
     i,igstat,ldsmin,n,ng
!
!  local arrays
  real &
     gstat(10,4),tag(20),y(20),z(10),ztag(10)
!
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr

  data   y( 1),   y( 2),   y( 3),   y( 4),   y( 5) &
     /    61.0e0,    61.0e0,    67.0e0,    67.0e0,    64.0e0/
  data   y( 6),   y( 7),   y( 8),   y( 9),   y(10) &
     /    78.0e0,    71.0e0,    75.0e0,    72.0e0,    74.0e0/
  data   y(11),   y(12),   y(13),   y(14),   y(15) &
     /    83.0e0,    81.0e0,    76.0e0,    78.0e0,    79.0e0/
  data   y(16),   y(17) &
     /    72.0e0,   72.0e0/

  data tag( 1), tag( 2), tag( 3), tag( 4), tag( 5) &
     /    11.5e0,    11.5e0,    11.5e0,    11.5e0,    11.5e0/
  data tag( 6), tag( 7), tag( 8), tag( 9), tag(10) &
     /    12.0e0,    12.0e0,    12.0e0,    12.0e0,    12.0e0/
  data tag(11), tag(12), tag(13), tag(14), tag(15) &
     /    11.0e0,    11.0e0,    11.0e0,    11.0e0,    11.0e0/
  data tag(16), tag(17) &
     /   -11.0e0,   11.0e0/
!
!  set various dimensions and program variables
!
  n = 17
  igstat = 10
!
! test with correct call statements.
!
  write ( *,1060)
!
!  test aov1
!
  write ( *,1000)
  call aov1( y, tag, n, ldstak )
  call msgx ( 0 )
  write ( *,1010)
  write ( *,1020) (tag(i),i=1,n)
!
!  test of aov1s
!
!  printout not supressed
!
  write ( *,1040)
  write ( *,1030)
  call aov1s(y, tag, n, ldstak, 1, gstat, igstat, ng)
  call msgx(0 )
!
!  print storage and zero vectors before using again
!
  call aov1xp(gstat, igstat, ng)
!
!  printout supressed
!
  write ( *,1050)
  call aov1s(y, tag, n, ldstak, 0, gstat, igstat, ng)
  call msgx(0 )
!
!  print storage and zero vectors before using again
!
  call aov1xp(gstat, igstat, ng)
!
!  number of observations less than 2
!
  write ( *,1090)
  call aov1(y, tag, 1, ldstak)
  call msgx(1 )

  call aov1s(y, tag, -14, ldstak, 0, gstat, igstat, ng)
  call msgx(1 )
!
!  all observations the same value
!
  write ( *,1100)
  write ( *,1000)
  call setrv(z, 10, 0.0e0)
  call aov1(z, tag, 10, ldstak)
  call msgx(0 )

  call setrv(z, 10, 2.0e0)
  write ( *,1100)
  write ( *,1030)
  call aov1s(z, tag, 10, ldstak, 1, gstat, igstat, ng)
  call msgx(0 )
!
!  test work area size handling
!
  call ldscmp(11, 0, 33, 0, 0, 0, 's', 40, ldsmin)
  write ( *,1070)
  write ( *,1000)
  call aov1(y, tag, n, 1)
  call msgx(1 )
  write ( *,1070)
  write ( *,1000)
  call aov1(y, tag, n, ldsmin-1)
  call msgx(1 )
  write ( *,1010)
  write ( *,1020) (tag(i),i=1,n)
  write ( *,1080)
  write ( *,1000)
  call aov1(y, tag, n, ldsmin)
  call msgx(0 )

  call ldscmp(11, 0, 33, 0, 0, 0, 's', 28, ldsmin)
  write ( *,1070)
  write ( *,1030)
  call aov1s(y, tag, n, 1, 0, gstat, igstat, ng)
  call msgx(1 )
  write ( *,1070)
  write ( *,1030)
  call aov1s(y, tag, n, ldsmin-1, 0, gstat, igstat, ng)
  call msgx(1 )
  write ( *,1010)
  write ( *,1020) (tag(i),i=1,n)
  write ( *,1080)
  write ( *,1030)
  call aov1s(y, tag, n, ldsmin, 1, gstat, igstat, ng)
  call msgx(0 )
!
!  same number of groups as non-zero tags
!
  write ( *,1120)
  do i=1,10
     fi = i
     ztag(i) = fi
     z(i) = 13.0e0 - fi
  end do
  call aov1(z, ztag, 10, ldstak)
  call msgx(1 )
  call aov1s(z, ztag, 10, ldstak, 1, gstat, igstat, ng)
  call msgx(1 )
!
!  less than 2 different tag groups
!
  write ( *,1130)
  call setrv(ztag, 10, 1.0e0)
  call aov1(z, ztag, 10, ldstak)
  call msgx(1 )
  call aov1s(z, ztag, 10, ldstak, 1, gstat, igstat, ng)
  call msgx(1 )
!
!  less than 2 tags
!
  call setrv(ztag, 9, 0.0e0)
  write ( *,1140)
  call aov1(z, ztag, 10, ldstak)
  call msgx(1 )
  call aov1s(z, ztag, 10, ldstak, 1, gstat, igstat, ng)
  call msgx(1 )
!
!  incorrect dimension of gstat
!
  write ( *,1150)
  call aov1s(y, tag, n, ldstak, 1, gstat, 2, ng)
  call msgx(1 )
!
!  all observations within a group same value
!
  z(1) = 53.0e0
  ztag(1) = 1.0e0
  z(2) = 62.0e0
  ztag(2) = 3.0e0
  z(3) = 53.0e0
  ztag(3) = 1.0e0
  z(4) = 71.0e0
  ztag(4) = 4.0e0
  z(5) = 89.0e0
  ztag(5) = 2.0e0
  z(6) = 71.0e0
  ztag(6) = 4.0e0
  z(7) = 89.0e0
  ztag(7) = 2.0e0
  z(8) = 62.0e0
  ztag(8) = 3.0e0
  z(9) = 71.0e0
  ztag(9) = 4.0e0
  z(10) = 62.0e0
  ztag(10) = 3.0e0
  write ( *,1160)
  call aov1(z, ztag, 10, ldstak)
  call msgx(0 )
!
!  2 tags
!
  write ( *,1170)
  call aov1(z, ztag, 3, ldstak)
  call msgx(0 )
!
!  all groups(except for 1) with 1 observation
!
  write ( *,1180)
  call aov1(z, ztag, 5, ldstak)
  call msgx(0 )

  return

 1000 format(' test of aov1 ')
 1010 format(' check to see if tags have been changed')
 1020 format(4f12.6)
 1030 format(' test of aov1s ')
 1040 format('1printout not supressed.')
 1050 format(' printout supressed.')
 1060 format('1****test routines with correct call****')
 1070 format('1****test with insufficient work area****')
 1080 format('1****test with exactly the right amount of work area****')
 1090 format('1****number of observations less than 2****')
 1100 format('1****all observations with same value****')
 1120 format('1****same number of groups as non-zero tags****')
 1130 format(' ****less than 2 different tag groups****')
 1140 format(' ****less than 2 tags****')
 1150 format('1****incorrect dimension of gstat****')
 1160 format('1****all observations within a group same value****')
 1170 format('1****test with 2 tags****')
 1180 format('1****all groups except for 1 with 1 observation ****')
end
