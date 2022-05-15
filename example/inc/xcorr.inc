subroutine xcorr ( ldstak ,typ)
use M_starpac, only : corr,corrs,corrxp,genr,iprint,ldscmp,msgx,setra
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XCORR exercises all aspects of the correlation family routines.
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
!  Reference:
!
!    Draper, Smith,
!    Applied Regression Analysis
!    page 216
!
!  Parameters:
!
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer i
!         an index variable.
!     integer ierr
!        common flag indicating if any errors were detected
!        if ierr = 0, then no errors were found
!     integer ivcv
!        the row dimension of vcv
!     integer iym
!        the row dimension of ym
!     integer j
!        an index variable.
!     integer ldsmin
!        the smallest acceptable size of common area cstak
!     integer ldstak
!        the size of the common area cstak
!     integer m
!        the number of variables
!     integer n
!        the number of observations
!     real vcv(4,4)
!        the variance covariance matrix
!     real ym(10,4)
!        general data set, from draper and smith
!     real z(10,4)
!        test observation matrix
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
     i,ivcv,iym,j,ldsmin,m,n
!
!  local arrays
  real &
     vcv(4,4),ym(10,4),z(10,4)
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

  ivcv = 4
  iym = 10
  m = 4
  n = 9
  ierr = 0
!
!  test routines with correct call statement.
!
  write ( *,1000)
  write ( *,1010)
!
!  test corr
!
  write ( *,1020)
  write ( *,1060)
  call corr(ym, n, m, iym, ldstak)
  call msgx(0 )
!
!  test corrs
!
!  printout suppressed
!
  write ( *,1030)
  write ( *,1040)
  write ( *,1060)
  call corrs(ym, n, m, iym, ldstak, 0, vcv, ivcv)
  call msgx ( 0 )
!
!  print stored output and zero arrays
!
  call corrxp ( m, vcv, ivcv )
!
!  with printout
!
  write ( *,1050)
  write ( *,1060)
  call corrs(ym, n, m, iym, ldstak, 1, vcv, ivcv)
  call msgx ( 0 )
!
!  print stored output
!
  call corrxp ( m, vcv, ivcv )
!
!  special 2 column matrix.
!
  write ( *,1070)
  write ( *,1060)
  call corr(ym, n, 2, iym, ldstak)
  call msgx ( 0 )
!
!  test work area requirements.
!
!  test corr
!
  call ldscmp(12, 0, max(n,m), 0, 0, 0, 's', &
     m*m + (max(n,m)+m+n*(m+3)+6*m*m), ldsmin)
  write ( *,1090)
  call corr(ym, n, m, iym, ldsmin-1)
  call msgx ( 1 )
  write ( *,1100)
  call corr(ym, n, m, iym, ldsmin)
  call msgx ( 0 )
!
!  test corrs with printout
!
  call ldscmp(12, 0, max(n,m), 0, 0, 0, 's', &
     max(n,m)+m+n*(m+3)+6*m*m, ldsmin)
  write ( *,1090)
  call corrs(ym, n, m, iym, ldsmin-1, 1, vcv, ivcv)
  call msgx ( 1 )
  write ( *,1100)
  call corrs(ym, n, m, iym, ldsmin, 1, vcv, ivcv)
  call corrxp(m, vcv, ivcv )
  call msgx ( 0 )
!
!  test corrs without printout
!
  call ldscmp(12, 0, 0, 0, 0, 0, 's', 0, ldsmin)
  write ( *,1090)
  call corrs(ym, n, m, iym, ldsmin-1, 0, vcv, ivcv)
  call msgx ( 1 )
  write ( *,1100)
  call corrs(ym, n, m, iym, ldsmin, 0, vcv, ivcv)
  call corrxp(m, vcv, ivcv )
  call msgx ( 0 )
!
!  number of variables less than 2.
!
  write ( *,1110)
!
!  test corr
!
  call corr(ym, n, 1, iym, ldstak)
  call msgx ( 1 )
!
!  test corrs
!
  call corrs(ym, n, 1, iym, ldstak, 1, vcv, ivcv)
  call msgx ( 1 )
!
!  number of observations less than 3.
!
  write ( *,1120)
!
!  test corr
!
  call corr(ym, 2, 4, iym, ldstak)
  call msgx ( 1 )
!
!  test corrs
!
  call corrs(ym, 2, 4, iym, ldstak, 1, vcv, ivcv)
  call msgx ( 1 )
!
!  observation matrix dimensioned less than n.
!
  write ( *,1150)
!
!  test corr
!
  call corr(ym, n, m, 8, ldstak)
  call msgx ( 1 )
!
!  test corrs
!
  call corrs(ym, n, m, 8, ldstak, 1, vcv, ivcv)
  call msgx ( 1 )
!
!  vcv matrix dimensioned less than m.
!
  write ( *,1130)
  call corrs(ym, n, m, iym, ldstak, 1, vcv, 2)
  call msgx ( 1 )
!
!  all observations on a single variable equal to zero.
!
  write ( *,1140)
  call setra(z, 10, 4, 10, 0.0e0)
  call corr(z, 9, 4, 10, ldstak)
  call msgx ( 1 )
  call corrs(z, 9, 4, 10, ldstak, 1, vcv, ivcv)
  call corrxp(m, vcv, ivcv )
  call msgx ( 1 )

  z(1:10,1) = i
  z(1:10,2) = 0.0e0
  call corr(z, 10, 4, 10, ldstak)
  call msgx ( 1 )
!
!  array filled with a single value.
!
  write ( *,1160)
  call setra(z, 10, 4, 10, 4.0e0)
  call corr(z, 4, 10, 4, ldstak)
  call msgx ( 1 )
!
!  2 columns the same.
!
  do i=1,3
     call genr(z(1,i), 5, 5.0e0*i, 5.0e0*i)
  end do
  z(1:5,4) = z(1:5,3)

  write ( *,1170)
  call corr(z, 5, 4, 10, ldstak)
  call msgx ( 1 )
!
!  2 columns inversely related.
!
  j = 5
  do i=1,5
     j = j - 1
     z(j,4) = z(i,3)
  end do
  write ( *,1170)
  call corr(z, 5, 4, 10, ldstak)
  call msgx ( 1 )

  return

 1000 format('1')
 1010 format(' ****test routines with correct call****')
 1020 format(' test of corr')
 1030 format('1test of corrs')
 1040 format(' printout supressed.')
 1050 format('1printout not supressed.')
 1060 format(' draper and smith data set (page 216).')
 1070 format('1****special case 2 column matrix****')
 1090 format('1****test with insufficient work area****')
 1100 format('1****test with exactly the right amount of work area****')
 1110 format('1****number of variables less than 2****')
 1120 format(' ****number of observations less than 3****')
 1130 format(' ****inadequate space in storage arrays****')
 1140 format('1****all observations on a variable equal to zero****')
 1150 format(' ****observation matrix dimensioned less than number', &
         ' of observations designated****')
 1160 format('1****array containing a single value****')
 1170 format('1****2 columns related****')
end
