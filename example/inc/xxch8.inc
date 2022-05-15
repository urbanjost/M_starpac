subroutine xxch8 ( ldstak ,typ)
use M_starpac, only : iprint,lls,llsp
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XXCH8 tests the linear least squares family of routines.
!
!  Discussion:
!
!    LLS problem is from daniel and wood [1971], pages 61-65.
!
!    LLSP problem is from miller and freund [1977], page 311.
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
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer ierr
!        the integer value designating whether any errors were
!        detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr == 1, errors were detected.
!     integer ixm
!        the first dimension of the matrix x.
!     integer ldstak
!        the length of dstak in common /cstak/.
!     integer ndeg
!        the degree of the polynomial model to be fit.
!     integer npar
!        the number of parameters to be estimated.
!     integer n1, n2
!        the number of observations in each problem.
!     real res(25)
!        the residuals.
!     real x(25)
!        the independent variable.
!     real xm(25,5)
!        the independent variable.
!     real y1(25), y2(25)
!        the dependent variable.
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
      integer ixm,n1,n2,ndeg,npar
!
!  local arrays
      real &
         res(25),x(25),xm(25,5),y1(25),y2(25)
!
!
!  common blocks
      common /cstak/dstak
      common /errchk/ierr

      data      xm(1,1),  xm(1,2),  xm(1,3),  xm(1,4) &
          /      1.0e0, 80.0e0, 27.0e0, 89.0e0/
      data      xm(2,1),  xm(2,2),  xm(2,3),  xm(2,4) &
          /      1.0e0, 80.0e0, 27.0e0, 88.0e0/
      data      xm(3,1),  xm(3,2),  xm(3,3),  xm(3,4) &
          /      1.0e0, 75.0e0, 25.0e0, 90.0e0/
      data      xm(4,1),  xm(4,2),  xm(4,3),  xm(4,4) &
          /      1.0e0, 62.0e0, 24.0e0, 87.0e0/
      data      xm(5,1),  xm(5,2),  xm(5,3),  xm(5,4) &
          /      1.0e0, 62.0e0, 22.0e0, 87.0e0/
      data      xm(6,1),  xm(6,2),  xm(6,3),  xm(6,4) &
          /      1.0e0, 62.0e0, 23.0e0, 87.0e0/
      data      xm(7,1),  xm(7,2),  xm(7,3),  xm(7,4) &
          /      1.0e0, 62.0e0, 24.0e0, 93.0e0/
      data      xm(8,1),  xm(8,2),  xm(8,3),  xm(8,4) &
          /      1.0e0, 62.0e0, 24.0e0, 93.0e0/
      data      xm(9,1),  xm(9,2),  xm(9,3),  xm(9,4) &
          /      1.0e0, 58.0e0, 23.0e0, 87.0e0/
      data     xm(10,1), xm(10,2), xm(10,3), xm(10,4) &
          /      1.0e0, 58.0e0, 18.0e0, 80.0e0/
      data     xm(11,1), xm(11,2), xm(11,3), xm(11,4) &
          /      1.0e0, 58.0e0, 18.0e0, 89.0e0/
      data     xm(12,1), xm(12,2), xm(12,3), xm(12,4) &
          /      1.0e0, 58.0e0, 17.0e0, 88.0e0/
      data     xm(13,1), xm(13,2), xm(13,3), xm(13,4) &
          /      1.0e0, 58.0e0, 18.0e0, 82.0e0/
      data     xm(14,1), xm(14,2), xm(14,3), xm(14,4) &
          /      1.0e0, 58.0e0, 19.0e0, 93.0e0/
      data     xm(15,1), xm(15,2), xm(15,3), xm(15,4) &
          /      1.0e0, 50.0e0, 18.0e0, 89.0e0/
      data     xm(16,1), xm(16,2), xm(16,3), xm(16,4) &
          /      1.0e0, 50.0e0, 18.0e0, 86.0e0/
      data     xm(17,1), xm(17,2), xm(17,3), xm(17,4) &
          /      1.0e0, 50.0e0, 19.0e0, 72.0e0/
      data     xm(18,1), xm(18,2), xm(18,3), xm(18,4) &
          /      1.0e0, 50.0e0, 19.0e0, 79.0e0/
      data     xm(19,1), xm(19,2), xm(19,3), xm(19,4) &
          /      1.0e0, 50.0e0, 20.0e0, 80.0e0/
      data     xm(20,1), xm(20,2), xm(20,3), xm(20,4) &
          /      1.0e0, 56.0e0, 20.0e0, 82.0e0/
      data     xm(21,1), xm(21,2), xm(21,3), xm(21,4) &
          /      1.0e0, 70.0e0, 20.0e0, 91.0e0/
!
      data        y1(1),    y1(2),    y1(3) &
          /     42.0e0, 37.0e0, 37.0e0/
      data        y1(4),    y1(5),    y1(6) &
          /     28.0e0, 18.0e0, 18.0e0/
      data        y1(7),    y1(8),    y1(9) &
          /     19.0e0, 20.0e0, 15.0e0/
      data       y1(10),   y1(11),   y1(12) &
          /     14.0e0, 14.0e0, 13.0e0/
      data       y1(13),   y1(14),   y1(15) &
          /     11.0e0, 12.0e0,  8.0e0/
      data       y1(16),   y1(17),   y1(18) &
          /      7.0e0,  8.0e0,  8.0e0/
      data       y1(19),   y1(20),   y1(21) &
          /      9.0e0, 15.0e0, 15.0e0/
!
      data         x(1),     x(2),     x(3) &
          /      0.0e0,  1.0e0,  2.0e0/
      data         x(4),     x(5),     x(6) &
          /      3.0e0,  4.0e0,  5.0e0/
      data         x(7),     x(8),     x(9) &
          /      6.0e0,  7.0e0,  8.0e0/

      data        y2(1),    y2(2),    y2(3) &
          /     12.0e0, 10.5e0, 10.0e0/
      data        y2(4),    y2(5),    y2(6) &
          /      8.0e0,  7.0e0,  8.0e0/
      data        y2(7),    y2(8),    y2(9) &
          /      7.5e0,  8.5e0,  9.0e0/
!
!     set parameters necessary for the computations
!
      ixm = 25
      n1 = 21
      n2 = 9
      npar = 4
      ndeg = 2
!
!     print header
!
      write ( *,1000)
!
!     run simple example of lls
!
      write ( *,1100)
      call lls(y1, xm, n1, ixm, npar, res, ldstak)
      write ( *,2000) ierr
!
!     run simple example of llsp
!
      write ( *,1200)
      call llsp(y2, x, n2, ndeg, res, ldstak)
      write ( *,2000) ierr

      return

 1000 format ('1*ch8')
 1100 format (' simple test of lls')
 1200 format ('1simple test of llsp')
 2000 format (/' the value of ierr is ', i4)
end
