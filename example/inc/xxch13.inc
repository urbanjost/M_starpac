subroutine xxch13 ( ldstak ,typ)
use M_starpac, only : aime,aimf,iprint
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XXCH13 tests the arima modeling and forecasting family of routines.
!
!  Discussion:
!
!    Data is the airline data listed on page 531 of box and jenkins.
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
!    George Box, Gwilym Jenkins,
!    Time Series Analysis: Forecasting and Control,
!    Holden-Day, 1970,
!    QA280.B67.
!
!  Parameters:
!
!     real air(200)
!        the airline data.
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     integer i
!        an index variable.
!     integer ierr
!        the integer value returned by this routine designating
!        whether any errors were detected in the parameter list.
!        if ierr == 0, no errors were detected.
!        if ierr .ge. 1, errors were detected.
!     integer ldstak
!        the length of dstak in common /cstak/.
!     integer mspec(4,10)
!        the array containing the values of p, d, q, and s for each
!        factor.
!     integer n
!        the number of observations.
!     integer nfac
!        the number of factors in the model
!     integer npar
!        the number of unknown parameters in the model.
!     real par(10)
!        the array in which the current estimates of the unknown
!        parameters are stored.
!     real res(200)
!        the residuals from the fit.
!     real y(200)
!        the array of the dependent variable.
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
      integer n,nfac,npar
!
!  local arrays
      real &
         air(200),par(10),res(200),y(200)
      integer &
         mspec(4,10)
!
!
!
!  common blocks
      common /cstak/dstak
      common /errchk/ierr
!
!  define constants
!
      data mspec(1,1), mspec(2,1), mspec(3,1), mspec(4,1) &
         /          0,          1,          1,          1/
      data mspec(1,2), mspec(2,2), mspec(3,2), mspec(4,2) &
         /          0,          1,          1,         12/

      data  air(  1), air(  2), air(  3), air(  4), air(  5), air(  6) &
          / 112.0e0, 118.0e0, 132.0e0, 129.0e0, 121.0e0, 135.0e0/
      data  air(  7), air(  8), air(  9), air( 10), air( 11), air( 12) &
          / 148.0e0, 148.0e0, 136.0e0, 119.0e0, 104.0e0, 118.0e0/
      data  air( 13), air( 14), air( 15), air( 16), air( 17), air( 18) &
          / 115.0e0, 126.0e0, 141.0e0, 135.0e0, 125.0e0, 149.0e0/
      data  air( 19), air( 20), air( 21), air( 22), air( 23), air( 24) &
          / 170.0e0, 170.0e0, 158.0e0, 133.0e0, 114.0e0, 140.0e0/
      data  air( 25), air( 26), air( 27), air( 28), air( 29), air( 30) &
          / 145.0e0, 150.0e0, 178.0e0, 163.0e0, 172.0e0, 178.0e0/
      data  air( 31), air( 32), air( 33), air( 34), air( 35), air( 36) &
          / 199.0e0, 199.0e0, 184.0e0, 162.0e0, 146.0e0, 166.0e0/
      data  air( 37), air( 38), air( 39), air( 40), air( 41), air( 42) &
          / 171.0e0, 180.0e0, 193.0e0, 181.0e0, 183.0e0, 218.0e0/
      data  air( 43), air( 44), air( 45), air( 46), air( 47), air( 48) &
          / 230.0e0, 242.0e0, 209.0e0, 191.0e0, 172.0e0, 194.0e0/
      data  air( 49), air( 50), air( 51), air( 52), air( 53), air( 54) &
          / 196.0e0, 196.0e0, 236.0e0, 235.0e0, 229.0e0, 243.0e0/
      data  air( 55), air( 56), air( 57), air( 58), air( 59), air( 60) &
          / 264.0e0, 272.0e0, 237.0e0, 211.0e0, 180.0e0, 201.0e0/
      data  air( 61), air( 62), air( 63), air( 64), air( 65), air( 66) &
          / 204.0e0, 188.0e0, 235.0e0, 227.0e0, 234.0e0, 264.0e0/
      data  air( 67), air( 68), air( 69), air( 70), air( 71), air( 72) &
          / 302.0e0, 293.0e0, 259.0e0, 229.0e0, 203.0e0, 229.0e0/
      data  air( 73), air( 74), air( 75), air( 76), air( 77), air( 78) &
          / 242.0e0, 233.0e0, 267.0e0, 269.0e0, 270.0e0, 315.0e0/
      data  air( 79), air( 80), air( 81), air( 82), air( 83), air( 84) &
          / 364.0e0, 347.0e0, 312.0e0, 274.0e0, 237.0e0, 278.0e0/
      data  air( 85), air( 86), air( 87), air( 88), air( 89), air( 90) &
          / 284.0e0, 277.0e0, 317.0e0, 313.0e0, 318.0e0, 374.0e0/
      data  air( 91), air( 92), air( 93), air( 94), air( 95), air( 96) &
          / 413.0e0, 405.0e0, 355.0e0, 306.0e0, 271.0e0, 306.0e0/
      data  air( 97), air( 98), air( 99), air(100), air(101), air(102) &
          / 315.0e0, 301.0e0, 356.0e0, 348.0e0, 355.0e0, 422.0e0/
      data  air(103), air(104), air(105), air(106), air(107), air(108) &
          / 465.0e0, 467.0e0, 404.0e0, 347.0e0, 305.0e0, 336.0e0/
      data  air(109), air(110), air(111), air(112), air(113), air(114) &
          / 340.0e0, 318.0e0, 362.0e0, 348.0e0, 363.0e0, 435.0e0/
      data  air(115), air(116), air(117), air(118), air(119), air(120) &
          / 491.0e0, 505.0e0, 404.0e0, 359.0e0, 310.0e0, 337.0e0/
      data  air(121), air(122), air(123), air(124), air(125), air(126) &
          / 360.0e0, 342.0e0, 406.0e0, 396.0e0, 420.0e0, 472.0e0/
      data  air(127), air(128), air(129), air(130), air(131), air(132) &
          / 548.0e0, 559.0e0, 463.0e0, 407.0e0, 362.0e0, 405.0e0/
      data  air(133), air(134), air(135), air(136), air(137), air(138) &
          / 417.0e0, 391.0e0, 419.0e0, 461.0e0, 472.0e0, 535.0e0/
      data  air(139), air(140), air(141), air(142), air(143), air(144) &
          / 622.0e0, 606.0e0, 508.0e0, 461.0e0, 390.0e0, 432.0e0/

      nfac = 2
      n = 144

      npar = 3
      par(1) = 0.000
      par(2) = 0.395
      par(3) = 0.615

      y(1:n) = log(air(1:n))
!
!  run simple test of aime
!
      write ( *,1000)
      write ( *,1100)
      call aime (y, n, mspec, nfac, par, npar, res, ldstak)
      write ( *,2000) ierr
!
!  run simple test of aimf
!
      write ( *,1200)
      call aimf (y, n, mspec, nfac, par, npar, ldstak)
      write ( *,2000) ierr

      return
 1000 format ('1*ch13')
 1100 format (' simple test of aime')
 1200 format ('1simple test of aimf')
 2000 format (/' the value of ierr is ', i4)
      end subroutine xxch13
