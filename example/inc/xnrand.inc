subroutine xnrand ( ldstak ,typ)
use M_starpac, only : hist,iprint,nrand,nrandc
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XNRAND tests the NRAND family.
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
!     integer i
!        an index variable
!     integer ierr
!        flag to indicate presence of error detected by preceding
!        starpac call.  (0 is ok, 1 is error)
!     integer iseed
!        the seed for the random number generator.
!     integer ldstak
!        amount of work area.  size of dstak.
!     integer n
!        the length of the vector y.
!     real sigma
!        the s.d. of the sample.
!     real y(1000)
!        data vector for tests.
!     real ymean
!        the mean of the sample.
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
     sigma,ymean
  integer &
     i,iseed,n
!
!  local arrays
  real &
     y(1000)
!
!
!  common blocks
  common /cstak/dstak
  common /errchk/ierr
!
!  check for sufficient work area length.
!
  if (ldstak.lt.1000) then
    write ( *, 1000)
     return
  end if
!
!  heading.
!
  write ( *,1150)
!
!  test 1.  check all error messages.
!
  write ( *,1160)
!
!  error 1, zero or fewer elements or negative standard deviation.
!
  iseed = 0
  sigma = -1
  n = 0
  ymean = 0.0e0
  call nrand(y, n, iseed)
  write ( *,1170) ierr
  call nrandc(y, n, iseed, ymean, sigma)
  write ( *,1170) ierr
!
!  compare results
!
  iseed = 334
  n = 10
  ymean = 0.0e0
  sigma = 1.0e0

  write ( *, 1120) n, iseed
  call nrand (y, n, iseed)
  write ( *, 1100) (y(i),i=1,n)

  iseed = 333
  write ( *, 1130) n, ymean, sigma, iseed
  call nrandc (y, n, iseed, ymean, sigma)
  write ( *, 1100) (y(i),i=1,n)

  iseed = 13531
  n = 1000
  ymean = 0.0e0
  sigma = 1.0e0

  write ( *, 1120) n, iseed
  call nrand (y, n, iseed)
  call hist (y, n, ldstak)

  write ( *, 1130) n, ymean, sigma, iseed
  call nrandc (y, n, iseed, ymean, sigma)
  call hist (y, n, ldstak)

  iseed = 99999
  n = 1000
  ymean = 4.0e0
  sigma = 4.0e0

  write ( *, 1120) n, iseed
  call nrand (y, n, iseed)
  call hist (y, n, ldstak)

  write ( *, 1130) n, ymean, sigma, iseed
  call nrandc (y, n, iseed, ymean, sigma)
  call hist (y, n, ldstak)

  return

 1000 format ('1the dimension of dstak and the value of ldstak needed'/ &
    ' for nrandx must equal or exceed 1000.  change driver'/ &
    ' and recall nrandx.')
 1100 format (5e15.8)
 1120 format ('1generate ', i4, &
    ' standard normal numbers using iseed = ', i5)
 1130 format ('1generate ', i4, &
    ' normally distributed numbers with ymean = ', f5.2, &
    ' and sigma = ', f5.2, ' using iseed = ', i5)
 1150 format ('1test runs for the nrand family of routines.')
 1160 format(' test 1.  generate each of the possible ', &
     15herror messages.)
 1170 format(/22h the value of ierr is , i4//)
end
