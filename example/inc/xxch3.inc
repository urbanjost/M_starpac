subroutine xxch3 ( typ)
use M_starpac, only : iprint,mvp,nrand,nrandc
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XXCH3 tests the normal random number generator family of routines.
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
!     integer i
!        an index variable.
!     integer ierr
!        flag to indicate presence of error detected by preceding
!        starpac call.  (0 is ok, 1 is error)
!     integer iseed
!        the seed for the random number generator.
!     integer iym
!        the exact value of the first dimension of array ym.
!     integer ldstak
!        a dummy variable for this test subprogram.
!     integer m
!        the number of sets of numbers to be generated
!     integer n
!        the number of observations to be generated.
!     real sigma
!        the s.d. of the sample.
!     real ym(50,2)
!        data vector for tests.
!     real ymean
!        the mean of the sample.
!

  integer ierr
!
!  local scalars
      real &
         sigma,ymean
      integer &
         i,iseed,iym,m,n
!
!  local arrays
      real &
         ym(50,2)
!
!
!  common blocks
      common /errchk/ierr

!
!     data initialization
!
      iym = 50
      iseed = 531
      n = 50
      m = 2
      ymean = 4.0e0
      sigma = 0.5e0
!
!
!  write heading
!
      write ( *,1000)
!
!     generate standard normal pseudo-random numbers into column 1 of ym
!
      write ( *,1100)
      call nrand(ym(1,1), n, iseed)
      write ( *,2000) ierr
      write ( *, 1400) (ym(i,1),i=1,n)
!
!     generate normal pseudo-random numbers
!     with mean 4.0 and standard deviation 0.5 into column 2 of ym
!
      write ( *,1200)
      call nrandc(ym(1,2), n, iseed, ymean, sigma)
      write ( *,2000) ierr
      write ( *, 1400) (ym(i,2),i=1,n)
!
!     plot results, sampling every observation
!
      write ( *,1500)
      call mvp (ym, n, m, iym, 1)
!
      return

 1000 format ('1*ch3')
 1100 format (' simple test of nrand')
 1200 format ('1simple test of nrandc')
 2000 format (/' the value of ierr is ', i4)
 1400 format (/' generated results = '//(5e15.8))
 1500 format ('1mvp display of generated results')
end
