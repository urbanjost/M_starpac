subroutine xxch12 ( ldstak ,typ)
use M_starpac, only : acf,bfs,ccf,iprint,mdflt,pgms,ppl,taper,uas,ufs
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XXCH12 tests time series correlation and spectrum analysis routines.
!
!  Discussion:
!
!    data for acf is taken from p. 362 of jenkins and watts [1968]
!
!    data for ccf is taken from p. 361 of jenkins and watts [1968]
!
!    data for ufs is taken from p. 318 of jenkins and watts [1968]
!
!    data for uas is taken from p. 318 of jenkins and watts [1968]
!
!    data for taper, pgms, mdflt and ppl is
!    the wolf sunspot numbers for the years 1700 to 1960 as
!    tabulated by waldmeier [1961].
!
!    data for bfs is taken from pp. 387-388 of jenkins and watts [1968]
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
!    Gwilym Jenkins, Donald Watts,
!    Spectral Analysis and its Applications,
!    Holden-Day 1968.
!
!    Max Waldmeier,
!    The Sunspot-Activity in the Years 1610-1960,
!    Shulthess, Zurich, 1961.
!
!  Parameters:
!
!     double precision dstak(12)
!        the double precision version of the /cstak/ work area.
!     real freq(300)
!        the frequencies at which the periodogram is computed.
!     integer ierr
!        a common variable used as a flag to indicate whether
!        or not there are any errors, if =0 then no errors.
!     integer iextnd
!        the indicator variable used to designate whether zero or the
!        series mean is to be used to extend the series.
!     integer ilog
!        the indicator variable used to designate whether the plot is
!        to have logarithmic axis or not.
!     integer kmd(10)
!        the vector of modified daniel filter lengths.
!     integer ldstak
!        the length of dstak in common /cstak/.
!     integer lfreq
!        the length of vector freq.
!     integer lper
!        the length of vector per.
!     integer lyfft
!        the length of vector yfft.
!     integer nf
!        the number of frequencies.
!     integer nfft
!        the extended series length for the fft.
!     integer nk
!        the number of daniel filters to apply.
!     integer nprt
!        the print control variable.
!     integer ny1, ny2, ny3, ny4, ny5, ny6
!        the number of observations.
!     real per(300)
!        the periodogram.
!     real perf(300)
!        the filtered periodogram.
!     real taperp
!        the percentage of the series to be tapered.
!     real yfft(600)
!        an array for the fft computations.
!     real y1(100)
!        the data from page 362 of jenkins and watts.
!     real y2a(100), y2b(100)
!        the data from page 361 of jenkins and watts.
!     real y3(50), y4(50)
!        the data from page 318 of jenkins and watts.
!     real y5(300)
!        the wolf sunspot data.
!     real y6a(100), y6b(100)
!        the data from page 387 and 388 of jenkins and watts.
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
      real &
         taperp
      integer &
         iextnd,ilog,lfreq,lper,lyfft,nf,nfft,nk,nprt,ny1,ny2, &
         ny3,ny4,ny5,ny6
!
!  local arrays
      real &
         freq(300),per(300),perf(300),y1(100),y2a(100),y2b(100), &
         y3(50),y4(50),y5(300),y6a(100),y6b(100),yfft(600)
      integer &
         kmd(10)
!
!
!  common blocks
      common /cstak/dstak
      common /errchk/ierr
!
!  equivalences
      equivalence (y3(1),y4(1))
      equivalence (y2a(1),y6a(1))
      equivalence (y2b(1),y6b(1))

      data   y1(  1),  y1(  2),  y1(  3),  y1(  4),  y1(  5),  y1(  6) &
          / -2.07e0, -1.15e0,  0.69e0, -0.46e0, -1.49e0, -0.70e0/
      data   y1(  7),  y1(  8),  y1(  9),  y1( 10),  y1( 11),  y1( 12) &
          / -1.07e0, -0.69e0, -0.68e0,  1.27e0, -1.05e0, -0.05e0/
      data   y1( 13),  y1( 14),  y1( 15),  y1( 16),  y1( 17),  y1( 18) &
          / -0.84e0, -0.62e0, -0.49e0, -1.29e0, -0.49e0, -1.06e0/
      data   y1( 19),  y1( 20),  y1( 21),  y1( 22),  y1( 23),  y1( 24) &
          / -0.38e0, -0.52e0, -0.13e0,  1.30e0, -1.51e0, -0.43e0/
      data   y1( 25),  y1( 26),  y1( 27),  y1( 28),  y1( 29),  y1( 30) &
          / -1.33e0, -0.78e0,  0.31e0, -0.95e0, -0.90e0, -0.30e0/
      data   y1( 31),  y1( 32),  y1( 33),  y1( 34),  y1( 35),  y1( 36) &
          / -1.02e0, -0.53e0,  0.15e0,  1.40e0,  1.22e0,  0.59e0/
      data   y1( 37),  y1( 38),  y1( 39),  y1( 40),  y1( 41),  y1( 42) &
          /  0.70e0,  1.70e0,  2.78e0,  1.98e0,  1.39e0,  1.85e0/
      data   y1( 43),  y1( 44),  y1( 45),  y1( 46),  y1( 47),  y1( 48) &
          /  2.60e0,  0.51e0,  2.77e0,  1.16e0,  1.07e0, -0.48e0/
      data   y1( 49),  y1( 50),  y1( 51),  y1( 52),  y1( 53),  y1( 54) &
          / -0.52e0,  0.37e0,  0.00e0, -1.99e0, -1.75e0,  0.70e0/
      data   y1( 55),  y1( 56),  y1( 57),  y1( 58),  y1( 59),  y1( 60) &
          /  0.73e0,  1.16e0,  0.06e0, -0.02e0,  1.10e0, -0.35e0/
      data   y1( 61),  y1( 62),  y1( 63),  y1( 64),  y1( 65),  y1( 66) &
          / -1.67e0, -1.57e0,  1.16e0,  1.84e0,  3.35e0,  0.40e0/
      data   y1( 67),  y1( 68),  y1( 69),  y1( 70),  y1( 71),  y1( 72) &
          /  0.45e0,  1.30e0,  0.93e0,  1.17e0, -1.74e0, -1.28e0/
      data   y1( 73),  y1( 74),  y1( 75),  y1( 76),  y1( 77),  y1( 78) &
          / -0.07e0,  1.50e0,  0.53e0,  0.20e0, -0.42e0,  1.18e0/
      data   y1( 79),  y1( 80),  y1( 81),  y1( 82),  y1( 83),  y1( 84) &
          /  0.82e0,  1.50e0,  2.92e0,  1.18e0,  1.23e0,  3.16e0/
      data   y1( 85),  y1( 86),  y1( 87),  y1( 88),  y1( 89),  y1( 90) &
          /  0.79e0,  0.68e0,  1.14e0,  1.02e0,  1.02e0, -0.71e0/
      data   y1( 91),  y1( 92),  y1( 93),  y1( 94),  y1( 95),  y1( 96) &
          / -0.17e0, -1.50e0, -0.26e0, -0.38e0,  0.93e0, -0.33e0/
      data   y1( 97),  y1( 98),  y1( 99),  y1(100) &
          / -1.12e0, -2.95e0, -2.09e0, -1.11e0                    /
!
      data  y2a(  1), y2a(  2), y2a(  3), y2a(  4), y2a(  5), y2a(  6) &
          /-0.88e0, -0.16e0, -1.87e0, -1.12e0,  1.38e0,  2.13e0/
      data  y2a(  7), y2a(  8), y2a(  9), y2a( 10), y2a( 11), y2a( 12) &
          / 2.76e0,  0.56e0, -0.69e0, -1.79e0, -3.82e0, -2.38e0/
      data  y2a( 13), y2a( 14), y2a( 15), y2a( 16), y2a( 17), y2a( 18) &
          / 1.00e0,  0.70e0, -0.15e0,  0.98e0,  0.11e0, -0.35e0/
      data  y2a( 19), y2a( 20), y2a( 21), y2a( 22), y2a( 23), y2a( 24) &
          /-0.73e0,  0.89e0, -1.63e0, -0.44e0, -1.37e0, -1.71e0/
      data  y2a( 25), y2a( 26), y2a( 27), y2a( 28), y2a( 29), y2a( 30) &
          /-1.22e0, -2.00e0, -0.22e0,  0.38e0,  1.31e0,  0.71e0/
      data  y2a( 31), y2a( 32), y2a( 33), y2a( 34), y2a( 35), y2a( 36) &
          / 0.32e0,  0.48e0, -1.88e0, -0.94e0, -1.54e0, -0.13e0/
      data  y2a( 37), y2a( 38), y2a( 39), y2a( 40), y2a( 41), y2a( 42) &
          / 1.02e0,  0.02e0, -0.77e0,  0.11e0, -0.60e0, -0.52e0/
      data  y2a( 43), y2a( 44), y2a( 45), y2a( 46), y2a( 47), y2a( 48) &
          /-0.09e0,  1.23e0,  1.46e0,  0.61e0,  0.42e0,  2.16e0/
      data  y2a( 49), y2a( 50), y2a( 51), y2a( 52), y2a( 53), y2a( 54) &
          / 3.18e0,  2.10e0,  0.37e0, -0.24e0,  0.57e0, -0.53e0/
      data  y2a( 55), y2a( 56), y2a( 57), y2a( 58), y2a( 59), y2a( 60) &
          / 2.44e0,  1.02e0, -0.53e0, -2.49e0, -2.12e0, -1.04e0/
      data  y2a( 61), y2a( 62), y2a( 63), y2a( 64), y2a( 65), y2a( 66) &
          /-0.12e0, -1.88e0, -1.50e0,  1.54e0,  3.33e0,  3.08e0/
      data  y2a( 67), y2a( 68), y2a( 69), y2a( 70), y2a( 71), y2a( 72) &
          / 1.71e0,  0.79e0,  1.55e0,  0.89e0, -0.89e0, -1.18e0/
      data  y2a( 73), y2a( 74), y2a( 75), y2a( 76), y2a( 77), y2a( 78) &
          / 0.89e0,  1.71e0,  3.05e0,  0.15e0, -1.04e0,  0.12e0/
      data  y2a( 79), y2a( 80), y2a( 81), y2a( 82), y2a( 83), y2a( 84) &
          / 0.08e0,  0.11e0, -2.62e0, -1.28e0,  1.07e0,  3.20e0/
      data  y2a( 85), y2a( 86), y2a( 87), y2a( 88), y2a( 89), y2a( 90) &
          / 1.92e0,  0.53e0, -1.08e0,  0.49e0, -0.58e0,  0.17e0/
      data  y2a( 91), y2a( 92), y2a( 93), y2a( 94), y2a( 95), y2a( 96) &
          / 1.15e0, -0.97e0, -1.63e0,  1.14e0, -0.67e0, -0.88e0/
      data  y2a( 97), y2a( 98), y2a( 99), y2a(100) &
          /-0.07e0,  0.24e0,  0.55e0, -2.16e0/
!
      data  y2b(  1), y2b(  2), y2b(  3), y2b(  4), y2b(  5), y2b(  6) &
          / 0.79e0,  1.12e0, -1.10e0, -2.39e0, -1.75e0, -0.82e0/
      data  y2b(  7), y2b(  8), y2b(  9), y2b( 10), y2b( 11), y2b( 12) &
          /-0.36e0,  1.27e0,  1.75e0,  2.44e0,  0.36e0, -2.10e0/
      data  y2b( 13), y2b( 14), y2b( 15), y2b( 16), y2b( 17), y2b( 18) &
          /-1.93e0, -1.30e0, -1.75e0, -0.34e0,  0.74e0,  0.49e0/
      data  y2b( 19), y2b( 20), y2b( 21), y2b( 22), y2b( 23), y2b( 24) &
          / 0.70e0,  0.71e0,  0.09e0,  0.59e0,  1.54e0,  0.14e0/
      data  y2b( 25), y2b( 26), y2b( 27), y2b( 28), y2b( 29), y2b( 30) &
          / 0.55e0, -1.40e0, -2.55e0, -1.66e0, -0.43e0,  0.58e0/
      data  y2b( 31), y2b( 32), y2b( 33), y2b( 34), y2b( 35), y2b( 36) &
          / 2.18e0, -0.24e0,  0.58e0, -0.18e0, -1.55e0, -0.64e0/
      data  y2b( 37), y2b( 38), y2b( 39), y2b( 40), y2b( 41), y2b( 42) &
          /-1.09e0,  0.90e0, -0.66e0, -0.35e0,  0.48e0,  0.50e0/
      data  y2b( 43), y2b( 44), y2b( 45), y2b( 46), y2b( 47), y2b( 48) &
          / 0.05e0, -0.68e0,  0.24e0,  0.58e0, -1.26e0, -0.25e0/
      data  y2b( 49), y2b( 50), y2b( 51), y2b( 52), y2b( 53), y2b( 54) &
          / 0.25e0,  2.18e0,  2.96e0,  1.56e0, -0.36e0, -0.59e0/
      data  y2b( 55), y2b( 56), y2b( 57), y2b( 58), y2b( 59), y2b( 60) &
          /-0.12e0,  3.03e0,  2.11e0,  0.78e0,  0.89e0, -1.45e0/
      data  y2b( 61), y2b( 62), y2b( 63), y2b( 64), y2b( 65), y2b( 66) &
          /-0.36e0, -0.37e0, -1.39e0, -4.19e0, -0.73e0, -0.98e0/
      data  y2b( 67), y2b( 68), y2b( 69), y2b( 70), y2b( 71), y2b( 72) &
          / 0.36e0,  0.06e0, -1.94e0, -0.08e0,  0.17e0,  1.00e0/
      data  y2b( 73), y2b( 74), y2b( 75), y2b( 76), y2b( 77), y2b( 78) &
          /-0.05e0,  0.43e0,  0.15e0,  2.69e0,  0.57e0,  0.29e0/
      data  y2b( 79), y2b( 80), y2b( 81), y2b( 82), y2b( 83), y2b( 84) &
          / 1.10e0,  0.48e0, -1.06e0, -2.28e0, -2.03e0, -0.75e0/
      data  y2b( 85), y2b( 86), y2b( 87), y2b( 88), y2b( 89), y2b( 90) &
          / 1.00e0,  1.71e0,  0.58e0,  1.97e0,  0.99e0,  1.94e0/
      data  y2b( 91), y2b( 92), y2b( 93), y2b( 94), y2b( 95), y2b( 96) &
          / 2.18e0,  3.14e0,  0.60e0,  0.51e0,  1.35e0,  0.56e0/
      data  y2b( 97), y2b( 98), y2b( 99), y2b(100) &
          / 0.11e0,  0.00e0,  2.34e0,  1.88e0/
!
      data  y3(  1),y3(  2),y3(  3),y3(  4),y3(  5),y3(  6) &
          /-0.88e0, -0.12e0, -0.89e0, -1.38e0, -0.07e0,  1.03e0/
      data  y3(  7),y3(  8),y3(  9),y3( 10),y3( 11),y3( 12) &
          / 2.14e0,  0.35e0, -1.10e0, -1.78e0, -2.76e0, -1.77e0/
      data  y3( 13),y3( 14),y3( 15),y3( 16),y3( 17),y3( 18) &
          / 0.98e0,  1.00e0, -0.70e0, -1.01e0, -1.30e0, -0.85e0/
      data  y3( 19),y3( 20),y3( 21),y3( 22),y3( 23),y3( 24) &
          /-0.46e0,  1.63e0,  0.06e0, -0.17e0, -1.01e0, -1.04e0/
      data  y3( 25),y3( 26),y3( 27),y3( 28),y3( 29),y3( 30) &
          /-0.66e0, -1.12e0, -0.51e0, -0.71e0, -0.20e0, -0.13e0/
      data  y3( 31),y3( 32),y3( 33),y3( 34),y3( 35),y3( 36) &
          / 0.14e0,  1.59e0, -0.76e0, -1.08e0, -1.77e0, -1.20e0/
      data  y3( 37),y3( 38),y3( 39),y3( 40),y3( 41),y3( 42) &
          / 0.45e0, -0.07e0, -0.63e0, -0.35e0, -0.87e0, -0.62e0/
      data  y3( 43),y3( 44),y3( 45),y3( 46),y3( 47),y3( 48) &
          / 0.28e0,  1.90e0,  2.14e0,  1.05e0,  0.31e0,  1.07e0/
      data  y3( 49),y3( 50) &
          / 2.67e0,  2.44e0/
!
      data  y5(  1), y5(  2), y5(  3), y5(  4), y5(  5), y5(  6) &
          /     5.0e0, 11.0e0, 16.0e0, 23.0e0, 36.0e0, 58.0e0/
      data  y5(  7), y5(  8), y5(  9), y5( 10), y5( 11), y5( 12) &
          /    29.0e0, 20.0e0, 10.0e0,  8.0e0,  3.0e0,  0.0e0/
      data  y5( 13), y5( 14), y5( 15), y5( 16), y5( 17), y5( 18) &
          /     0.0e0, 2.0e0, 11.0e0, 27.0e0, 47.0e0, 63.0e0/
      data  y5( 19), y5( 20), y5( 21), y5( 22), y5( 23), y5( 24) &
          /    60.0e0, 39.0e0, 28.0e0, 26.0e0, 22.0e0, 11.0e0/
      data  y5( 25), y5( 26), y5( 27), y5( 28), y5( 29), y5( 30) &
          /    21.0e0, 40.0e0, 78.0e0,122.0e0,103.0e0, 73.0e0/
      data  y5( 31), y5( 32), y5( 33), y5( 34), y5( 35), y5( 36) &
          /    47.0e0, 35.0e0, 11.0e0,  5.0e0, 16.0e0, 34.0e0/
      data  y5( 37), y5( 38), y5( 39), y5( 40), y5( 41), y5( 42) &
          /    70.0e0, 81.0e0,111.0e0,101.0e0, 73.0e0, 40.0e0/
      data  y5( 43), y5( 44), y5( 45), y5( 46), y5( 47), y5( 48) &
          /    20.0e0, 16.0e0,  5.0e0, 11.0e0, 22.0e0, 40.0e0/
      data  y5( 49), y5( 50), y5( 51), y5( 52), y5( 53), y5( 54) &
          /    60.0e0, 80.9e0, 83.4e0, 47.7e0, 47.8e0, 30.7e0/
      data  y5( 55), y5( 56), y5( 57), y5( 58), y5( 59), y5( 60) &
          /    12.2e0,  9.6e0, 10.2e0, 32.4e0, 47.6e0, 54.0e0/
      data  y5( 61), y5( 62), y5( 63), y5( 64), y5( 65), y5( 66) &
          /    62.9e0, 85.9e0, 61.2e0, 45.1e0, 36.4e0, 20.9e0/
      data  y5( 67), y5( 68), y5( 69), y5( 70), y5( 71), y5( 72) &
          /    11.4e0, 37.8e0, 69.8e0,106.1e0,100.8e0, 81.6e0/
      data  y5( 73), y5( 74), y5( 75), y5( 76), y5( 77), y5( 78) &
          /    66.5e0, 34.8e0, 30.6e0,  7.0e0, 19.8e0, 92.5e0/
      data  y5( 79), y5( 80), y5( 81), y5( 82), y5( 83), y5( 84) &
          /   154.4e0,125.9e0, 84.8e0, 68.1e0, 38.5e0, 22.8e0/
      data  y5( 85), y5( 86), y5( 87), y5( 88), y5( 89), y5( 90) &
          /    10.2e0, 24.1e0, 82.9e0,132.0e0,130.9e0,118.1e0/
      data  y5( 91), y5( 92), y5( 93), y5( 94), y5( 95), y5( 96) &
          /    89.9e0, 66.6e0, 60.0e0, 46.9e0, 41.0e0, 21.3e0/
      data  y5( 97), y5( 98), y5( 99), y5(100), y5(101), y5(102) &
          /    16.0e0,  6.4e0,  4.1e0,  6.8e0, 14.5e0, 34.0e0/
      data  y5(103), y5(104), y5(105), y5(106), y5(107), y5(108) &
          /    45.0e0, 43.1e0, 47.5e0, 42.2e0, 28.1e0, 10.1e0/
      data  y5(109), y5(110), y5(111), y5(112), y5(113), y5(114) &
          /     8.1e0,  2.5e0,  0.0e0,  1.4e0,  5.0e0, 12.2e0/
      data  y5(115), y5(116), y5(117), y5(118), y5(119), y5(120) &
          /    13.9e0, 35.4e0, 45.8e0, 41.1e0, 30.1e0, 23.9e0/
      data  y5(121), y5(122), y5(123), y5(124), y5(125), y5(126) &
          /    15.6e0,  6.6e0,  4.0e0,  1.8e0,  8.5e0, 16.6e0/
      data  y5(127), y5(128), y5(129), y5(130), y5(131), y5(132) &
          /    36.3e0, 49.6e0, 64.2e0, 67.0e0, 70.9e0, 47.8e0/
      data  y5(133), y5(134), y5(135), y5(136), y5(137), y5(138) &
          /    27.5e0,  8.5e0, 13.2e0, 56.9e0,121.5e0,138.3e0/
      data  y5(139), y5(140), y5(141), y5(142), y5(143), y5(144) &
          /   103.2e0, 85.7e0, 64.6e0, 36.7e0, 24.2e0, 10.7e0/
      data  y5(145), y5(146), y5(147), y5(148), y5(149), y5(150) &
          /    15.0e0, 40.1e0, 61.5e0, 98.5e0,124.7e0, 96.3e0/
      data  y5(151), y5(152), y5(153), y5(154), y5(155), y5(156) &
          /    66.6e0, 64.5e0, 54.1e0, 39.0e0, 20.6e0,  6.7e0/
      data  y5(157), y5(158), y5(159), y5(160), y5(161), y5(162) &
          /     4.3e0, 22.7e0, 54.8e0, 93.8e0, 95.8e0, 77.2e0/
      data  y5(163), y5(164), y5(165), y5(166), y5(167), y5(168) &
          /    59.1e0, 44.0e0, 47.0e0, 30.5e0, 16.3e0,  7.3e0/
      data  y5(169), y5(170), y5(171), y5(172), y5(173), y5(174) &
          /    37.6e0, 74.0e0,139.0e0,111.2e0,101.6e0, 66.2e0/
      data  y5(175), y5(176), y5(177), y5(178), y5(179), y5(180) &
          /    44.7e0, 17.0e0, 11.3e0, 12.4e0,  3.4e0,  6.0e0/
      data  y5(181), y5(182), y5(183), y5(184), y5(185), y5(186) &
          /    32.3e0, 54.3e0, 59.7e0, 63.7e0, 63.5e0, 52.2e0/
      data  y5(187), y5(188), y5(189), y5(190), y5(191), y5(192) &
          /    25.4e0, 13.1e0,  6.8e0,  6.3e0,  7.1e0, 35.6e0/
      data  y5(193), y5(194), y5(195), y5(196), y5(197), y5(198) &
          /    73.0e0, 85.1e0, 78.0e0, 64.0e0, 41.8e0, 26.2e0/
      data  y5(199), y5(200), y5(201), y5(202), y5(203), y5(204) &
          /    26.7e0, 12.1e0,  9.5e0,  2.7e0,  5.0e0, 24.4e0/
      data  y5(205), y5(206), y5(207), y5(208), y5(209), y5(210) &
          /    42.0e0, 63.5e0, 53.8e0, 62.0e0, 48.5e0, 43.9e0/
      data  y5(211), y5(212), y5(213), y5(214), y5(215), y5(216) &
          /    18.6e0,  5.7e0,  3.6e0,  1.4e0,  9.6e0, 47.4e0/
      data  y5(217), y5(218), y5(219), y5(220), y5(221), y5(222) &
          /    57.1e0,103.9e0, 80.6e0, 63.6e0, 37.6e0, 26.1e0/
      data  y5(223), y5(224), y5(225), y5(226), y5(227), y5(228) &
          /    14.2e0,  5.8e0, 16.7e0, 44.3e0, 63.9e0, 69.0e0/
      data  y5(229), y5(230), y5(231), y5(232), y5(233), y5(234) &
          /    77.8e0, 64.9e0, 35.7e0, 21.2e0, 11.1e0,  5.7e0/
      data  y5(235), y5(236), y5(237), y5(238), y5(239), y5(240) &
          /     8.7e0, 36.1e0, 79.7e0,114.4e0,109.6e0, 88.8e0/
      data  y5(241), y5(242), y5(243), y5(244), y5(245), y5(246) &
          /    67.8e0, 47.5e0, 30.6e0, 16.3e0,  9.6e0, 33.2e0/
      data  y5(247), y5(248), y5(249), y5(250), y5(251), y5(252) &
          /    92.6e0,151.6e0,136.3e0,134.7e0, 83.9e0, 69.4e0/
      data  y5(253), y5(254), y5(255), y5(256), y5(257), y5(258) &
          /    31.5e0, 13.9e0,  4.4e0, 38.0e0,141.7e0,190.2e0/
      data  y5(259), y5(260), y5(261) &
          /   184.8e0,159.0e0,112.3e0/
!
!  define constants
!
      lper = 300
      lfreq = 300
      lyfft = 600

      ny1 = 100
      ny2 = 50
      ny3 = 50
      ny4 = 50
      ny5 = 261
      ny6 = 100

      nk = 3
      kmd(1) = 8
      kmd(2) = 8
      kmd(3) = 8

      taperp = 0.10

      nfft = 514
      iextnd = 0
      nprt = -1

      ilog = 1
!
!  write header
!
      write ( *, 1000)
!
!  Run simple test of acf
!
      write ( *, 1100)
      call acf ( y1, ny1 )
      write ( *,2000) ierr
!
!  Run simple test of ccf
!
      write ( *, 1200)
      call ccf ( y2a, y2b, ny2 )
      write ( *,2000) ierr
!
!  Run simple test of ufs
!
      write ( *, 1300)
      call ufs ( y3, ny3 )
      write ( *,2000) ierr
!
!  Run simple test of uas
!
      write ( *, 1400)
      call uas ( y4, ny4 )
      write ( *,2000) ierr
!
!  Run simple test of taper
!
      write ( *, 1510)
      call taper ( y5, ny5, taperp, yfft )
      write ( *,2000) ierr
!
!  Run simple test of pgms
!
      write ( *, 1520)
      call pgms ( yfft, ny5, nfft, lyfft, &
                iextnd, nf, per, lper, freq, lfreq, nprt )
      write ( *,2000) ierr
!
!  Run simple test of mdflt
!
      write ( *, 1530)
      call mdflt ( per, nf, nk, kmd, perf, ldstak )
      write ( *,2000) ierr
!
!  display results of mdflt
!
      write ( *, 1540)
      call ppl ( perf, freq, nf, ilog )
      write ( *,2000) ierr
!
!  Run simple test of bfs
!
      write ( *, 1600)
      call bfs ( y6a, y6b, ny6 )
      write ( *,2000) ierr

      return

 1000 format ('1*ch12')
 1100 format (' simple test of acf')
 1200 format ('1simple test of ccf')
 1300 format ('1simple test of ufs')
 1400 format ('1simple test of uas')
 1510 format ('1simple test of taper (no output unless error found)')
 1520 format ('1simple test of pgms')
 1530 format ('1simple test of mdflt (no output unless error found)')
 1540 format ('1display results of mdflt')
 1600 format ('1simple test of bfs')
 2000 format (/' the value of ierr is ', i4)
end
