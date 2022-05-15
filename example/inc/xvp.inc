subroutine xvp ( typ)
use M_starpac, only : iprint,mvp,mvpc,mvpl,mvpm,mvpmc,mvpml,scopy,setrv,svp, &
     svpc,svpl,svpm,svpmc,svpml,vp,vpc,vpl,vpm,vpmc,vpml
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XVP tests the plotting subroutines.
!
!  Discussion:
!
!    series y is the airline data listed on page 531 of box and jenkins.
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
!     real air(144)
!        the airline data.
!     integer ibar
!        the variable used to determine if single points (ibar .ne. 0)
!        or bars (ibar == 0) are to be plotted.
!     integer ierr
!        a common variable used as a flag to indicate whether
!        or not there are any errors, if =0 then no errors.
!     integer ilog
!        the two digit integer, pq, used to select axis scale, where
!        p designates the x-axis and q designates the y-axis.
!        if p==0 (q==0), then the x-axis (y-axis) is linear.
!        if p.ne.0 (q.ne.0), then the x-axis (y-axis) is log.
!     integer irlin
!        the indicator variable used to designate whether zero or the
!        series mean is to be plotted as a reference line, or whether
!        no reference line is to be plotted.
!        if irlin <= -1, no reference line is plotted.
!        if irlin ==  0, zero is plotted as the reference line.
!        if irlin .ge.  1, the series mean is plotted.
!     integer isize
!        the two digit integer, pq, used to select axis size, where
!        p designates the x-axis and q designates the y-axis.
!        if p==0 (q==0), then the x-axis (y-axis) is the maximum.
!        if p.ne.0 (q.ne.0), then the x-axis (y-axis) is half the maximu
!     integer isym(144)
!        vector containing symbol designations for plotting
!     integer itest
!        the number of the test.
!     integer iym
!        actual dimension of ym in users main program
!     integer ldstak
!        *
!     integer m
!        the number of vectors in ym
!     integer ns
!        the sampling frequency,
!        where if ns <= 1, every point is plotted,
!                       = 2, every other point is plotted,
!                       = 3, every third point is plotted, etc.
!     integer ny, nym
!        the number of observations in arrays y and ym, respectively.
!     real xinc
!        the increment for the x axis.
!     real xlb
!        the lower bound for the x-axis.
!     real y(144)
!        vector of observations for the y (vertical) coordinates
!     real ylb
!        the lower bound for the y-axis.  (ylb=yub indicates limits are
!        to be determined from the range of the data.)
!     real ym(12,12)
!        multivariate observations for the y (vertical) coordinates.
!     real ymiss
!        the missing value code for the y-axis.
!     real ymmiss(144)
!        the missing value codes for each column of ym.
!     real yub
!        the upper bound for the y-axis.  (ylb=yub indicates limits are
!        to be determined from the range of the data.)
!
!
!
!  scalars in common
  integer &
     ierr
!
!  local scalars
  real &
     xinc,xlb,ylb,ymiss,yub
  integer &
     ibar,ilog,irlin,isize,itest,iym,m,nout,ns,ny,nym
!
!  local arrays
  real &
     air(144),y(144),ym(12,12),ymmiss(144)
  integer &
     isym(144)
!
!
!  common blocks
  common /errchk/ierr
!
!  equivalences
  equivalence (y(1),ym(1,1))

  data ymiss/180.0e0/
!
  data isym(  1),isym(  2),isym(  3),isym(  4),isym(  5),isym(  6) &
      /    -5000,     6000,        7,        8,        9,       10/
  data isym(  7),isym(  8),isym(  9),isym( 10),isym( 11),isym( 12) &
      /       11,       12,       13,       14,       15,       16/
  data isym( 13),isym( 14),isym( 15),isym( 16),isym( 17),isym( 18) &
      /        5,        6,        7,        8,        9,       10/
  data isym( 19),isym( 20),isym( 21),isym( 22),isym( 23),isym( 24) &
      /       11,       12,       13,       14,       15,       16/
  data isym( 25),isym( 26),isym( 27),isym( 28),isym( 29),isym( 30) &
      /        5,        6,        7,        8,        9,       10/
  data isym( 31),isym( 32),isym( 33),isym( 34),isym( 35),isym( 36) &
      /       11,       12,       13,       14,       15,       16/
  data isym( 37),isym( 38),isym( 39),isym( 40),isym( 41),isym( 42) &
      /        5,        6,        7,        8,        9,       10/
  data isym( 43),isym( 44),isym( 45),isym( 46),isym( 47),isym( 48) &
      /       11,       12,       13,       14,       15,       16/
  data isym( 49),isym( 50),isym( 51),isym( 52),isym( 53),isym( 54) &
      /        5,        6,        7,        8,        9,       10/
  data isym( 55),isym( 56),isym( 57),isym( 58),isym( 59),isym( 60) &
      /       11,       12,       13,       14,       15,       16/
  data isym( 61),isym( 62),isym( 63),isym( 64),isym( 65),isym( 66) &
      /        5,        6,        7,        8,        9,       10/
  data isym( 67),isym( 68),isym( 69),isym( 70),isym( 71),isym( 72) &
      /       11,       12,       13,       14,       15,       16/
  data isym( 73),isym( 74),isym( 75),isym( 76),isym( 77),isym( 78) &
      /        5,        6,        7,        8,        9,       10/
  data isym( 79),isym( 80),isym( 81),isym( 82),isym( 83),isym( 84) &
      /       11,       12,       13,       14,       15,       16/
  data isym( 85),isym( 86),isym( 87),isym( 88),isym( 89),isym( 90) &
      /        5,        6,        7,        8,        9,       10/
  data isym( 91),isym( 92),isym( 93),isym( 94),isym( 95),isym( 96) &
      /       11,       12,       13,       14,       15,       16/
  data isym( 97),isym( 98),isym( 99),isym(100),isym(101),isym(102) &
      /        5,        6,        7,        8,        9,       10/
  data isym(103),isym(104),isym(105),isym(106),isym(107),isym(108) &
      /       11,       12,       13,       14,       15,       16/
  data isym(109),isym(110),isym(111),isym(112),isym(113),isym(114) &
      /        5,        6,        7,        8,        9,       10/
  data isym(115),isym(116),isym(117),isym(118),isym(119),isym(120) &
      /       11,       12,       13,       14,       15,       16/
  data isym(121),isym(122),isym(123),isym(124),isym(125),isym(126) &
      /        5,        6,        7,        8,        9,       10/
  data isym(127),isym(128),isym(129),isym(130),isym(131),isym(132) &
      /       11,       12,       13,       14,       15,       16/
  data isym(133),isym(134),isym(135),isym(136),isym(137),isym(138) &
      /        5,        6,        7,        8,        9,       10/
  data isym(139),isym(140),isym(141),isym(142),isym(143),isym(144) &
      /       11,       12,       13,       14,       15,       16/
!
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
!
  call setrv(ymmiss, 144, ymiss)
  call scopy(144, air, 1, y, 1)

  itest = 0
!
!  short calls
!
  ny = 144
  nym = 12
  iym = 12
  m = 12
  ns = 1
  ilog = -1
  isize = -1
  isize = -1
  irlin = -1
  ibar = -1
  ylb = 0.0e0
  yub = 0.0e0
  xlb = 0.0e0
  xinc = 0.0e0

   10 continue
!
!  test of vp
!
  write ( *, 2000)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3030) ns
  call vp ( y, ny, ns )
  write ( *, 3000) ierr
!
!  test of vpm
!
  write ( *, 2030)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3030) ns
  call vpm ( y, ymiss, ny, ns )
  write ( *, 3000) ierr
!
!  test of svp
!
  write ( *, 2120)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3030) ns
  call svp (y, ny, ns, isym)
  write ( *, 3000) ierr
!
!  test of svpm
!
  write ( *, 2150)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3030) ns
  call svpm (y, ymiss, ny, ns, isym)
  write ( *, 3000) ierr
!
!  test of mvp
!
  write ( *, 2060)
  write ( *, 3100) itest
  write ( *, 3010) nym
  write ( *, 3020) m, iym
  write ( *, 3030) ns
  call mvp (ym, nym, m, iym, ns)
  write ( *, 3000) ierr
!
!  test of mvpm
!
  write ( *, 2090)
  write ( *, 3100) itest
  write ( *, 3010) nym
  write ( *, 3020) m, iym
  write ( *, 3030) ns
  call mvpm (ym, ymmiss, nym, m, iym, ns)
  write ( *, 3000) ierr
!
!     log option calls
!
   20 continue
!
!     test of vpl
!
  write ( *, 2010)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3030) ns
  write ( *, 3040) ilog
  call vpl (y, ny, ns, ilog)
  write ( *, 3000) ierr
!
!     test of vpml
!
  write ( *, 2040)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3030) ns
  write ( *, 3040) ilog
  call vpml (y, ymiss, ny, ns, ilog)
  write ( *, 3000) ierr
!
!     test of svpl
!
  write ( *, 2130)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3030) ns
  write ( *, 3040) ilog
  call svpl (y, ny, ns, isym, ilog)
  write ( *, 3000) ierr
!
!     test of svpml
!
  write ( *, 2160)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3030) ns
  write ( *, 3040) ilog
  call svpml (y, ymiss, ny, ns, isym, ilog)
  write ( *, 3000) ierr
!
!     test of mvpl
!
  write ( *, 2070)
  write ( *, 3100) itest
  write ( *, 3010) nym
  write ( *, 3020) m, iym
  write ( *, 3030) ns
  write ( *, 3040) ilog
  call mvpl (ym, nym, m, iym, ns, ilog)
  write ( *, 3000) ierr
!
!     test of mvpml
!
  write ( *, 2100)
  write ( *, 3100) itest
  write ( *, 3010) nym
  write ( *, 3020) m, iym
  write ( *, 3030) ns
  write ( *, 3040) ilog
  call mvpml(ym, ymmiss, nym, m, iym, ns, ilog)
  write ( *, 3000) ierr
!
!     test of long calls
!
   30 continue
!
!     test of vpc
!
  write ( *, 2020)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3030) ns
  write ( *, 3040) ilog
  write ( *, 3060) isize, irlin, ibar
  write ( *, 3070) ylb, yub, xlb
  write ( *, 3090) xinc
  call vpc (y, ny, ns, ilog, isize, irlin, ibar, ylb, &
     yub, xlb, xinc)
  write ( *, 3000) ierr
!
!     test of vpmc
!
  write ( *, 2050)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3030) ns
  write ( *, 3040) ilog
  write ( *, 3060) isize, irlin, ibar
  write ( *, 3070) ylb, yub, xlb
  write ( *, 3090) xinc
  call vpmc (y, ymiss, ny, ns, ilog, isize, irlin, ibar, ylb, &
     yub, xlb, xinc)
  write ( *, 3000) ierr
!
!     test of svpc
!
  write ( *, 2140)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3030) ns
  write ( *, 3040) ilog
  write ( *, 3060) isize, irlin, ibar
  write ( *, 3070) ylb, yub, xlb
  write ( *, 3090) xinc
  call svpc (y, ny, ns, isym, ilog, isize, irlin, ibar, ylb, &
     yub, xlb, xinc)
  write ( *, 3000) ierr
!
!     test of svpmc
!
  write ( *, 2170)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3030) ns
  write ( *, 3040) ilog
  write ( *, 3060) isize, irlin, ibar
  write ( *, 3070) ylb, yub, xlb
  write ( *, 3090) xinc
  call svpmc(y, ymiss, ny, ns, isym, ilog, isize, irlin, ibar, &
     ylb, yub, xlb, xinc)
  write ( *, 3000) ierr
!
!     test of mvpc
!
   40 write ( *, 2080)
  write ( *, 3100) itest
  write ( *, 3010) nym
  write ( *, 3020) m, iym
  write ( *, 3030) ns
  write ( *, 3040) ilog
  write ( *, 3060) isize, irlin, ibar
  write ( *, 3070) ylb, yub, xlb
  write ( *, 3090) xinc
  call mvpc(ym, nym, m, iym, ns, ilog, isize, ylb, &
     yub, xlb, xinc)
  write ( *, 3000) ierr
!
!     test of mvpmc
!
   50 write ( *, 2110)
  write ( *, 3100) itest
  write ( *, 3010) nym
  write ( *, 3020) m, iym
  write ( *, 3030) ns
  write ( *, 3040) ilog
  write ( *, 3060) isize, irlin, ibar
  write ( *, 3070) ylb, yub, xlb
  write ( *, 3090) xinc
  call mvpmc(ym, ymmiss, nym, m, iym, ns, ilog, isize, ylb, &
     yub, xlb, xinc)
  write ( *, 3000) ierr
!
  itest = itest + 1
!
  go to (110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 300), &
    itest
!
!     test valid options
!
  110 ilog = 0
  isize = 0
  ylb = 100.0e0
  yub = 700.0e0
  xlb = 4.0e0
  xinc = 16.0e0
  go to 20
!
  120 ilog = 2
  isize = 2
  nout = 5
  xinc = -1.0e0
  go to 20
!
  130 ilog = 20
  isize = 20
  nout = 55
  yub = 300.0e0
  go to 30
!
  140 ilog = 22
  isize = 22
  go to 40
!
  150 ny = 1
  nym = 1
  m = 144
  iym = 1
  go to 40
!
  160 call setrv(y, 144, 1.0e0)
  nym = 6
  iym = 12
  m = 6
  ny = 36
  ylb = 0.0e0
  yub = 0.0e0
  xlb = 0.0e0
  xinc = 0.0e0
  go to 30
!
!     test error response
!
  170 ny = 0
  nym = 0
  m = 0
  iym = -1
  go to 10

  180 ny = 144
  nym = 12
  m = 12
  iym = -1
  xlb = -1.0e0
  ylb = -1.0e0
  go to 40

  190 iym = 12
  y(1) = 0.0e0
  go to 50

  200 call setrv(y, 144, ymiss)
  xlb = xinc
  ylb = yub
  go to 50

  300 continue

  return

 2000 format ('test of vp')
 2010 format ('test of vpl')
 2020 format ('test of vpc')
 2030 format ('test of vpm')
 2040 format ('test of vpml')
 2050 format ('test of vpmc')
 2060 format ('test of mvp')
 2070 format ('test of mvpl')
 2080 format ('test of mvpc')
 2090 format ('test of mvpm')
 2100 format ('test of mvpml')
 2110 format ('test of mvpmc')
 2120 format ('test of svp')
 2130 format ('test of svpl')
 2140 format ('test of svpc')
 2150 format ('test of svpm')
 2160 format ('test of svpml')
 2170 format ('test of svpmc')
 3000 format (/' ierr = ', i4)
 3010 format (' ', 5x, '   n     =', i5)
 3020 format ('+', 20x, 10h / m     =, i5, 10h / iym   =, i5)
 3030 format ('+', 50x, 10h / ns    =, i5)
 3040 format ('+', 65x, 10h / ilog  =, i5)
 3060 format (' ',  5x, '   isize=', i5, ' / irlin=', i5, &
     10h / ibar  =, i5)
 3070 format ('+', 50x, 10h / ylb   =, f10.4, 10h / yub   =, f10.4, &
     10h / xlb   =, f10.4)
 3090 format ('+', 110x, 10h / xinc  =, f10.4)
 3100 format (' test number ', i5)
end
