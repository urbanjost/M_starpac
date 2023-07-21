subroutine xpp ( typ)
use M_starpac, only : iprint,mpp,mppc,mppl,mppm,mppmc,mppml,pp,ppc,ppl,ppm, &
     ppmc,ppml,scopy,setrv,spp,sppc,sppl,sppm,sppmc,sppml
integer,parameter :: wp=kind(0.0e0)
real(kind=wp),intent(in) :: typ

!*****************************************************************************80
!
!! XPP tests the plotting subroutines.
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
!     integer ierr
!        a common variable used as a flag to indicate whether
!        or not there are any errors, if =0 then no errors.
!     integer ilog
!        the two digit integer, pq, used to select axis scale, where
!        p designates the x-axis and q designates the y-axis.
!        if p==0 (q==0), then the x-axis (y-axis) is linear.
!        if p.ne.0 (q.ne.0), then the x-axis (y-axis) is log.
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
!     integer m
!        the number of vectors in ym
!     integer nout
!        used to indicate how many of the points outside the bounds
!        of the plot are to be listed.
!     integer ny, nym
!        the number of observations in arrays y and ym, respectively.
!     real time(144)
!        the time values for the airline data.
!     real x(144)
!        vector of observations for x(horizontal) coordinates
!     real xlb
!        the lower bound for the x-axis.  (xlb=xub indicates limits are
!        to be determined from the range of the data.)
!     real xmiss
!        the missing value code for the x-axis.
!     real xub
!        the upper bound for the x-axis.  (xlb=xub indicates limits are
!        to be determined from the range of the data.)
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
!  scalars in common
  integer &
     ierr
!
!  local scalars
  real &
     xlb,xmiss,xub,ylb,ymiss,yub
  integer &
     ilog,isize,itest,iym,m,nout,ny,nym
!
!  local arrays
  real &
     air(144),time(144),x(144),y(144),ym(12,12),ymmiss(144)
  integer &
     isym(144)
!
!
!  common blocks
  common /errchk/ierr
!
!  equivalences
  equivalence (y(1),ym(1,1))

  data     xmiss,    ymiss &
      /      7.0e0,    180.0e0/

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
  data time(  1),time(  2),time(  3),time(  4),time(  5),time(  6) &
      /   1.0e0,    2.0e0,    3.0e0,    4.0e0,    5.0e0,    6.0e0/
  data time(  7),time(  8),time(  9),time( 10),time( 11),time( 12) &
      /   7.0e0,    8.0e0,    9.0e0,   10.0e0,   11.0e0,   12.0e0/
  data time( 13),time( 14),time( 15),time( 16),time( 17),time( 18) &
      /  13.0e0,   14.0e0,   15.0e0,   16.0e0,   17.0e0,   18.0e0/
  data time( 19),time( 20),time( 21),time( 22),time( 23),time( 24) &
      /  19.0e0,   20.0e0,   21.0e0,   22.0e0,   23.0e0,   24.0e0/
  data time( 25),time( 26),time( 27),time( 28),time( 29),time( 30) &
      /  25.0e0,   26.0e0,   27.0e0,   28.0e0,   29.0e0,   30.0e0/
  data time( 31),time( 32),time( 33),time( 34),time( 35),time( 36) &
      /  31.0e0,   32.0e0,   33.0e0,   34.0e0,   35.0e0,   36.0e0/
  data time( 37),time( 38),time( 39),time( 40),time( 41),time( 42) &
      /  37.0e0,   38.0e0,   39.0e0,   40.0e0,   41.0e0,   42.0e0/
  data time( 43),time( 44),time( 45),time( 46),time( 47),time( 48) &
      /  43.0e0,   44.0e0,   45.0e0,   46.0e0,   47.0e0,   48.0e0/
  data time( 49),time( 50),time( 51),time( 52),time( 53),time( 54) &
      /  49.0e0,   50.0e0,   51.0e0,   52.0e0,   53.0e0,   54.0e0/
  data time( 55),time( 56),time( 57),time( 58),time( 59),time( 60) &
      /  55.0e0,   56.0e0,   57.0e0,   58.0e0,   59.0e0,   60.0e0/
  data time( 61),time( 62),time( 63),time( 64),time( 65),time( 66) &
      /  61.0e0,   62.0e0,   63.0e0,   64.0e0,   65.0e0,   66.0e0/
  data time( 67),time( 68),time( 69),time( 70),time( 71),time( 72) &
      /  67.0e0,   68.0e0,   69.0e0,   70.0e0,   71.0e0,   72.0e0/
  data time( 73),time( 74),time( 75),time( 76),time( 77),time( 78) &
      /  73.0e0,   74.0e0,   75.0e0,   76.0e0,   77.0e0,   78.0e0/
  data time( 79),time( 80),time( 81),time( 82),time( 83),time( 84) &
      /  79.0e0,   80.0e0,   81.0e0,   82.0e0,   83.0e0,   84.0e0/
  data time( 85),time( 86),time( 87),time( 88),time( 89),time( 90) &
      /  85.0e0,   86.0e0,   87.0e0,   88.0e0,   89.0e0,   90.0e0/
  data time( 91),time( 92),time( 93),time( 94),time( 95),time( 96) &
      /  91.0e0,   92.0e0,   93.0e0,   94.0e0,   95.0e0,   96.0e0/
  data time( 97),time( 98),time( 99),time(100),time(101),time(102) &
      /  97.0e0,   98.0e0,   99.0e0,  100.0e0,  101.0e0,  102.0e0/
  data time(103),time(104),time(105),time(106),time(107),time(108) &
      / 103.0e0,  104.0e0,  105.0e0,  106.0e0,  107.0e0,  108.0e0/
  data time(109),time(110),time(111),time(112),time(113),time(114) &
      / 109.0e0,  110.0e0,  111.0e0,  112.0e0,  113.0e0,  114.0e0/
  data time(115),time(116),time(117),time(118),time(119),time(120) &
      / 115.0e0,  116.0e0,  117.0e0,  118.0e0,  119.0e0,  120.0e0/
  data time(121),time(122),time(123),time(124),time(125),time(126) &
      / 121.0e0,  122.0e0,  123.0e0,  124.0e0,  125.0e0,  126.0e0/
  data time(127),time(128),time(129),time(130),time(131),time(132) &
      / 127.0e0,  128.0e0,  129.0e0,  130.0e0,  131.0e0,  132.0e0/
  data time(133),time(134),time(135),time(136),time(137),time(138) &
      / 133.0e0,  134.0e0,  135.0e0,  136.0e0,  137.0e0,  138.0e0/
  data time(139),time(140),time(141),time(142),time(143),time(144) &
      / 139.0e0,  140.0e0,  141.0e0,  142.0e0,  143.0e0,  144.0e0/
!
  data  air(  1), air(  2), air(  3), air(  4), air(  5), air(  6) &
      / 112.0e0,  118.0e0,  132.0e0,  129.0e0,  121.0e0,  135.0e0/
  data  air(  7), air(  8), air(  9), air( 10), air( 11), air( 12) &
      / 148.0e0,  148.0e0,  136.0e0,  119.0e0,  104.0e0,  118.0e0/
  data  air( 13), air( 14), air( 15), air( 16), air( 17), air( 18) &
      / 115.0e0,  126.0e0,  141.0e0,  135.0e0,  125.0e0,  149.0e0/
  data  air( 19), air( 20), air( 21), air( 22), air( 23), air( 24) &
      / 170.0e0,  170.0e0,  158.0e0,  133.0e0,  114.0e0,  140.0e0/
  data  air( 25), air( 26), air( 27), air( 28), air( 29), air( 30) &
      / 145.0e0,  150.0e0,  178.0e0,  163.0e0,  172.0e0,  178.0e0/
  data  air( 31), air( 32), air( 33), air( 34), air( 35), air( 36) &
      / 199.0e0,  199.0e0,  184.0e0,  162.0e0,  146.0e0,  166.0e0/
  data  air( 37), air( 38), air( 39), air( 40), air( 41), air( 42) &
      / 171.0e0,  180.0e0,  193.0e0,  181.0e0,  183.0e0,  218.0e0/
  data  air( 43), air( 44), air( 45), air( 46), air( 47), air( 48) &
      / 230.0e0,  242.0e0,  209.0e0,  191.0e0,  172.0e0,  194.0e0/
  data  air( 49), air( 50), air( 51), air( 52), air( 53), air( 54) &
      / 196.0e0,  196.0e0,  236.0e0,  235.0e0,  229.0e0,  243.0e0/
  data  air( 55), air( 56), air( 57), air( 58), air( 59), air( 60) &
      / 264.0e0,  272.0e0,  237.0e0,  211.0e0,  180.0e0,  201.0e0/
  data  air( 61), air( 62), air( 63), air( 64), air( 65), air( 66) &
      / 204.0e0,  188.0e0,  235.0e0,  227.0e0,  234.0e0,  264.0e0/
  data  air( 67), air( 68), air( 69), air( 70), air( 71), air( 72) &
      / 302.0e0,  293.0e0,  259.0e0,  229.0e0,  203.0e0,  229.0e0/
  data  air( 73), air( 74), air( 75), air( 76), air( 77), air( 78) &
      / 242.0e0,  233.0e0,  267.0e0,  269.0e0,  270.0e0,  315.0e0/
  data  air( 79), air( 80), air( 81), air( 82), air( 83), air( 84) &
      / 364.0e0,  347.0e0,  312.0e0,  274.0e0,  237.0e0,  278.0e0/
  data  air( 85), air( 86), air( 87), air( 88), air( 89), air( 90) &
      / 284.0e0,  277.0e0,  317.0e0,  313.0e0,  318.0e0,  374.0e0/
  data  air( 91), air( 92), air( 93), air( 94), air( 95), air( 96) &
      / 413.0e0,  405.0e0,  355.0e0,  306.0e0,  271.0e0,  306.0e0/
  data  air( 97), air( 98), air( 99), air(100), air(101), air(102) &
      / 315.0e0,  301.0e0,  356.0e0,  348.0e0,  355.0e0,  422.0e0/
  data  air(103), air(104), air(105), air(106), air(107), air(108) &
      / 465.0e0,  467.0e0,  404.0e0,  347.0e0,  305.0e0,  336.0e0/
  data  air(109), air(110), air(111), air(112), air(113), air(114) &
      / 340.0e0,  318.0e0,  362.0e0,  348.0e0,  363.0e0,  435.0e0/
  data  air(115), air(116), air(117), air(118), air(119), air(120) &
      / 491.0e0,  505.0e0,  404.0e0,  359.0e0,  310.0e0,  337.0e0/
  data  air(121), air(122), air(123), air(124), air(125), air(126) &
      / 360.0e0,  342.0e0,  406.0e0,  396.0e0,  420.0e0,  472.0e0/
  data  air(127), air(128), air(129), air(130), air(131), air(132) &
      / 548.0e0,  559.0e0,  463.0e0,  407.0e0,  362.0e0,  405.0e0/
  data  air(133), air(134), air(135), air(136), air(137), air(138) &
      / 417.0e0,  391.0e0,  419.0e0,  461.0e0,  472.0e0,  535.0e0/
  data  air(139), air(140), air(141), air(142), air(143), air(144) &
      / 622.0e0,  606.0e0,  508.0e0,  461.0e0,  390.0e0,  432.0e0/

  call setrv(ymmiss, 144, ymiss)
  call scopy(144, air, 1, y, 1)
  call scopy(144, time, 1, x, 1)
!
!     define constants
!

  itest = 0
!
!  short calls
!
  ny = 144
  nym = 12
  iym = 12
  m = 12
  ilog = -1
  isize = -1
  nout = -1
  ylb = 0.0e0
  yub = 0.0e0
  xlb = 0.0e0
  xub = 0.0e0

   10 continue
!
!  test of pp
!
  write ( *, 1000)
  write ( *, 3100) itest
  write ( *, 3010) ny
  call pp(y, x, ny)
  write ( *, 3000) ierr
!
!  test of ppm
!
  write ( *, 1030)
  write ( *, 3100) itest
  write ( *, 3010) ny
  call ppm(y, ymiss, x, xmiss, ny)
  write ( *, 3000) ierr
!
!  test of spp
!
  write ( *, 1120)
  write ( *, 3100) itest
  write ( *, 3010) ny
  call spp(y, x, ny, isym)
  write ( *, 3000) ierr
!
!  test of sppm
!
  write ( *, 1150)
  write ( *, 3100) itest
  write ( *, 3010) ny
  call sppm(y, ymiss, x, xmiss, ny, isym)
  write ( *, 3000) ierr
!
!  test of mpp
!
  write ( *, 1060)
  write ( *, 3100) itest
  write ( *, 3010) nym
  write ( *, 3020) m, iym
  call mpp(ym, x, nym, m, iym)
  write ( *, 3000) ierr
!
!  test of mppm
!
  write ( *, 1090)
  write ( *, 3100) itest
  write ( *, 3010) nym
  write ( *, 3020) m, iym
  call mppm(ym, ymmiss, x, xmiss, nym, m, iym)
  write ( *, 3000) ierr
!
!  log option calls
!
   20 continue
!
!  test of ppl
!
  write ( *, 1010)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3040) ilog
  call ppl(y, x, ny, ilog)

  write ( *, 3000) ierr
!
!  test of ppml
!
  write ( *, 1040)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3040) ilog
  call ppml(y, ymiss, x, xmiss, ny, ilog)
  write ( *, 3000) ierr
!
!  test of sppl
!
  write ( *, 1130)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3040) ilog
  call sppl(y, x, ny, isym, ilog)
  write ( *, 3000) ierr
!
!  test of sppml
!
  write ( *, 1160)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3040) ilog
  call sppml(y, ymiss, x, xmiss, ny, isym, ilog)
  write ( *, 3000) ierr
!
!  test of mppl
!
  write ( *, 1070)
  write ( *, 3100) itest
  write ( *, 3010) nym
  write ( *, 3020) m, iym
  write ( *, 3040) ilog
  call mppl(ym, x, nym, m, iym, ilog)
  write ( *, 3000) ierr
!
!  test of mppml
!
  write ( *, 1100)
  write ( *, 3100) itest
  write ( *, 3010) nym
  write ( *, 3020) m, iym
  write ( *, 3040) ilog
  call mppml(ym, ymmiss, x, xmiss, nym, m, iym, ilog)
  write ( *, 3000) ierr
!
!  test of long calls
!
   30 continue
!
!  test of ppc
!
  write ( *, 1020)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3040) ilog
  write ( *, 3050) isize, nout
  write ( *, 3070) ylb, yub, xlb
  write ( *, 3080) xub
  call ppc(y, x, ny, ilog, isize, nout, ylb, yub, xlb, xub)
  write ( *, 3000) ierr
!
!  test of ppmc
!
  write ( *, 1050)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3040) ilog
  write ( *, 3050) isize, nout
  write ( *, 3070) ylb, yub, xlb
  write ( *, 3080) xub
  call ppmc(y, ymiss, x, xmiss, ny, ilog, isize, nout, ylb, yub, &
     xlb, xub)
  write ( *, 3000) ierr
!
!  test of sppc
!
  write ( *, 1140)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3040) ilog
  write ( *, 3050) isize, nout
  write ( *, 3070) ylb, yub, xlb
  write ( *, 3080) xub
  call sppc(y, x, ny, isym, ilog, isize, nout, ylb, yub, xlb, &
     xub)
  write ( *, 3000) ierr
!
!  test of sppmc
!
  write ( *, 1170)
  write ( *, 3100) itest
  write ( *, 3010) ny
  write ( *, 3040) ilog
  write ( *, 3050) isize, nout
  write ( *, 3070) ylb, yub, xlb
  write ( *, 3080) xub
  call sppmc(y, ymiss, x, xmiss, ny, isym, ilog, isize, nout, &
     ylb, yub, xlb, xub)
  write ( *, 3000) ierr
!
!  test of mppc
!
   40 write ( *, 1080)
  write ( *, 3100) itest
  write ( *, 3010) nym
  write ( *, 3020) m, iym
  write ( *, 3040) ilog
  write ( *, 3050) isize, nout
  write ( *, 3070) ylb, yub, xlb
  write ( *, 3080) xub
  call mppc(ym, x, nym, m, iym, ilog, isize, nout, ylb, yub, &
     xlb, xub)
  write ( *, 3000) ierr
!
!  test of mppmc
!
   50 write ( *, 1110)
  write ( *, 3100) itest
  write ( *, 3010) nym
  write ( *, 3020) m, iym
  write ( *, 3040) ilog
  write ( *, 3050) isize, nout
  write ( *, 3070) ylb, yub, xlb
  write ( *, 3080) xub
  call mppmc(ym, ymmiss, x, xmiss, nym, m, iym, ilog, isize, nout, &
     ylb, yub, xlb, xub)
  write ( *, 3000) ierr

  itest = itest + 1

  go to (110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 300), &
    itest
!
!  test valid options
!
  110 ilog = 0
  isize = 0
  nout = 0
  ylb = 100.0e0
  yub = 700.0e0
  xlb = 4.0e0
  xub = 16.0e0
  go to 20

  120 ilog = 2
  isize = 2
  nout = 5
  go to 20

  130 ilog = 20
  isize = 20
  nout = 55
  yub = 300.0e0
  go to 30

  140 ilog = 22
  isize = 22
  go to 40

  150 ny = 1
  nym = 1
  m = 144
  iym = 1
  x(1) = 10.0e0
  go to 40
!
  160 call setrv(y, 144, 1.0e0)
  call setrv(x, 144, 1.0e0)
  nym = 6
  iym = 12
  m = 6
  ny = 36
  ylb = 0.0e0
  yub = 0.0e0
  xlb = 0.0e0
  xub = 0.0e0
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
  x(1) = 0.0e0
  y(1) = 0.0e0
  go to 50

  200 call setrv(x, 144, xmiss)
  call setrv(y, 144, ymiss)
  xlb = xub
  ylb = yub
  go to 50

  300 continue

  return

 1000 format ('1', 10htest of pp)
 1010 format ('1', 11htest of ppl)
 1020 format ('1', 11htest of ppc)
 1030 format ('1', 11htest of ppm)
 1040 format ('1', 12htest of ppml)
 1050 format ('1', 12htest of ppmc)
 1060 format ('1', 11htest of mpp)
 1070 format ('1', 12htest of mppl)
 1080 format ('1', 12htest of mppc)
 1090 format ('1', 12htest of mppm)
 1100 format ('1', 13htest of mppml)
 1110 format ('1', 13htest of mppmc)
 1120 format ('1', 11htest of spp)
 1130 format ('1', 12htest of sppl)
 1140 format ('1', 12htest of sppc)
 1150 format ('1', 12htest of sppm)
 1160 format ('1', 13htest of sppml)
 1170 format ('1', 13htest of sppmc)
 3000 format (/8h ierr = , i4)
 3010 format (' ', 5x, 10h   n     =, i5)
 3020 format ('+', 20x, 10h / m     =, i5, 10h / iym   =, i5)
 3040 format ('+', 65x, 10h / ilog  =, i5)
 3050 format (' ',  5x, 10h   isize =, i5, 10h / nout  =, i5)
 3070 format ('+', 50x, 10h / ylb   =, f10.4, 10h / yub   =, f10.4, &
     10h / xlb   =, f10.4)
 3080 format ('+', 110x, 10h / xub   =, f10.4)
 3100 format (' ', 13h test number , i5)
end
