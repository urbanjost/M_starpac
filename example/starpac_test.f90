program main
use M_starpac !, only : acf,acfd,acff,acffs,acfm,acfms,acfs,iprint
!*****************************************************************************80
!
!! MAIN is the main program for STARPAC_TEST.
!
!  Discussion:
!
!    STARPAC_TEST tests the STARPAC library.
!
!  Modified:
!
!    03 December 2006
!
  implicit none
  integer,parameter :: wp=kind(0.0e0)

  integer ldstak
  parameter ( ldstak = 3000 )

  double precision dstak(ldstak)
  integer ierr
  integer iflag
  integer mbo
  integer mbol
  integer mspect
  integer nfact
  integer nparar
  integer npardf
  integer nparma
  integer nrests
  integer parar
  integer pardf
  integer parma
  real q
  integer t
  integer temp

  save / cstak /
  save / errchk /
  save / mdltsc /
  save / notopt /

  common / cstak / dstak
  common / errchk / ierr
  common / mdltsc / mspect,nfact,pardf,npardf,parar,nparar,parma, &
     nparma,mbo,mbol,t,temp,nrests,iflag
  common / notopt / q

  !!call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STARPAC_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the STARPAC library.'

  call xacf ( ldstak ,typ=0.0_wp)
  call xaimd ( ldstak ,typ=0.0_wp)
  call xaimt ( ldstak ,typ=0.0_wp)
  call xaov1 ( ldstak ,typ=0.0_wp)
  call xbfs ( ldstak ,typ=0.0_wp)
  call xccf ( ldstak ,typ=0.0_wp)
  call xcorr ( ldstak ,typ=0.0_wp)
  call xdckld ( ldstak ,typ=0.0_wp)
  call xdckle ( ldstak ,typ=0.0_wp)
  call xdcklt ( ldstak ,typ=0.0_wp)
  call xdemod ( ldstak ,typ=0.0_wp)
  call xdflt ( ldstak ,typ=0.0_wp)
  call xhist ( ldstak ,typ=0.0_wp)
  call xlls ( ldstak ,typ=0.0_wp)
  call xnlsd ( ldstak ,typ=0.0_wp)
  call xnlse ( ldstak ,typ=0.0_wp)
  call xnlst ( ldstak ,typ=0.0_wp)
  call xnrand ( ldstak ,typ=0.0_wp)
  call xpgm ( ldstak ,typ=0.0_wp)
  call xpp ( ldstak ,typ=0.0_wp )
  call xstat ( ldstak ,typ=0.0_wp)
  call xstpld ( ldstak ,typ=0.0_wp)
  call xstple ( ldstak ,typ=0.0_wp)
  call xstplt ( ldstak ,typ=0.0_wp)
  call xuas ( ldstak ,typ=0.0_wp)
  call xufs ( ldstak ,typ=0.0_wp)
  call xvp ( ldstak, typ=0.0_wp )
  call xxch1 ( ldstak ,typ=0.0_wp)
  call xxch2 ( ldstak, typ=0.0_wp)
  call xxch3 ( ldstak, typ=0.0_wp)
  call xxch4 ( ldstak ,typ=0.0_wp)
  call xxch5 ( ldstak ,typ=0.0_wp)
  call xxch6 ( ldstak ,typ=0.0_wp)
  call xxch7 ( ldstak ,typ=0.0_wp)
  call xxch8 ( ldstak ,typ=0.0_wp)
  call xxch9 ( ldstak ,typ=0.0_wp)
  call xxch10 ( ldstak, typ=0.0_wp )
  call xxch11 ( ldstak ,typ=0.0_wp)
  call xxch12 ( ldstak ,typ=0.0_wp)
  call xxch13 ( ldstak ,typ=0.0_wp)
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STARPAC_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  !!call timestamp ( )

  stop 0
end program main
