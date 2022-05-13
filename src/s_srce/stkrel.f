*STKREL
      SUBROUTINE STKREL(NUMBER)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C  DE-ALLOCATES THE LAST (NUMBER) ALLOCATIONS MADE IN THE STACK
C  BY STKGET.
C
C  ERROR STATES -
C
C    1 - NUMBER .LT. 0
C    2 - LNOW, LUSED, LMAX OR LBOOK OVERWRITTEN
C    3 - ATTEMPT TO DE-ALLOCATE NON-EXISTENT ALLOCATION
C    4 - THE POINTER AT ISTAK(LNOW) OVERWRITTEN
C
C     THIS FUNCTION WAS ADAPTED FROM THE FRAMEWORK FUNCTION ISTKGT
C
C     ADAPTED BY - JANET R. DONALDSON
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  NOVEMBER 26, 1980
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   NUMBER
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  ARRAYS IN COMMON
      DOUBLE PRECISION DSTAK(12)
C
C  LOCAL SCALARS
      INTEGER
     +   IN,IPRT,LBOOK,LMAX,LNOW,LOUT,LUSED
C
C  LOCAL ARRAYS
      INTEGER
     +   ISTAK(12)
C
C  EXTERNAL SUBROUTINES
      EXTERNAL IPRINT
C
C  COMMON BLOCKS
      COMMON /CSTAK/DSTAK
      COMMON /ERRCHK/IERR
C
C  EQUIVALENCES
      EQUIVALENCE (DSTAK(1),ISTAK(1))
      EQUIVALENCE (ISTAK(1),LOUT)
      EQUIVALENCE (ISTAK(2),LNOW)
      EQUIVALENCE (ISTAK(3),LUSED)
      EQUIVALENCE (ISTAK(4),LMAX)
      EQUIVALENCE (ISTAK(5),LBOOK)
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED
C        IF IERR .EQ. 1, ERRORS HAVE BEEN DETECTED
C     INTEGER IN
C        ...
C     INTEGER IPRT
C        THE UNIT NUMBER FOR PRINTED OUTPUT.
C     INTEGER ISTAK(12)
C        THE INTEGER VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER LBOOK
C        THE NUMBER OF WORDS USED FOR BOOKEEPING.
C     INTEGER LMAX
C        THE MAXIMUM LENGTH OF THE STACK.
C     INTEGER LNOW
C        THE CURRENT ACTIVE LENGTH OF THE STACK.
C     INTEGER LOUT
C        THE NUMBER OF CURRENT ALLOCATIONS.
C     INTEGER LUSED
C        THE MAXIMUM VALUE OF ISTAK(2) ACHEIVED.
C     INTEGER NUMBER
C        THE NUMBER OF ALLOCATIONS TO BE FREED FROM THE STACK.
C
C
      IF (LNOW.LT.LBOOK.OR.LNOW.GT.LUSED.OR.LUSED.GT.LMAX) GO TO 20
C
      IN = NUMBER
 10      IF (IN.EQ.0) RETURN
C
         IF (LNOW.LE.LBOOK) GO TO 30
C
C     CHECK TO MAKE SURE THE BACK POINTERS ARE MONOTONE.
C
         IF (ISTAK(LNOW).LT.LBOOK.OR.ISTAK(LNOW).GE.LNOW-1) GO TO 40
C
         LOUT = LOUT-1
         LNOW = ISTAK(LNOW)
         IN = IN-1
         GO TO 10
C
C     PRINT ERROR MESSAGES
C
   20 IERR = 1
      CALL IPRINT(IPRT)
      WRITE (IPRT, 1000)
      RETURN
C
   30 IERR = 1
      CALL IPRINT(IPRT)
      WRITE (IPRT, 1010)
      RETURN
C
   40 IERR = 1
      CALL IPRINT(IPRT)
      WRITE (IPRT, 1020) LOUT
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (///18H ***** ERROR *****//
     +   50H DSTAK BOOKKEEPING ELEMENTS HAVE BEEN OVERWRITTEN.)
 1010 FORMAT (///18H ***** ERROR *****//
     +   52H ATTEMPT HAS BEEN MADE TO DE-ALLOCATE A NON-EXISTANT,
     +   21H ALLOCATION IN DSTAK.)
 1020 FORMAT (///18H ***** ERROR *****//
     +   35H THE POINTER FOR ALLOCATION NUMBER , I3, 9H HAS BEEN,
     +   13H OVERWRITTEN.)
C
      END
