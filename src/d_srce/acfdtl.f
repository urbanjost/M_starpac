*ACFDTL
      SUBROUTINE ACFDTL (NDF, ND, IOD, NTIMES)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE PRINTS TITLING FOR ACORRD.
C
C     WRITTEN BY - JANET R. DONALDSON
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  NOVEMBER 21, 1980
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   NDF,NTIMES
C
C  ARRAY ARGUMENTS
      INTEGER
     +   IOD(*),ND(*)
C
C  LOCAL SCALARS
      INTEGER
     +   I,IPRT,ISTOP
      CHARACTER
     +   ICOM*1,IPER*1,IPUNCT*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL IPRINT
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     INTEGER I
C        AN INDEX VARIABLE.
C     CHARACTER*1 ICOM
C        THE HOLLERITH VALUE -,- (COMMA)
C     INTEGER IOD(NDF)
C        THE ORDER OF EACH OF THE DIFFERENCE FACTORS.
C     CHARACTER*1 IPER
C        THE HOLLERITH VALUE -.- (PERIOD)
C     INTEGER IPRT
C        THE UNIT NUMBER OF THE DEVICE USED FOR PRINTED
C        OUTPUT.
C     CHARACTER*1 IPUNCT
C        THE HOLLERITH VALUE OF EITHER COMMA OR PERIOD.
C     INTEGER ISTOP
C        ONE LESS THAN THE NUMBER OF DIFFERENCE FACTORS.
C     INTEGER ND(NDF)
C        THE ARRAY CONTAINING THE NUMBER OF TIMES THE DIFFERENCE
C        FACTORS ARE TO BE APPLIED.
C     INTEGER NDF
C        THE NUMBER OF DIFFERENCE FACTORS.
C     INTEGER NTIMES
C        THE NUMBER OF TIMES THE DIFFERENCING FACTOR HAS BEEN APPLIED.
C
      DATA ICOM/','/, IPER/'.'/
C
      CALL IPRINT (IPRT)
C
      IF (NDF .LE. 1) GO TO 10
C
      ISTOP = NDF - 1
      IPUNCT = IPER
      IF (NTIMES .GE. 1) IPUNCT = ICOM
      WRITE(IPRT, 1000)
      IF (NDF .EQ. 2)  WRITE(IPRT, 1001) ND(2), IOD(2), IPER
      IF (NDF .GE. 3) WRITE(IPRT, 1001)
     +   (ND(I), IOD(I), ICOM, I = 1, ISTOP), ND(NDF), IOD(NDF), IPUNCT
      GO TO 20
C
   10 WRITE(IPRT, 1002)
C
   20 IF (NTIMES .EQ. 0) RETURN
C
      IF (NDF .GE. 2) WRITE(IPRT, 1003) NTIMES, IOD(1)
      IF (NDF .EQ. 1) WRITE(IPRT, 1004) NTIMES, IOD(1)
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT(//47H SERIES ANALYZED IS INPUT SERIES DIFFERENCED BY/)
 1001 FORMAT(3X, 3(I3, ' FACTOR(S) OF ORDER ', I3, A1, 1X)/)
 1002 FORMAT(//' SERIES ANALYZED IS ORIGINAL INPUT SERIES'/)
 1003 FORMAT(4X, 34H AND, IN ADDITION, DIFFERENCED BY , I3,
     +   18H FACTORS OF ORDER , I3, '.'//)
 1004 FORMAT(4X, 16H DIFFERENCED BY , I3, 18H FACTORS OF ORDER ,
     +   I3, '.'//)
      END
