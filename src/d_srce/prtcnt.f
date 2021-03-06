*PRTCNT
      SUBROUTINE PRTCNT(NPRT, NDIGIT, IPTOUT)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE SETS UP THE PRINT CONTROL PARAMETERS.
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  DECEMBER 29, 1982
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   NDIGIT,NPRT
C
C  ARRAY ARGUMENTS
      INTEGER
     +   IPTOUT(NDIGIT)
C
C  LOCAL SCALARS
      INTEGER
     +   I,IFAC1,IFAC2
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MOD
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     INTEGER I, IFAC1, IFAC2
C     INTEGER IPTOUT(NDIGIT)
C        THE VARIABLE USED TO CONTROL PRINTED OUTPUT FOR EACH SECTION.
C     INTEGER NDIGIT
C        THE NUMBER OF DIGITS IN THE PRINT CONTROL VALUE.
C     INTEGER NPRT
C        THE PARAMETER USED TO INDICATE HOW MUCH PRINTED OUTPUT IS
C        TO BE PROVIDED.
C
C
      IF (NPRT.LE.-1) GO TO 20
C
      IFAC1 = 10 ** (NDIGIT)
      DO 10 I = 1, NDIGIT
         IFAC2 = IFAC1/10
         IPTOUT(I) = MOD(NPRT, IFAC1) / IFAC2
         IFAC1 = IFAC2
   10 CONTINUE
      RETURN
C
   20 DO 30 I = 1, NDIGIT
         IPTOUT(I) = 1
   30 CONTINUE
      IPTOUT (NDIGIT) = 2
C
      RETURN
C
      END
