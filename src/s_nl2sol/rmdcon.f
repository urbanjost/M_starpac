*RMDCON
      REAL FUNCTION RMDCON(K)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   K
C
C  LOCAL SCALARS
      REAL
     +   BIG,ETA,MACHEP,ONE001,PT999
C
C  EXTERNAL FUNCTIONS
      REAL
     +   R1MACH
      EXTERNAL R1MACH
C
C  INTRINSIC FUNCTIONS
      INTRINSIC SQRT
C
C
C  ***  RETURN MACHINE DEPENDENT CONSTANTS USED BY NL2SOL  ***
C
C +++  COMMENTS BELOW CONTAIN DATA STATEMENTS FOR VARIOUS MACHINES.  +++
C +++  TO CONVERT TO ANOTHER MACHINE, PLACE A C IN COLUMN 1 OF THE   +++
C +++  DATA STATEMENT LINE(S) THAT CORRESPOND TO THE CURRENT MACHINE +++
C +++  AND REMOVE THE C FROM COLUMN 1 OF THE DATA STATEMENT LINE(S)  +++
C +++  THAT CORRESPOND TO THE NEW MACHINE.                           +++
C
C     INTEGER K
C
C  ***  THE CONSTANT RETURNED DEPENDS ON K...
C
C  ***        K = 1... SMALLEST POS. ETA SUCH THAT -ETA EXISTS.
C  ***        K = 2... SQUARE ROOT OF 1.001*ETA.
C  ***        K = 3... UNIT ROUNDOFF = SMALLEST POS. NO. MACHEP SUCH
C  ***                 THAT 1 + MACHEP .GT. 1 .AND. 1 - MACHEP .LT. 1.
C  ***        K = 4... SQUARE ROOT OF 0.999*MACHEP.
C  ***        K = 5... SQUARE ROOT OF 0.999*BIG (SEE K = 6).
C  ***        K = 6... LARGEST MACHINE NO. BIG SUCH THAT -BIG EXISTS.
C
C
      DATA ONE001/1.001E0/, PT999/0.999E0/
C
      BIG = R1MACH(2)
      ETA = R1MACH(1)
      MACHEP = R1MACH(4)
C
C-------------------------------  BODY  --------------------------------
C
      GO TO (10, 20, 30, 40, 50, 60), K
C
 10   RMDCON = ETA
      GO TO 999
C
 20   RMDCON = SQRT(ONE001*ETA)
      GO TO 999
C
 30   RMDCON = MACHEP
      GO TO 999
C
 40   RMDCON = SQRT(PT999*MACHEP)
      GO TO 999
C
 50   RMDCON = SQRT(PT999*BIG)
      GO TO 999
C
 60   RMDCON = BIG
C
 999  RETURN
C  ***  LAST CARD OF RMDCON FOLLOWS  ***
      END
