*STPAMO
      SUBROUTINE STPAMO(HEAD, N, EXM, NEXMPT, NETA, J, PAR, NPAR, STP,
     +   NFAIL, IFAIL, SCALE, LSCALE, HDR, PAGE, WIDE, ISUBHD, NPRT,
     +   PRTFXD, IFIXD)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS A DUMMY ROUTINE FOR THE ARIMA ESTIMATION ROUTINES
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  APRIL 2, 1981
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      REAL
     +   EXM
      INTEGER
     +   ISUBHD,J,LSCALE,N,NETA,NEXMPT,NPAR,NPRT
      LOGICAL
     +   HEAD,PAGE,PRTFXD,WIDE
C
C  ARRAY ARGUMENTS
      REAL
     +   PAR(NPAR),SCALE(LSCALE),STP(NPAR)
      INTEGER
     +   IFAIL(N),IFIXD(NPAR),NFAIL(NPAR)
C
C  SUBROUTINE ARGUMENTS
      EXTERNAL HDR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     REAL EXM
C        THE PROPORTION OF OBSERVATIONS ACTUALLY USED FOR WHICH THE
C        COMPUTED NUMERICAL DERIVATIVES WRT A GIVEN PARAMETER ARE
C        EXEMPTED FROM MEETING THE DERIVATIVE ACCEPTANCE CRITERIA.
C     EXTERNAL HDR
C        THE NAME OF THE ROUTINE WHICH PRODUCES THE HEADING
C     LOGICAL HEAD
C        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
C        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
C        OF HEAD WILL BE CHANGED TO FALSE.
C     INTEGER IFAIL(N)
C        THE ARRAY OF INDICATOR VARIABLES DESIGNATING WHETHER
C        THE STEP SIZE SELECTED WAS SATISFACTORY FOR A GIVEN
C        OBSERVATION AND PARAMETER.
C     INTEGER IFIXD(NPAR)
C        THE INDICATOR VALUES USED TO DESIGNATE WHETHER THE
C        PARAMETERS ARE TO BE OPTIMIZED OR ARE TO BE HELD FIXED.  IF
C        IFIXD(I).NE.0, THEN PAR(I) WILL BE OPTIMIZED.  IF
C        IFIXD(I).EQ.0, THEN PAR(I) WILL BE HELD FIXED.
C     INTEGER ISUBHD
C     INTEGER J
C        THE INDEX OF THE PARAMETER BEING EXAMINED.
C     INTEGER LSCALE
C        THE LENGTH OF VECTOR SCALE.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS.
C     INTEGER NPAR
C        THE NUMBER OF PARAMETERS IN THE MODEL.
C     INTEGER NETA
C        THE NUMBER OF RELIABLE DIGITS IN THE MODEL.
C     INTEGER NEXMPT
C        THE NUMBER OF OBSERVATIONS FOR WHICH A GIVEN STEP SIZE
C        DOES NOT HAVE TO BE SATISFACTORY AND THE SELECTED STEP
C        SIZE STILL BE CONSIDERED OK.
C     INTEGER NFAIL(NPAR)
C        THE NUMBER OF OBSERVATIONS FOR WHICH THE SELECTED STEP
C        SIZE DOES NOT MEET THE CRITERIA.
C     INTEGER NPRT
C        THE INDICATOR VARIABLE USED TO SPECIFY WHETHER OR NOT
C        PRINTED OUTPUT IS TO BE PROVIDED, WHERE IF THE VALUE OF
C        NPRT IS ZERO, NO PRINTED OUTPUT IS GIVEN.
C     LOGICAL PAGE
C        THE VARIABLE USED TO INDICATE WHETHER OR NOT THE OUTPUT
C        IS TO BEGIN ON A NEW PAGE.
C     REAL PAR(NPAR)
C        THE ARRAY IN WHICH THE CURRENT ESTIMATES OF THE
C        PARAMETERS ARE STORED.
C     LOGICAL PRTFXD
C        THE INDICATOR VALUE USED TO DESIGNATE WHETHER THE
C        OUTPUT IS TO INCLUDE INFORMATION ON WHETHER THE
C        PARAMETER IS FIXED (TRUE) OR NOT (FALSE).
C     REAL SCALE(LSCALE)
C        THE TYPICAL SIZE OF THE PARAMETERS.
C     REAL STP(NPAR)
C        THE SELECTED STEP SIZE.
C     LOGICAL WIDE
C        THE VARIABLE USED TO INDICATE WHETHER THE HEADING SHOULD
C        BE FULL WIDTH (TRUE) OR NOT (FALSE).
C
      RETURN
C
      END
