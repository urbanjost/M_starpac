*DCKDRV
      SUBROUTINE DCKDRV (NMSUB, LDSTAK, XM, N, M, IXM, MDL,
     +   DRV, PAR, NPAR, NETA, NTAU, SCALE, LSCALE, NROW, NPRT)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE IS THE DRIVER FOR THE DERIVATIVE CHECKING ROUTINES.
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
      INTEGER
     +   IXM,LDSTAK,LSCALE,M,N,NETA,NPAR,NPRT,NROW,NTAU
C
C  ARRAY ARGUMENTS
      REAL
     +   PAR(*),SCALE(*),XM(*)
      CHARACTER
     +   NMSUB(6)*1
C
C  SUBROUTINE ARGUMENTS
      EXTERNAL DRV,MDL
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
     +   ISUBHD,LIFIXD
      LOGICAL
     +   HLFRPT,PAGE,PRTFXD,WIDE
C
C  LOCAL ARRAYS
      INTEGER
     +   IFIXED(1)
C
C  EXTERNAL SUBROUTINES
      EXTERNAL DCKCNT,DCKER,DCKHDR,STKSET
C
C  COMMON BLOCKS
      COMMON /CSTAK/DSTAK
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     EXTERNAL DCKHDR
C        THE NAME OF THE ROUTINE WHICH PRODUCES THE HEADING
C     EXTERNAL DRV
C        THE NAME OF THE USER SUPPLIED SUBROUTINE WHICH COMPUTES THE
C        ANALYTIC DERIVATIVES (JACOBIAN MATRIX) OF THE MODEL.
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     LOGICAL HLFRPT
C        THE VARIABLE WHICH INDICATES WHETHER THE DERIVATIVE
C        CHECKING ROUTINE HAS ALREADY PRINTED PART OF THE
C        INITIAL SUMMARY (TRUE) OR NOT (FALSE).
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED
C        IF IERR .EQ. 1, ERRORS HAVE BEEN DETECTED
C     INTEGER IFIXED(1)
C        THE INDICATOR VALUES USED TO DESIGNATE WHETHER THE
C        PARAMETERS ARE TO BE OPTIMIZED OR ARE TO BE HELD FIXED.  IF
C        IFIXED(I).NE.0, THEN PAR(I) WILL BE OPTIMIZED.  IF
C        IFIXED(I).EQ.0, THEN PAR(I) WILL BE HELD FIXED.
C     INTEGER ISUBHD
C        AN INDICATOR VALUE SPECIFYING SUBHEADINGS TO BE PRINTED.
C     INTEGER IXM
C        THE FIRST DIMENSION OF THE INDEPENDENT VARIABLE ARRAY.
C     INTEGER LDSTAK
C        THE LENGTH OF THE ARRAY DSTAK.
C     INTEGER LIFIXD
C        THE LENGTH OF THE VECTOR IFIXED.
C     INTEGER LSCALE
C        THE LENGTH OF VECTOR SCALE.
C     INTEGER M
C        THE NUMBER OF INDEPENDENT VARIABLES.
C     EXTERNAL MDL
C        THE NAME OF THE USER SUPPLIED SUBROUTINE WHICH COMPUTES THE
C        PREDICTED VALUES BASED ON THE CURRENT PARAMETER ESTIMATES.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS OF DATA.
C     INTEGER NETA
C        THE NUMBER OF ACCURATE DIGITS IN THE MODEL.
C     CHARACTER*1 NMSUB(6)
C        THE NAME OF THE SUBROUTINE CALLING THE ERROR CHECKING
C        SUBROUTINES.
C     INTEGER NPAR
C        THE NUMBER OF UNKNOWN PARAMETERS IN THE MODEL.
C     INTEGER NPRT
C        THE INDICATOR VARIABLE USED TO SPECIFY WHETHER OR NOT
C        PRINTED OUTPUT IS TO BE PROVIDED, WHERE IF THE VALUE OF
C        NPRT IS ZERO, NO PRINTED OUTPUT IS GIVEN.
C     INTEGER NROW
C        THE USER-SUPPLIED NUMBER OF THE ROW OF THE INDEPENDENT
C        VARIABLE ARRAY AT WHICH THE DERIVATIVE IS TO BE CHECKED.
C     INTEGER NTAU
C        THE NUMBER OF DIGITS OF AGREEMENT REQUIRED BETWEEN THE
C        NUMERICAL DERIVATIVES AND THE USER SUPPLIED DERIVATIVES.
C     LOGICAL PAGE
C        THE VARIABLE USED TO INDICATE WHETHER A GIVEN SECTION OF
C        THE OUTPUT IS TO BEGIN ON A NEW PAGE (TRUE) OR NOT (FALSE).
C     REAL PAR(NPAR)
C        THE ARRAY IN WHICH THE CURRENT ESTIMATES OF THE UNKNOWN
C        PARAMETERS ARE STORED.
C     LOGICAL PRTFXD
C        THE INDICATOR VALUE USED TO DESIGNATE WHETHER THE
C        OUTPUT IS TO INCLUDE INFORMATION ON WHETHER THE
C        PARAMETER IS FIXED (TRUE) OR NOT (FALSE).
C     REAL SCALE(LSCALE)
C        THE TYPICAL SIZE OF THE UNKNOWN PARAMETERS.
C     LOGICAL WIDE
C        THE VARIABLE USED TO INDICATE WHETHER THE HEADING SHOULD
C        BE FULL WIDTH (TRUE) OR NOT (FALSE).
C     REAL XM(IXM,M)
C        THE INDEPENDENT VARIABLE ARRAY
C
C
C     CHECK FOR ERRORS IN INPUT PARAMETERS
C
      CALL DCKER(NMSUB, N, M, IXM, NPAR, LDSTAK, SCALE, LSCALE)
C
      IF (IERR.NE.0) RETURN
C
      PAGE = .FALSE.
      WIDE = .TRUE.
      ISUBHD = 0
C
      PRTFXD = .FALSE.
      IFIXED(1) = -1
      LIFIXD = 1
C
      CALL STKSET(LDSTAK, 4)
C
C     PASS CONTROL OF DERIVATIVE CHECKING TO DCKCNT
C
      CALL DCKCNT (XM, N, M, IXM, MDL, DRV, PAR, NPAR, NETA,
     +   NTAU, SCALE, LSCALE, NROW, NPRT, DCKHDR, PAGE, WIDE, ISUBHD,
     +   HLFRPT, PRTFXD, IFIXED, LIFIXD)
C
      RETURN
C
      END
