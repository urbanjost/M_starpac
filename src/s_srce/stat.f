*STAT
      SUBROUTINE STAT(Y, N, LDSTAK)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS SUBROUTINE COMPUTES 53 DIFFERENT STATISTICS FOR A VECTOR
C     Y, WITH NO WEIGHTS SPECIFIED.  ONE PAGE OF AUTOMATIC
C     PRINTOUT IS PRODUCED.
C
C     WRITTEN BY - JANET R. DONALDSON, JOHN E. KOONTZ
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  MAY 17, 1982
C        (EXTENSIVE REVISION OF OLDER VERSION)
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   LDSTAK,N
C
C  ARRAY ARGUMENTS
      REAL
     +   Y(*)
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  ARRAYS IN COMMON
      DOUBLE PRECISION DSTAK(12)
C
C  LOCAL SCALARS
      REAL
     +   ALPHA,SUMD2,SUMD3,SUMD4,SUMDA,SUMDI,SUMT1
      INTEGER
     +   IDP,IINT,IPRT,LSORT,MID,NALL0,NNZW
      LOGICAL
     +   STACK,WTS
C
C  LOCAL ARRAYS
      REAL
     +   STS(53),WT(1)
      INTEGER
     +   ISTAK(12)
      CHARACTER
     +   NMSUB(6)*1
C
C  EXTERNAL FUNCTIONS
      INTEGER
     +   STKGET,STKST
      EXTERNAL STKGET,STKST
C
C  EXTERNAL SUBROUTINES
      EXTERNAL GENI,IPRINT,SRTIR,SRTRI,STAT1,STAT2,STATER,STKCLR,STKSET,
     +   SUMBS,SUMDS,SUMID,SUMOT,SUMSS,SUMTS
C
C  COMMON BLOCKS
      COMMON /CSTAK/DSTAK
      COMMON /ERRCHK/IERR
C
C  EQUIVALENCES
      EQUIVALENCE (DSTAK(1),ISTAK(1))
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     REAL ALPHA
C        THE PERCENTAGE TO BE TRIMMED OFF EACH END OF Y FOR THE
C        TRIMMED MEANS CALCULATIONS.
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER IDP
C        THE CODE VALUE FOR DOUBLE PRECISION FOR FRAMEWORK.
C     INTEGER IERR
C        THE CODE INDICATING WHETHER OR NOT AN ERROR HAS
C        BEEN DISCOVERED.  0 MEANS NO ERROR, NOT 0 MEANS
C        SOME ERROR EXISTS.
C     INTEGER IINT
C        THE CODE VALUE FOR INTEGER FOR FRAMEWORK
C     INTEGER IPRT
C        THE NUMBER OF THE STANDARD OUTPUT UNIT.
C     INTEGER ISTAK(12)
C        THE INTEGER VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER LDSTAK
C        INPUT PARAMETER.  THE NUMBER OF DOUBLE PRECISION
C        ELEMENTS DIMENSIONED FOR DSTAK IN THE USER PROGRAM.
C     INTEGER LSORT
C        THE STARTING LOCATION IN ISTAK OF THE PERMUTATION
C        VECTOR.
C     INTEGER MID
C        THE INDEX OF THE (AN) ELEMENT OF Y CLOSEST TO ZERO, WHEN
C        Y HAS BEEN SORTED.
C     INTEGER N
C        INPUT PARAMETER.  THE LENGTH OF Y.
C     INTEGER NALL0
C        THE NUMBER OF ALLOCATIONS OUTSTANDING AT THE TIME THIS ROUTINE
C        WAS CALLED.
C     CHARACTER*1 NMSUB(6)
C        THE ARRAY CONTAINING THE NAME OF THIS SUBROUTINE.
C     INTEGER NNZW
C        NUMBER OF NONZERO WEIGHTS.
C     LOGICAL STACK
C        A FLAG INDICATING WHETHER THIS ROUTINE USES THE STACK (TRUE)
C        OR NOT (FALSE).
C     REAL STS(53)
C        THE VECTOR OF THE 53 STATISTICS COMPUTED.
C     REAL SUMDA
C        THE SUM OF THE ABSOLUTE DIFFERENCES FROM THE MEAN.
C     REAL SUMDI
C        THE SUM OF THE PRODUCTS OF THE INDICES AND THE DIFFERENCES.
C     REAL SUMD2
C        THE SUM OF THE SQUARES OF THE DIFFERENCES.
C     REAL SUMD3
C        THE SUM OF THE CUBES OF THE DIFFERENCES.
C     REAL SUMD4
C        THE SUM OF THE 4TH POWERS OF THE DIFFERENCES.
C     REAL SUMT1
C        THE SUM OF THE ALPHA TRIMMED ARRAY Y.
C     REAL WT(1)
C        THE DUMMY WEIGHTS VECTOR.
C     LOGICAL WTS
C        A FLAG INDICATING WHETHER THERE ARE WEIGHTS (TRUE)
C        OR NOT (FALSE).
C     REAL Y(N)
C        INPUT PARAMETER.  THE VECTOR OF DATA POINTS ON WHICH
C        THE STATISTICS ARE COMPUTED.  Y IS SORTED, BUT RESTORED
C        TO ITS ORIGINAL ORDER AFTERWARDS.
C
C
C     INITIALIZE NAME VECTORS
C
      DATA  NMSUB(1), NMSUB(2), NMSUB(3), NMSUB(4), NMSUB(5), NMSUB(6)
     +   /       'S',      'T',      'A',      'T',      ' ',      ' '/
C
      DATA ALPHA /0.25E0/
      DATA IDP /4/
      DATA IINT /2/
      DATA WTS /.FALSE./
      DATA STACK /.TRUE./
C
C     CHECK FOR ERRORS IN THE INPUT PARAMETERS
C
      CALL STATER(NMSUB, WT, N, LDSTAK, WTS, NNZW, STACK, IERR)
      IF (IERR.NE.0) THEN
C
C     PRINT ERROR MESSAGE.
C
         CALL IPRINT(IPRT)
         WRITE (IPRT,1000)
         RETURN
      END IF
C
C     SET UP FRAMEWORK AREA.
C
      CALL STKSET (LDSTAK, IDP)
      NALL0 = STKST(1)
C
C     SET UP LSORT, THE PERMUTATION VECTOR.
C
      LSORT = STKGET(N,IINT)
      CALL GENI(ISTAK(LSORT), N, 1, 1)
C
C     SORT THE VECTOR Y.
C
      CALL SRTIR(ISTAK(LSORT), N, Y)
C
C     COMPUTE THE STATISTICS WHICH USE A SORTED ARRAY.
C
      CALL STAT1(Y, N, STS(5), STS(34), STS(35), STS(6),
     +   STS(11), 10, 0.0E0, 0.0E0, STS(44))
C
C     CALCULATE SUMS OF THE SORTED ARRAY.
C
      CALL SUMBS(Y, N, 1, MID, N)
      CALL SUMSS(Y, N, 1, MID, N, STS(38), STS(39), STS(42),
     +   STS(3))
      STS(4) = STS(3)
      CALL SUMTS(Y, N, ALPHA, SUMT1, STS(7))
      STS(8) = STS(7)
      CALL SUMDS(Y, N, 1, MID, N, STS(3), SUMDA, SUMD2, SUMD3, SUMD4)
C
C     RESTORE THE VECTOR Y TO ITS ORIGINAL ORDER.
C
      CALL SRTRI(Y, N, ISTAK(LSORT))
C
C     COMPUTE REST OF STATISTICS.
C
      CALL SUMID(Y, N, STS(3), SUMDI)
      CALL STAT2(Y, N, STS, SUMDA, SUMDI, SUMD2, SUMD3, SUMD4)
      CALL SUMOT(STS, N, N, WTS)
C
C     RETURN THE VECTOR LSORT.
C
      CALL STKCLR(NALL0)
      RETURN
C
C     FORMAT STATEMENTS.
C
 1000 FORMAT (/42H THE CORRECT FORM OF THE CALL STATEMENT IS//
     +   '       CALL STAT (Y, N, LDSTAK)')
      END
