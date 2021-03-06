*STATW
      SUBROUTINE STATW(Y, WT, N, LDSTAK)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS SUBROUTINE COMPUTES 53 DIFFERENT STATISTICS FOR A
C     VECTOR Y, WITH WEIGHTS SPECIFIED.  ONE PAGE OF AUTOMATIC
C     PRINTOUT IS PRODUCED.
C
C     WRITTEN BY - JANET R. DONALDSON, JOHN E. KOONTZ
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS
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
      DOUBLE PRECISION
     +   WT(*),Y(*)
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  ARRAYS IN COMMON
      DOUBLE PRECISION DSTAK(12)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   ALPHA,SUM1,SUMD2,SUMD3,SUMD4,SUMDA,SUMDI,SUMT1,SUMW,SUMWD2,
     +   SUMWT1
      INTEGER
     +   IDP,IINT,IPRT,LSORT,MID,NALL0,NNZW
      LOGICAL
     +   STACK,WTS
C
C  LOCAL ARRAYS
      DOUBLE PRECISION
     +   STS(53)
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
      EXTERNAL GENI,IPRINT,SRTIRR,SRTRRI,STAT1W,STAT2W,STATER,STKCLR,
     +   STKSET,SUMBS,SUMIDW,SUMOT,SUMWDS,SUMWSS,SUMWTS
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
C     DOUBLE PRECISION ALPHA
C        THE PERCENTAGE OF POINTS TO BE TRIMMED FROM EITHER END OF
C        Y IN CALCULATING THE TRIMMED MEANS.
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER IDP
C        FRAMEWORK CODE VALUE FOR DOUBLE PRECISION NUMBERS.
C     INTEGER IERR
C        THE CODE INDICATING WHETHER OR NOT AN ERROR HAS
C        BEEN DISCOVERED.  0 MEANS NO ERROR, NOT 0 MEANS
C        SOME ERROR EXISTS.
C     INTEGER IINT
C        THE CODE VALUE FOR INTEGER FOR FRAMEWORK.
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
C        THE INDEX OF A ZERO ELEMENT IN THE SORTED Y, OR OF THE
C        ELEMENT CLOSEST TO ZERO.
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
C     DOUBLE PRECISION STS(53)
C        THE VECTOR OF THE 53 STATISTICS COMPUTED.
C     DOUBLE PRECISION SUMDA
C        THE SUM OF THE ABSOLUTE DIFFERENCES FROM THE MEAN.
C     DOUBLE PRECISION SUMDI
C        THE SUM OF THE PRODUCTS OF THE INDICES AND THE
C        DIFFERENCES.
C     DOUBLE PRECISION SUMD2
C        THE SUM OF THE SQUARES OF THE DIFFERENCES.
C     DOUBLE PRECISION SUMD3
C        THE SUM OF THE CUBES OF THE DIFFERENCES.
C     DOUBLE PRECISION SUMD4
C        THE SUM OF THE 4TH POWERS OF THE DIFFERENCES.
C     DOUBLE PRECISION SUMT1
C        THE SUM OF THE ALPHA TRIMMED ARRAY Y.
C     DOUBLE PRECISION SUMW
C        THE SUM OF THE WEIGHTS VECTOR WT.
C     DOUBLE PRECISION SUMWD2
C        THE WEIGHTED SUM OF THE SQUARES OF THE DIFFERENCES.
C     DOUBLE PRECISION SUMWT1
C        THE WEIGHTED SUM OF THE ALPHA TRIMMED ARRAY.
C     DOUBLE PRECISION SUM1
C        THE UNWEIGHTED SUM OF THE ELEMENTS OF Y.
C     DOUBLE PRECISION WT(N)
C        INPUT PARAMETER.  THE WEIGHTS VECTOR.
C     LOGICAL WTS
C        A FLAG INDICATING WHETHER THERE ARE WEIGHTS (TRUE)
C        OR NOT (FALSE).
C     DOUBLE PRECISION Y(N)
C        INPUT PARAMETER.  THE VECTOR OF DATA POINTS ON WHICH
C        THE STATISTICS ARE COMPUTED.  Y IS SORTED, BUT RESTORED
C        TO ITS ORIGINAL ORDER AFTERWARDS.
C
C
C     INITIALIZE NAME VECTORS
C
      DATA  NMSUB(1), NMSUB(2), NMSUB(3), NMSUB(4), NMSUB(5), NMSUB(6)
     +   /       'S',      'T',      'A',      'T',      'W',      ' '/
C
      DATA ALPHA /0.25D0/
      DATA IDP /4/
      DATA IINT /2/
      DATA WTS /.TRUE./
      DATA STACK /.TRUE./
C
C     CHECK FOR ERRORS IN THE INPUT PARAMETERS.
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
      CALL SRTIRR(ISTAK(LSORT), WT, N, Y)
C
C     COMPUTE THE STATISTICS WHICH USE A SORTED ARRAY.
C
      CALL STAT1W(Y, WT, N, STS(5), STS(34), STS(35), STS(6),
     +   STS(11), 10, 0.0D0, 0.0D0, STS(44), NNZW)
C
C     COMPUTED VARIOUS SUMS IN THE SORTED ARRAY Y.
C
      CALL SUMBS(Y, N, 1, MID, N)
      CALL SUMWSS(Y, WT, N, 1, MID, N, NNZW, SUM1, STS(38), STS(39),
     +   STS(42), SUMW, STS(3), STS(4))
      CALL SUMWTS(Y, WT, N, NNZW, ALPHA, SUMT1, SUMWT1, STS(7),
     +   STS(8))
      CALL SUMWDS(Y, WT, N, 1, MID, N, STS(4), SUMDA, SUMWD2, SUMD2,
     +   SUMD3, SUMD4)
C
C     RESTORE THE VECTOR Y TO ITS ORIGINAL ORDER.
C
      CALL SRTRRI(Y, WT, N, ISTAK(LSORT))
C
C     COMPUTE REST OF STATISTICS.
C
      CALL SUMIDW(Y, WT, N, STS(4), SUMDI)
      CALL STAT2W(Y, WT, N, NNZW, STS, SUMDA, SUMDI, SUMWD2, SUMD2,
     +   SUMD3, SUMD4, SUMW)
      CALL SUMOT(STS, N, NNZW, WTS)
C
C     RETURN THE VECTOR LSORT.
C
      CALL STKCLR(NALL0)
      RETURN
C
C     FORMAT STATEMENTS.
C
 1000 FORMAT (/42H THE CORRECT FORM OF THE CALL STATEMENT IS//
     +   '       CALL STATW (Y, WT, N, LDSTAK)')
      END
