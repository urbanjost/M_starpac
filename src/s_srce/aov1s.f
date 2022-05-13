*AOV1S
      SUBROUTINE AOV1S(Y, TAG, N, LDSTAK, NPRT, GSTAT, IGSTAT, NG)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE -
C     1. CALLS OTHER ROUTINES TO CHECK THE INPUT PARAMETERS
C     2. SETS UP NEEDED STORAGE LOCATIONS AND
C     3. CALLS AOV1MN TO COMPUTE A COMPREHENSIVE SET OF RESULTS FOR A
C         ONEWAY ANALYSIS OF VARIANCE WITH OPTIONAL OUTPUT.
C
C     WRITTEN BY -
C       LINDA MITCHELL
C       STATISTICAL ENGINEERING DIVISION
C       NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  MAY 17, 1982
C                       BASED ON EARLIER VERSION BY J. R. DONALDSON
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   IGSTAT,LDSTAK,N,NG,NPRT
C
C  ARRAY ARGUMENTS
      REAL
     +   GSTAT(*),TAG(*),Y(*)
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
     +   B10,GPMAX,GPMIN,IFP,INDEX,INT,IPRT,ITEMP,NALL0,NZTAGS,
     +   RANKS,SRANK
C
C  LOCAL ARRAYS
      REAL
     +   RSTAK(12)
      INTEGER
     +   ISTAK(12)
      CHARACTER
     +   NMSUB(6)*1
C
C  EXTERNAL FUNCTIONS
      INTEGER
     +   STKGET
      EXTERNAL STKGET
C
C  EXTERNAL SUBROUTINES
      EXTERNAL AOV1ER,AOV1HD,AOV1MN,IPRINT,STKCLR
C
C  COMMON BLOCKS
      COMMON /CSTAK/DSTAK
      COMMON /ERRCHK/IERR
C
C  EQUIVALENCES
      EQUIVALENCE (DSTAK(1),RSTAK(1))
      EQUIVALENCE (DSTAK(1),ISTAK(1))
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     INTEGER B10
C        STARTING LOCATION IN THE STACK AREA FOR B10
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER GPMAX
C        THE STARTING LOCATION IN THE STACK AREA OF MAXIMUM
C        OBSERVATION
C     INTEGER GPMIN
C        THE STARTING LOCATION IN THE STACK AREA OF THE MINUMUM
C        OBSERVATION
C     REAL GSTAT(IGSTAT,4)
C        THE GROUP STATISTICS.  COLUMNS CORRESPOND TO THE TAG
C        VALUE, SAMPLE SIZE, GROUP MEAN, AND GROUP STANDARD DEVIATION.
C     INTEGER IERR
C        A COMMON VARIABLE USED AS A FLAG INDICATING WHETHER THERE
C        ARE ANY ERRORS, IF = 0 THEN NO ERRORS
C     INTEGER IFP
C        AN INDICATOR FOR STACK ALLOCATION TYPE, WHERE IFP=3 INDICATES
C        SINGLE PRECISION AND IFP=4 INDICATES DOUBLE PRECISION.
C     INTEGER IGSTAT
C        THE FIRST DIMENSION OF GSTAT.
C     INTEGER INDEX
C        THE STARTING LOCATION IN THE STACK ARRAY OF THE INDEX FOR
C        THE SORTED TAGS
C     INTEGER INT
C        FRAMEWORK CODE VALUE FOR INTEGER NUMBERS
C     INTEGER IPRT
C        THE OUTPUT LOGICAL UNIT NUMBER
C     INTEGER ISTAK(12)
C        THE INTEGER VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER ITEMP
C        STARTING LOCATION IN THE STACK FOR THE
C        TEMPORARY STORAGE ARRAY
C     INTEGER LDSTAK
C         SIZE OF THE STACK AREA ALLOCATED IN THE USERS MAIN PROGRAM
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS TO BE ANALYZED
C     INTEGER NALL0
C        THE NUMBER OF ALLOCATIONS OUTSTANDING AT THE TIME THAT THIS
C        ROUTINE WAS CALLED.
C     INTEGER NG
C        THE COMPUTED NUMBER OF GROUPS WITH
C        DIFFERENT POSITIVE TAG VALUES
C     CHARACTER*1 NMSUB(6)
C        SUBROUTINE NAME
C     INTEGER NPRT
C        THE VARIABLE CONTROLLING AUTOMATIC PRINTOUT
C        IF =0, PRINTOUT IS SUPRESSED
C        OTHERWISE PRINTOUT IS PROVIDED
C     INTEGER NZTAGS
C        THE NUMBER OF OBSERVATIONS WITH POSITIVE NON-ZERO WIEGHTS
C     INTEGER RANKS
C        THE STARTING LOCATION IN STACK AREA FOR THE RANKS OF Y
C     REAL RSTAK(12)
C        THE REAL VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER SRANK
C        THE STARTING LOCATION IN STACK FOR THE SUM OF RANKS
C     REAL TAG(N)
C        THE VECTOR OF TAG VALUES
C     REAL Y(N)
C        THE VECTOR OF OBSERVATIONS
C
      DATA   NMSUB(1), NMSUB(2), NMSUB(3), NMSUB(4), NMSUB(5), NMSUB(6)
     +     /      'A',      'O',      'V',      '1',      'S',      ' '/
C
C     SET UP FRAMEWORK VARIABLES FOR NUMBER TYPES
C
      INT = 2
      IFP = 3
C
      CALL IPRINT(IPRT)
C
C     CHECK FOR ERRORS IN PARAMETERS, INITIALIZE STACK, AND SET
C     NALL0.
C
      CALL AOV1ER(Y, TAG, N, IGSTAT, NZTAGS, NG, LDSTAK, NMSUB, INDEX,
     +   0, NALL0)
C
      IF (IERR.EQ.0) GO TO 20
C
C     PRINT CORRECT FORM OF CALL STATEMENT AND RETURN TO CALLER
C
      IERR = 1
      WRITE (IPRT,1000)
      RETURN
C
C     PRINT HEADING IF DESIRED
C
   20 IF (NPRT.EQ.0) GO TO 30
      CALL AOV1HD(IPRT)
C
C     SET UP ADDITIONAL WORK VECTORS FOR AOV1MN AS CALLED FROM AOV1S
C
   30 SRANK = STKGET(NG,IFP)
      GPMIN = STKGET(NG,IFP)
      GPMAX = STKGET(NG,IFP)
      B10 = STKGET(NG,IFP)
      RANKS = STKGET(NZTAGS,IFP)
      ITEMP = STKGET(NZTAGS,INT)
C
      CALL AOV1MN(Y, TAG, N,
     +            GSTAT(1), GSTAT(IGSTAT+1),
     +            GSTAT(2*IGSTAT+1), GSTAT(3*IGSTAT+1),
     +            NPRT, ISTAK(INDEX), RSTAK(SRANK), RSTAK(GPMIN),
     +            RSTAK(GPMAX), RSTAK(B10), RSTAK(RANKS),
     +            ISTAK(ITEMP), NG, NZTAGS)
C
C     RELEASE THE STACK AREA
C
      CALL STKCLR(NALL0)
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT(/42H THE CORRECT FORM OF THE CALL STATEMENT IS//
     + '       CALL AOV1S (Y, TAG, N, LDSTAK, NPRT, GSTAT, IGSTAT, NG)')
      END
