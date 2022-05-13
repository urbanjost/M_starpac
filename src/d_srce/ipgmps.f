*IPGMPS
      SUBROUTINE IPGMPS (PER, FREQ, NF, N, LDSTAK, PERI, NPRT)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE USER CALLABLE ROUTINE FOR COMPUTING
C     THE INTEGRATED PERIODOGRAM OF A SERIES (LONG CALL).
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  DECEMBER 7, 1981
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   LDSTAK,N,NF,NPRT
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   FREQ(*),PER(*),PERI(*)
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
     +   IPRT,ISYM,LDSMIN,LPCV,NALL0,XAXIS,YAXIS
      LOGICAL
     +   ERR01,ERR02,ERR03,HEAD
C
C  LOCAL ARRAYS
      DOUBLE PRECISION
     +   RSTAK(12)
      INTEGER
     +   ISTAK(12)
      CHARACTER
     +   LLDS(8)*1,LN(8)*1,LNF(8)*1,NMSUB(6)*1
C
C  EXTERNAL FUNCTIONS
      INTEGER
     +   STKGET,STKST
      EXTERNAL STKGET,STKST
C
C  EXTERNAL SUBROUTINES
      EXTERNAL EISGE,IPGDV,IPRINT,LDSCMP,STKCLR,STKSET
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
C     DOUBLE PRECISION DSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     LOGICAL ERR01, ERR02, ERR03
C        VALUE(S) INDICATING WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT
C        (FALSE).
C     DOUBLE PRECISION FREQ(NF)
C        THE ARRAY IN WHICH THE FREQUENCIES CORRESPONDING TO THE
C        INTEGRATED SPECTRUM VALUES ARE STORED.
C     LOGICAL HEAD
C        A VARIABLE USED TO INDICATE WHETHER A HEADING IS NEEDED FOR
C        ERROR MESSAGES (TRUE) OR NOT (FALSE).
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED
C        IF IERR .EQ. 1, ERRORS HAVE BEEN DETECTED
C     INTEGER IPRT
C        THE LOGICAL UNIT USED FOR PRINTED OUTPUT.
C     INTEGER ISTAK(12)
C        THE INTEGER VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER LDSMIN
C        THE MINIMUM LENGTH ALLOWED FOR DSTAK.
C     INTEGER LDSTAK
C        THE LENGTH OF THE VECTOR DSTAK IN COMMON CSTAK.
C     CHARACTER*1 LLDS(8), LN(8), LNF(8)
C        THE ARRAY(S) CONTAINING THE NAME(S) OF THE PARAMETER(S) CHECKED
C        FOR ERRORS.
C     INTEGER N
C        THE ACTUAL NUMBER OF OBSERVATIONS IN THE SERIES FROM WHICH
C        THE PERIODOGRAM WAS COMPUTED.
C     INTEGER NALL0
C        THE NUMBER OF OUTSTANDING ALLOCATIONS OF THE STACK AT THE
C        TIME OF THIS CALL.
C     INTEGER NF
C        THE NUMBER OF FREQUENCIES AT WHICH THE PERIODGRAM IS
C        COMPUTED.
C     CHARACTER*1 NMSUB(6)
C        THE ARRAY CONTAINING THE NAME OF THIS SUBROUTINE.
C     INTEGER NPRT
C        THE VARIABLE CONTROLING PRINTED OUTPUT, WHERE
C        IF NPRT .EQ.  0, THE OUTPUT IS SUPPRESSED,
C        IF NPRT .GE.  1, THE OUTPUT CONSISTS OF A PAGE PLOT.
C     DOUBLE PRECISION PER(NF)
C        THE INTEGRATED PERIODOGRAM.
C     DOUBLE PRECISION PERI(NF)
C        THE VECTOR IN WHICH THE INTEGRATED PERIODOGRAM IS STORED.
C     DOUBLE PRECISION RSTAK(12)
C        THE DOUBLE PRECISION VERSION OF THE /CSTAK/ WORK AREA.
C     INTEGER XAXIS
C        THE STARTING LOCATION IN THE STACK FOR
C        THE ARRAY IN WHICH THE X AXIS VALUES TO BE PLOTTED ARE STORED.
C     INTEGER YAXIS
C        THE STARTING LOCATION IN THE STACK FOR
C        THE ARRAY IN WHICH THE Y AXIS VALUES TO BE PLOTTED ARE STORED.
C
C     SET UP NAME ARRAYS
C
      DATA
     +  NMSUB(1),  NMSUB(2),  NMSUB(3),  NMSUB(4),  NMSUB(5),  NMSUB(6)
     + /     'I',       'P',       'G',       'M',       'P',       'S'/
      DATA
     + LLDS(1), LLDS(2), LLDS(3), LLDS(4), LLDS(5),
     +  LLDS(6), LLDS(7), LLDS(8) /'L','D','S','T','A','K',' ',' '/
      DATA
     + LN(1), LN(2), LN(3), LN(4), LN(5), LN(6), LN(7), LN(8)
     + /'N',' ',' ',' ',' ',' ',' ',' '/
      DATA
     + LNF(1), LNF(2), LNF(3), LNF(4), LNF(5), LNF(6), LNF(7), LNF(8)
     + /'N','F',' ',' ',' ',' ',' ',' '/
C
C     SET UP FOR ERROR CHECKING
C
      IERR = 0
      HEAD = .TRUE.
C
C     CALL ERROR CHECKING ROUTINES
C
      CALL EISGE(NMSUB, LN, N, 17, 1, HEAD, ERR01, LN)
C
      CALL EISGE(NMSUB, LNF, NF, (N+2)/2, 1, HEAD, ERR02, LNF)
C
      IF (ERR01) GO TO 5
C
      IF (NPRT .EQ. 0) THEN
        LDSMIN = 0
      ELSE
        CALL LDSCMP(3, 0, NF+103, 0, 0, 0, 'D', 2*NF+206, LDSMIN)
      END IF
C
      CALL EISGE(NMSUB, LLDS, LDSTAK, LDSMIN, 9, HEAD, ERR03, LLDS)
C
      IF (ERR02 .OR. ERR03) GO TO 5
      GO TO 10
C
    5 CONTINUE
      IERR = 1
      CALL IPRINT (IPRT)
      WRITE (IPRT, 1000)
      RETURN
C
   10 CONTINUE
C
C     SET THE SIZE OF THE WORK AREA
C
      CALL STKSET(LDSTAK, 4)
C
C     SET THE NUMBER OF OUTSTANDING ALLOCATIONS.
C
      NALL0 = STKST(1)
C
C     SET VARIOUS PROGRAM PARAMETERS.
C
      LPCV = NF + 103
C
C     SUBDIVIDE THE STACK.
C
      IF (NPRT .EQ. 0) THEN
         ISYM = 1
         XAXIS = 1
         YAXIS = 1
      ELSE
         ISYM = STKGET(LPCV, 2)
         XAXIS = STKGET(LPCV, 4)
         YAXIS = STKGET(LPCV, 4)
      END IF
C
C     CALL THE MAIN DRIVER FOR COMPUTING (AND PLOTTING) THE INTEGRATED
C     PERIODOGRAM.
C
      CALL IPGDV (PER, NF, N, PERI, FREQ, RSTAK(XAXIS),
     +   RSTAK(YAXIS), ISTAK(ISYM), LPCV, NPRT)
C
      CALL STKCLR(NALL0)
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (/42H THE CORRECT FORM OF THE CALL STATEMENT IS//
     +  '       CALL IPGMPS (PER, FREQ, NF, N, LDSTAK, PERI, NPRT)')
      END
