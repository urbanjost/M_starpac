*STPER
      SUBROUTINE STPER(NMSUB, N, M, IXM, NPAR, LDSTAK, SCALE, LSCALE)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE ERROR CHECKING ROUTINE FOR STEP SIZE SELECTION
C     ROUTINES.
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  OCTOBER 3, 1983
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   IXM,LDSTAK,LSCALE,M,N,NPAR
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   SCALE(*)
      CHARACTER
     +   NMSUB(6)*1
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      INTEGER
     +   I,LDSMIN,NV
      LOGICAL
     +   HEAD
C
C  LOCAL ARRAYS
      LOGICAL
     +   ERROR(10)
      CHARACTER
     +   LIXM(8)*1,LLDS(8)*1,LM(8)*1,LN(8)*1,LNPAR(8)*1,
     +   LONE(8)*1,LSCL(8)*1,LZERO(8)*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL EISGE,ERVGT,LDSCMP
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MAX
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     LOGICAL ERROR(10)
C        VALUE(S) INDICATING WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT
C        (FALSE).
C     LOGICAL HEAD
C        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
C        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
C        OF HEAD WILL BE CHANGED TO FALSE.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .EQ. 1, ERRORS WERE DETECTED.
C     INTEGER IXM
C        THE FIRST DIMENSION OF THE INDEPENDENT VARIABLE ARRAY XM.
C     INTEGER LDSMIN
C        THE MINIMUM LENGTH ALLOWED FOR THE ARRAY DSTAK.
C     INTEGER LDSTAK
C        THE LENGTH OF THE ARRAY DSTAK.
C     CHARACTER*1 LIXM(8), LLDS(8), LM(8), LN(8), LNPAR(8), LONE(8),
C    +            LSCL(8), LZERO(8)
C        THE ARRAY(S) CONTAINING THE NAME(S) OF INPUT PARAMETER(S)
C        CHECKED FOR ERRORS.
C     INTEGER LSCALE
C        THE LENGTH OF VECTOR SCALE.
C     INTEGER M
C        THE NUMBER OF INDEPENDENT VARIABLES.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS.
C     CHARACTER*1 NMSUB(6)
C        THE NAME OF THE SUBROUTINE CALLING THE ERROR CHECKING
C        SUBROUTINES.
C     INTEGER NPAR
C        THE NUMBER OF UNKNOWN PARAMETERS IN THE MODEL.
C     INTEGER NV
C        THE NUMBER OF VIOLATIONS FOUND.
C     DOUBLE PRECISION SCALE(LSCALE)
C        THE TYPICAL SIZE OF THE UNKNOWN PARAMETERS.
C
C     SET UP NAME ARRAYS
C
      DATA LIXM(1), LIXM(2), LIXM(3), LIXM(4), LIXM(5), LIXM(6),
     +   LIXM(7), LIXM(8) /'I','X','M',' ',' ',' ',' ',' '/
      DATA LLDS(1), LLDS(2), LLDS(3), LLDS(4), LLDS(5), LLDS(6),
     +   LLDS(7), LLDS(8) /'L','D','S','T','A','K',' ',' '/
      DATA LM(1), LM(2), LM(3), LM(4), LM(5), LM(6), LM(7), LM(8) /'M',
     +   ' ',' ',' ',' ',' ',' ',' '/
      DATA LN(1), LN(2), LN(3), LN(4), LN(5), LN(6), LN(7), LN(8) /'N',
     +   ' ',' ',' ',' ',' ',' ',' '/
      DATA LNPAR(1), LNPAR(2), LNPAR(3), LNPAR(4), LNPAR(5),
     +   LNPAR(6), LNPAR(7), LNPAR(8) /'N','P','A','R',' ',' ',' ',
     +   ' '/
      DATA LONE(1), LONE(2), LONE(3), LONE(4), LONE(5),
     +   LONE(6), LONE(7), LONE(8) /' ',' ','O','N','E',' ',' ',
     +   ' '/
      DATA LSCL(1), LSCL(2), LSCL(3), LSCL(4), LSCL(5),
     +   LSCL(6), LSCL(7), LSCL(8) /'S','C','A','L','E',' ',' ',
     +   ' '/
      DATA LZERO(1), LZERO(2), LZERO(3), LZERO(4), LZERO(5),
     +   LZERO(6), LZERO(7), LZERO(8) /'Z','E','R','O',' ',' ',' ',' '/
C
C     ERROR CHECKING
C
      DO 10 I=1,10
         ERROR(I) = .FALSE.
   10 CONTINUE
C
      IERR = 0
      HEAD = .TRUE.
C
      CALL EISGE(NMSUB, LN, N, 1, 2, HEAD, ERROR(1), LONE)
C
      CALL EISGE(NMSUB, LM, M, 1, 2, HEAD, ERROR(2), LONE)
C
      CALL EISGE(NMSUB, LIXM, IXM, N, 3, HEAD, ERROR(3), LN)
C
      CALL EISGE(NMSUB, LNPAR, NPAR, 1, 2, HEAD, ERROR(4), LONE)
C
      CALL LDSCMP(14, 0, 2*(N+NPAR), 0, 0, 0, 'D', 9*N + MAX(N,NPAR),
     +            LDSMIN)
C
      IF ((.NOT.ERROR(1)) .AND. (.NOT.ERROR(4)))
     +   CALL EISGE(NMSUB, LLDS, LDSTAK, LDSMIN, 9, HEAD, ERROR(5),
     +   LLDS)
C
      CALL ERVGT(NMSUB, LSCL, SCALE, LSCALE, 0.0D0, 0, HEAD, 6, NV,
     +   ERROR(9), LZERO)
C
C
      DO 20 I=1,10
         IF (ERROR(I)) GO TO 30
   20 CONTINUE
      RETURN
C
   30 CONTINUE
      IERR = 1
      RETURN
C
      END
