*LLER
      SUBROUTINE LLER(NMSUB, IXM, IVCV, N, NPAR, LPAR, LDSTAK, WT, LNWT,
     +   WEIGHT, NNZW, IFIT, SAVE)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE ERROR CHECKING ROUTINE FOR THE LINEAR LEAST
C     SQUARES LLSTING ROUTINES.
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
     +   IFIT,IVCV,IXM,LDSTAK,LNWT,LPAR,N,NNZW,NPAR
      LOGICAL
     +   SAVE,WEIGHT
C
C  ARRAY ARGUMENTS
      REAL
     +   WT(*)
      CHARACTER
     +   NMSUB(6)*1
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      INTEGER
     +   I,LDSMIN,NZW
      LOGICAL
     +   HEAD
C
C  LOCAL ARRAYS
      LOGICAL
     +   ERROR(10)
      CHARACTER
     +   LIVCV(8)*1,LIXM(8)*1,LLDS(8)*1,LLPAR(8)*1,LN(8)*1,
     +   LN1(8)*1,LNC(8)*1,LNDEG(8)*1,LNDEG1(8)*1,LNPAR(8)*1,
     +   LONE(8)*1,LWT(8)*1,LZERO(8)*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL EISGE,EISII,ERVWT,LDSCMP
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
C     INTEGER I
C        AN INDEX.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .GE. 1, ERRORS WERE DETECTED.
C     INTEGER IFIT
C        THE INDICATOR VALUE DESIGNATING WHETHER THE LLS IS OF A
C        GENERAL MODEL (IFIT=3) OR A POLYNOMIAL MODEL (IFIT=1).
C     INTEGER IVCV
C        THE FIRST DIMENSION OF THE VARIANCE COVARIANCE MATRIX VCV.
C     INTEGER IXM
C        THE FIRST DIMENSION OF THE INDEPENDENT VARIABLE ARRAY.
C     INTEGER LDSTAK
C        THE LENGTH OF THE ARRAY DSTAK.
C     INTEGER LDSMIN
C        THE MINIMUM LENGTH ALLOWED FOR THE ARRAY DSTAK.
C     CHARACTER*1 LIVCV(8), LIXM(8), LLPAR(8), LLDS(8), LN(8), LNC(8),
C    *   LNDEG(8), LNDEG1(8), LNPAR(8), LN1(8), LONE(8), LWT(8),
C    *   LZERO(8)
C        THE ARRAY(S) CONTAINING THE NAME(S) OF INPUT PARAMETER(S)
C        CHECKED FOR ERRORS.
C     INTEGER LPAR
C        THE ACTUAL LENGTH OF THE VECTOR P.
C     INTEGER LNWT
C        THE ACTUAL LENGTH OF THE VECTOR WT.
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS.
C     INTEGER NPAR
C        THE NUMBER OF UNKNOWN PARAMETERS IN THE MODEL.
C     CHARACTER*1 NMSUB(6)
C        THE NAME OF THE SUBROUTINE CALLING THE ERROR CHECKING
C        SUBROUTINES.
C     INTEGER NNZW
C        THE NUMBER OF NON ZERO WEIGHTS.
C     INTEGER NZW
C        THE NUMBER OF ZERO WEIGHTS.
C     LOGICAL SAVE
C        THE VARIABLE USED TO INDICATE WHETHER ANY RESULTS OTHER THAN
C        THE RESIDUALS ARE TO VE SAVED (TRUE) OR NOT (FALSE).
C     LOGICAL WEIGHT
C        THE VARIABLE USED TO INDICATE WHETHER WEIGHTED ANALYSIS IS TO
C        BE PERFORMED (TRUE) OR NOT (FALSE).
C     REAL WT(LNWT)
C        THE USER SUPPLIED WEIGHTS.
C
C     SET UP NAME ARRAYS
C
      DATA LIVCV(1), LIVCV(2), LIVCV(3), LIVCV(4), LIVCV(5), LIVCV(6),
     +   LIVCV(7), LIVCV(8) /'I','V','C','V',' ',' ',' ',' '/
      DATA LIXM(1), LIXM(2), LIXM(3), LIXM(4), LIXM(5), LIXM(6),
     +   LIXM(7), LIXM(8) /'I','X','M',' ',' ',' ',' ',' '/
      DATA LLDS(1), LLDS(2), LLDS(3), LLDS(4), LLDS(5), LLDS(6),
     +   LLDS(7), LLDS(8) /'L','D','S','T','A','K',' ',' '/
      DATA LLPAR(1), LLPAR(2), LLPAR(3), LLPAR(4), LLPAR(5), LLPAR(6),
     +   LLPAR(7), LLPAR(8) /'L','P','A','R',' ',' ',' ',' '/
      DATA LN(1), LN(2), LN(3), LN(4), LN(5), LN(6), LN(7), LN(8) /'N',
     +   ' ',' ',' ',' ',' ',' ',' '/
      DATA LNDEG(1), LNDEG(2), LNDEG(3), LNDEG(4), LNDEG(5), LNDEG(6),
     +   LNDEG(7), LNDEG(8) /'N','D','E','G',' ',' ',' ',' '/
      DATA LNDEG1(1), LNDEG1(2), LNDEG1(3), LNDEG1(4), LNDEG1(5),
     +   LNDEG1(6), LNDEG1(7), LNDEG1(8) /'N','D','E','G','+','1',
     +   ' ',' '/
      DATA LNPAR(1), LNPAR(2), LNPAR(3), LNPAR(4), LNPAR(5),
     +   LNPAR(6), LNPAR(7), LNPAR(8) /'N','P','A','R',' ',' ',' ',
     +   ' '/
      DATA LN1(1), LN1(2), LN1(3), LN1(4), LN1(5), LN1(6),
     +   LN1(7), LN1(8) /'N','-','1',' ',' ',' ',' ',' '/
      DATA LONE(1), LONE(2), LONE(3), LONE(4), LONE(5), LONE(6),
     +   LONE(7), LONE(8) /'O','N','E',' ',' ',' ',' ',' '/
      DATA LWT(1), LWT(2), LWT(3), LWT(4), LWT(5), LWT(6), LWT(7),
     +   LWT(8) /'W','T',' ',' ',' ',' ',' ',' '/
      DATA LZERO(1), LZERO(2), LZERO(3), LZERO(4), LZERO(5), LZERO(6),
     +   LZERO(7), LZERO(8) /'Z','E','R','O',' ',' ',' ',' '/
C
C     ERROR CHECKING
C
      IERR = 0
      HEAD = .TRUE.
C
      DO 10 I=1,10
         ERROR(I) = .FALSE.
   10 CONTINUE
C
      IF (IFIT.EQ.1) GO TO 30
C
      DO 20 I = 1, 8
         LNC(I) = LNPAR(I)
   20 CONTINUE
      GO TO 50
C
   30 CONTINUE
      DO 40 I = 1, 8
         LNC(I) = LNDEG1(I)
   40 CONTINUE
C
   50 CONTINUE
C
      CALL EISGE(NMSUB, LN, N, 1, 1, HEAD, ERROR(1), LN)
C
      IF (IFIT.EQ.3)
     +   CALL EISII(NMSUB, LNPAR, NPAR, 1, N, 1, HEAD, ERROR(2), LONE,
     +   LN)
      IF (IFIT.EQ.1)
     +   CALL EISII(NMSUB, LNDEG, NPAR-1, 0, N-1, 1, HEAD, ERROR(2),
     +      LZERO, LN1)
C
      CALL EISGE(NMSUB, LIXM, IXM, N, 3, HEAD, ERROR(4), LN)
C
      IF (SAVE .AND. (IFIT.EQ.1))
     +   CALL EISGE(NMSUB, LLPAR, LPAR, NPAR, 7, HEAD, ERROR(5), LNDEG1)
C
      IF (SAVE)
     +    CALL EISGE(NMSUB, LIVCV, IVCV, NPAR, 3, HEAD, ERROR(6), LNC)
C
      IF (ERROR(1) .OR. ERROR(2) .OR. ERROR(3)) GO TO 70
C
      NNZW = N
      IF (WEIGHT) CALL ERVWT(NMSUB, LWT, WT, N, NPAR, HEAD, NNZW,
     +   NZW, 2, ERROR(8), LNC)
C
      CALL LDSCMP(15, 0, 0, 0, 0, 0, 'S',
     +            6*N + NPAR*(N+2*NPAR+5) + 1, LDSMIN)
C
      CALL EISGE(NMSUB, LLDS, LDSTAK, LDSMIN, 9, HEAD, ERROR(9), LLDS)
C
      DO 60 I=1,10
         IF (ERROR(I)) GO TO 70
   60 CONTINUE
      RETURN
C
   70 CONTINUE
      IERR = 1
      RETURN
C
      END
