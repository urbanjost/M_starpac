*AMEER
      SUBROUTINE AMEER(NMSUB, N, NPAR, NPARE, LDSTAK, LDSMIN,
     +  STP, LSTP, SCALE, LSCALE, IVCV, SAVE, MSPEC, NFAC)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE ERROR CHECKING ROUTINE FOR NONLINEAR LEAST SQUARES
C     ESTIMATION ROUTINES.
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  DECEMBER 2, 1985
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   IVCV,LDSMIN,LDSTAK,LSCALE,LSTP,N,NFAC,NPAR,NPARE
      LOGICAL
     +   SAVE
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   SCALE(*),STP(*)
      INTEGER
     +   MSPEC(4,*)
      CHARACTER
     +   NMSUB(6)*1
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      INTEGER
     +   I,NP,NV
      LOGICAL
     +   HEAD
C
C  LOCAL ARRAYS
      LOGICAL
     +   ERROR(20)
      CHARACTER
     +   LIVCV(8)*1,LLDS(8)*1,LMSPEC(8)*1,LN(8)*1,LNFAC(8)*1,
     +   LNPAR(8)*1,LNPARE(8)*1,LONE(8)*1,LSCL(8)*1,LSTEP(8)*1,
     +   LZERO(8)*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL EIAGE,EISEQ,EISGE,ERVGT
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     LOGICAL ERROR(20)
C        VALUE(S) INDICATING WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT
C        (FALSE).
C     LOGICAL HEAD
C        THE VARIABLE USED TO INDICATE WHETHER A HEADING IS TO BE
C        PRINTED DURING A GIVEN CALL TO THE ITERATION REPORT (TRUE)
C        OR NOT (FALSE).
C     INTEGER I
C        AN INDEX VARIABLE.
C     INTEGER IERR
C        THE VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST.
C        IF IERR .EQ. 0, NO ERRORS WERE DETECTED.
C        IF IERR .GE. 1, ERRORS WERE DETECTED.
C     INTEGER IVCV
C        THE FIRST DIMENSION OF MATRIX VCV.
C     INTEGER LDSMIN
C        THE MINIMUM LENGTH ALLOWED FOR THE ARRAY DSTAK.
C     INTEGER LDSTAK
C        THE LENGTH OF THE ARRAY DSTAK.
C     CHARACTER*1 LIVCV(8), LLDS(8), LMSPEC(8), LN(8), LNFAC(8),
C    *   LNPAR(8), LNPARE(8), LONE(8), LSCL(8), LSTEP(8), LZERO(8)
C        THE ARRAY(S) CONTAINING THE NAME(S) OF INPUT PARAMETER(S)
C        CHECKED FOR ERRORS.
C     INTEGER LSCALE
C        THE DIMENSION OF VECTOR SCALE.
C     INTEGER LSTP
C        THE DIMENSION OF VECTOR STP.
C     INTEGER MSPEC(4,*)
C        INTEGER MSPEC(4,NFAC)
C        THE ARRAY CONTAINING THE VALUES OF P, D, Q, AND S FOR EACH FACT
C     INTEGER N
C        THE NUMBER OF OBSERVATIONS.
C     INTEGER NFAC
C        THE NUMBER OF FACTORS IN THE MODEL.
C     CHARACTER*1 NMSUB(6)
C        THE NAME OF THE ROUTINE CALLING THE ERROR CHECKING ROUTINE
C     INTEGER NP
C        THE NUMBER OF PARAMETERS SPECIFIED BY MSPEC.
C     INTEGER NPAR
C        THE NUMBER OF PARAMETERS IN THE MODEL.
C     INTEGER NPARE
C        THE NUMBER OF PARAMETERS TO BE OPTIMIZED.
C     INTEGER NV
C        *
C     LOGICAL SAVE
C        THE VARIABLE USED TO INDICATE WHETHER ANY RESULTS OTHER THAN
C        THE RESIDUALS AND PARAMETERS ARE TO BE SAVED (TRUE) OR NOT
C        (FALSE).
C     DOUBLE PRECISION SCALE(LSCALE)
C        THE TYPICAL SIZE OF THE PARAMETERS.
C     DOUBLE PRECISION STP(LSTP)
C        THE STEP SIZE ARRAY.
C
C     SET UP NAME ARRAYS
C
      DATA LIVCV(1), LIVCV(2), LIVCV(3), LIVCV(4), LIVCV(5),
     +   LIVCV(6), LIVCV(7), LIVCV(8) /'I','V','C','V',' ',' ',' ',' '/
      DATA LLDS(1), LLDS(2), LLDS(3), LLDS(4), LLDS(5), LLDS(6),
     +   LLDS(7), LLDS(8) /'L','D','S','T','A','K',' ',' '/
      DATA LMSPEC(1), LMSPEC(2), LMSPEC(3), LMSPEC(4), LMSPEC(5),
     +   LMSPEC(6), LMSPEC(7), LMSPEC(8)
     +  /'M','S','P','C',' ',' ',' ',' '/
      DATA LN(1), LN(2), LN(3), LN(4), LN(5), LN(6), LN(7), LN(8) /'N',
     +   ' ',' ',' ',' ',' ',' ',' '/
      DATA LNFAC(1), LNFAC(2), LNFAC(3), LNFAC(4), LNFAC(5),
     +   LNFAC(6), LNFAC(7), LNFAC(8) /'N','F','A','C',' ',' ',' ',' '/
      DATA LNPAR(1), LNPAR(2), LNPAR(3), LNPAR(4), LNPAR(5),
     +   LNPAR(6), LNPAR(7), LNPAR(8) /'N','P','A','R',' ',' ',' ',
     +   ' '/
      DATA LNPARE(1), LNPARE(2), LNPARE(3), LNPARE(4), LNPARE(5),
     +   LNPARE(6), LNPARE(7), LNPARE(8) /'N','P','A','R','E',' ',' ',
     +   ' '/
      DATA LONE(1), LONE(2), LONE(3), LONE(4), LONE(5),
     +   LONE(6), LONE(7), LONE(8) /'1',' ',' ',' ',' ',' ',' ',' '/
      DATA LSCL(1), LSCL(2), LSCL(3), LSCL(4), LSCL(5),
     +   LSCL(6), LSCL(7), LSCL(8) /'S','C','A','L','E',' ',' ',
     +   ' '/
      DATA LSTEP(1), LSTEP(2), LSTEP(3), LSTEP(4), LSTEP(5),
     +   LSTEP(6), LSTEP(7), LSTEP(8) /'S','T','P',' ',' ',' ',' ',' '/
      DATA LZERO(1), LZERO(2), LZERO(3), LZERO(4), LZERO(5),
     +   LZERO(6), LZERO(7), LZERO(8) /'Z','E','R','O',' ',' ',' ',' '/
C
C     ERROR CHECKING
C
      DO 10 I=1,20
         ERROR(I) = .FALSE.
   10 CONTINUE
C
      IERR = 0
      HEAD = .TRUE.
C
      CALL EISGE(NMSUB, LN, N, 1, 2, HEAD, ERROR(1), LONE)
C
      CALL EISGE(NMSUB, LNFAC, NFAC, 1, 2, HEAD, ERROR(2), LONE)
C
      IF (.NOT. ERROR(2))
     +  CALL EIAGE(NMSUB, LMSPEC, MSPEC, 4, NFAC, 4, 0, 0, HEAD, 1, NV,
     +  ERROR(3), LMSPEC)
C
      IF ((.NOT. ERROR(2)) .AND. (.NOT. ERROR(3))) THEN
        NP = 1
         DO 20 I = 1, NFAC
          NP = NP + MSPEC(1,I) + MSPEC(3,I)
   20   CONTINUE
        CALL EISEQ(NMSUB, LNPAR, NPAR, NP, 1, HEAD, ERROR(4), LNPAR)
C
        IF (.NOT.ERROR(4)) THEN
          CALL EISGE(NMSUB, LNPARE, NPARE, 1, 2, HEAD, ERROR(5), LONE)
          CALL ERVGT(NMSUB, LSTEP, STP, LSTP, 0.0D0, 0, HEAD, 6, NV,
     +      ERROR(8), LZERO)
          CALL ERVGT(NMSUB, LSCL, SCALE, LSCALE, 0.0D0, 0, HEAD, 6, NV,
     +      ERROR(12), LZERO)
          IF (SAVE .AND. (.NOT.ERROR(5)))
     +      CALL EISGE(NMSUB, LIVCV, IVCV, NPARE, 3, HEAD, ERROR(15),
     +      LNPARE)
        END IF
      END IF
C
      IF ((.NOT.ERROR(1)) .AND. (.NOT.ERROR(2)) .AND. (.NOT.ERROR(3))
     +   .AND. (.NOT.ERROR(4)) .AND. (.NOT.ERROR(5)))
     +   CALL EISGE(NMSUB, LLDS, LDSTAK, LDSMIN, 9, HEAD, ERROR(6),
     +   LLDS)
C
      DO 30 I=1,20
         IF (ERROR(I)) GO TO 40
   30 CONTINUE
      RETURN
C
   40 CONTINUE
      IERR = 1
      RETURN
C
      END
