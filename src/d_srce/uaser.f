*UASER
      SUBROUTINE UASER(NMSUB, N, ACOV, IAR, PHI, LAGMAX, LAG, LACOV,
     +   NF, LDSTAK, LDSMIN, LYFFT, NFFT, OPTION)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE ERROR CHECKING ROUTINE FOR THE TIME SERIES
C     FOURIER UNIVARIATE SPECTRUM ANALYSIS ROUTINES.
C
C     WRITTEN BY - JANET R. DONALDSON
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  DECEMBER 2, 1985  (JRD)
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   IAR,LACOV,LAG,LAGMAX,LDSMIN,LDSTAK,LYFFT,N,NF,NFFT
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   ACOV(*),PHI(*)
      LOGICAL
     +   OPTION(4)
      CHARACTER
     +   NMSUB(6)*1
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      INTEGER
     +   I,NV
      LOGICAL
     +   HEAD
C
C  LOCAL ARRAYS
      LOGICAL
     +   ERR(20)
      CHARACTER
     +   L1(8)*1,LACV(8)*1,LACV1M(8)*1,LACV1P(8)*1,LIAR(8)*1,
     +   LLACOV(8)*1,LLAG(8)*1,LLDS(8)*1,LLGMX(8)*1,LLGMX1(8)*1,
     +   LLGMXM(8)*1,LLGMXP(8)*1,LLYFFT(8)*1,LN(8)*1,LNF(8)*1,
     +   LNM1(8)*1,LPHI(8)*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL EISGE,EISII,ERVII
C
C  INTRINSIC FUNCTIONS
      INTRINSIC ABS,IABS
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     DOUBLE PRECISION ACOV(LACOV)
C        THE AUTOCOVARIANCE FUNCTION.
C     LOGICAL ERR(20)
C        VALUE(S) INDICATING WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT
C        (FALSE).
C     LOGICAL HEAD
C        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
C        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
C        OF HEAD WILL BE CHANGED TO FALSE.
C     INTEGER IAR
C        THE ORDER OF THE AUTOREGRESSIVE PROCESS CHOSEN.
C     INTEGER IERR
C        THE INTEGER VALUE RETURNED BY THIS ROUTINE DESIGNATING
C        WHETHER ANY ERRORS WERE DETECTED IN THE PARAMETER LIST
C        IF ERR01, NO ERRORS WERE DETECTED
C        IF IERR .EQ. 1, ERRORS HAVE BEEN DETECTED
C     INTEGER LACOV
C        THE LENGTH OF THE VECTOR ACOV.
C     INTEGER LAGMAX
C        THE MAXIMUM LAG VALUE TO BE USED.
C     INTEGER LAG
C        THE LAG WINDOW TRUNCATION POINT USED FOR A SPECIFIC WINDOW.
C     INTEGER LDSMIN
C        THE MINIMUM LENGTH ALLOWED FOR DSTAK.
C     INTEGER LDSTAK
C        THE LENGTH OF THE VECTOR DSTAK IN COMMON CSTAK.
C     CHARACTER*1 LACV(8), LACV1M(8), LACV1P(8),
C    *   LIAR(8), LLACOV(8), LLAG(8), LLGMX(8), LLGMXM(8),
C    *   LLGMXP(8),  LLGMX1(8), LLDS(8), LN(8), LNF(8), LNM1(8),
C    *   LLYFFT(8), LPHI(8), L1(8)
C        THE ARRAY(S) CONTAINING THE NAME(S) OF THE ARGUMENT(S)
C        CHECKED FOR ERRORS.
C     INTEGER LYFFT
C        THE LENGTH OF THE VECTOR YFFT.
C     INTEGER N
C        THE INTEGER NUMBER OF OBSERVATIONS IN THE SERIES.
C     INTEGER NF
C        THE NUMBER OF FREQUENCIES AT WHICH THE SPECTRUM IS
C        TO BE COMPUTED.
C     INTEGER NFFT
C        THE NUMBER OF OBSERVATIONS IN THE EXTENDED SERIES.
C     CHARACTER*1 NMSUB(6)
C        THE ARRAY CONTAINING THE NAME OF THE USER CALLED SUBROUTINE.
C     INTEGER NV
C        THE NUMBER OF VIOLATIONS FOUND WHEN CHECKING VECTOR LAGS.
C     LOGICAL OPTION(4)
C        AN INDICATOR ARRAY USED TO DESIGNATE WHETHER ANY OF THE
C        FOUR POSSIBLE OPTIONS (F, M, V, OR S) HAVE BEEN USED (TRUE)
C        OR NOT (FALSE).
C     DOUBLE PRECISION PHI(IAR)
C        THE ARRAY OF AUTOREGRESSIVE COEFFICIENTS FOR THE
C        SELECTED ORDER.
C
C
C     SET UP NAME ARRAYS
C
      DATA LACV(1), LACV(2), LACV(3), LACV(4), LACV(5), LACV(6),
     +   LACV(7), LACV(8) /'A','C','O','V',' ',' ',' ',' '/
      DATA LACV1M(1), LACV1M(2), LACV1M(3), LACV1M(4), LACV1M(5),
     +   LACV1M(6), LACV1M(7), LACV1M(8) /'-','A','C','O','V','(','1',
     +   ')'/
      DATA LACV1P(1), LACV1P(2), LACV1P(3), LACV1P(4), LACV1P(5),
     +   LACV1P(6), LACV1P(7), LACV1P(8) /'+','A','C','O','V','(','1',
     +   ')'/
      DATA LIAR(1), LIAR(2), LIAR(3), LIAR(4), LIAR(5),
     +   LIAR(6), LIAR(7), LIAR(8) /'I','A','R',' ',' ',' ',' ',
     +   ' '/
      DATA LLACOV(1), LLACOV(2), LLACOV(3), LLACOV(4), LLACOV(5),
     +   LLACOV(6), LLACOV(7), LLACOV(8) /'L','A','C','O','V',' ',' ',
     +   ' '/
      DATA LLAG(1), LLAG(2), LLAG(3), LLAG(4), LLAG(5), LLAG(6),
     +   LLAG(7), LLAG(8) /'L','A','G',' ',' ',' ',' ',' '/
      DATA LLGMX(1), LLGMX(2), LLGMX(3), LLGMX(4), LLGMX(5),
     +   LLGMX(6), LLGMX(7), LLGMX(8) /'L','A','G','M','A','X',' ',
     +   ' '/
      DATA LLGMXM(1), LLGMXM(2), LLGMXM(3), LLGMXM(4), LLGMXM(5),
     +   LLGMXM(6), LLGMXM(7), LLGMXM(8) /'-','L','A','G','M','A','X',
     +   ' '/
      DATA LLGMXP(1), LLGMXP(2), LLGMXP(3), LLGMXP(4), LLGMXP(5),
     +   LLGMXP(6), LLGMXP(7), LLGMXP(8) /'+','L','A','G','M','A','X',
     +   ' '/
      DATA LLGMX1(1), LLGMX1(2), LLGMX1(3), LLGMX1(4), LLGMX1(5),
     +   LLGMX1(6), LLGMX1(7), LLGMX1(8) /'L','A','G','M','A','X','+',
     +   '1'/
      DATA LLDS(1), LLDS(2), LLDS(3), LLDS(4), LLDS(5), LLDS(6),
     +   LLDS(7), LLDS(8) /'L','D','S','T','A','K',' ',' '/
      DATA LN(1), LN(2), LN(3), LN(4), LN(5), LN(6), LN(7), LN(8) /'N',
     +   ' ',' ',' ',' ',' ',' ',' '/
      DATA LNF(1), LNF(2), LNF(3), LNF(4), LNF(5), LNF(6), LNF(7),
     +   LNF(8) /'N','F',' ',' ',' ',' ',' ',' '/
      DATA LNM1(1), LNM1(2), LNM1(3), LNM1(4), LNM1(5), LNM1(6),
     +   LNM1(7), LNM1(8) /'N','-','1',' ',' ',' ',' ',' '/
      DATA LLYFFT(1), LLYFFT(2), LLYFFT(3), LLYFFT(4), LLYFFT(5),
     +   LLYFFT(6), LLYFFT(7), LLYFFT(8) /'L','Y','F','F','T',' ',' ',
     +   ' '/
      DATA LPHI(1), LPHI(2), LPHI(3), LPHI(4), LPHI(5), LPHI(6),
     +   LPHI(7), LPHI(8) /'P','H','I',' ',' ',' ',' ',' '/
      DATA L1(1), L1(2), L1(3), L1(4), L1(5), L1(6), L1(7), L1(8) /'1',
     +   ' ',' ',' ',' ',' ',' ',' '/
C
C     SET UP FOR ERROR CHECKING
C
C
      IERR = 0
      HEAD = .TRUE.
C
      DO 10 I=1,20
         ERR(I) = .FALSE.
   10 CONTINUE
C
C     CALL ERROR CHECKING ROUTINES
C
      CALL EISGE(NMSUB, LN, N, 17, 1, HEAD, ERR(1), LN)
C
      IF ((.NOT.OPTION(3))) GO TO 15
C
      CALL ERVII(NMSUB, LACV, ACOV, LAGMAX+1, -ABS(ACOV(1)),
     +   ABS(ACOV(1)), 0, HEAD, 4, NV, ERR(15), LACV1M, LACV1P)
C
      CALL EISII(NMSUB, LLGMX, LAGMAX, 1, N-1, 1, HEAD, ERR(2),
     +   L1, LNM1)
C
      IF (OPTION(2)) THEN
         CALL EISGE(NMSUB, LLACOV, LACOV, LAGMAX+1, 8, HEAD, ERR(3),
     +   LLGMX1)
      ELSE
         CALL EISGE(NMSUB, LLACOV, LACOV, LAGMAX+1, 7, HEAD, ERR(3),
     +   LLGMX1)
      END IF
C
   15 IF (OPTION(1) .AND. (.NOT.ERR(1)))
     +   CALL EISGE(NMSUB, LLYFFT, LYFFT, NFFT, 9, HEAD, ERR(4),
     +   LLYFFT)
C
      IF (OPTION(1) .AND. (.NOT.OPTION(4)))
     +   CALL EISGE(NMSUB, LLDS, LDSTAK, LDSMIN, 9, HEAD, ERR(5), LLDS)
C
      IF (OPTION(4)) GO TO 30
C
      DO 20 I=1,15
         IF (ERR(I)) GO TO 50
   20 CONTINUE
C
      RETURN
C
   30 CONTINUE
C
      CALL EISII(NMSUB, LIAR, IAR, -IABS(LAGMAX), IABS(LAGMAX), 1, HEAD,
     +   ERR(6), LLGMXM, LLGMXP)
C
      CALL ERVII(NMSUB, LPHI, PHI, IAR, -1.0D0, 1.0D0, 0, HEAD, 1, NV,
     +   ERR(7), L1, L1)
C
      IF (.NOT.OPTION(3))
     +   CALL EISII(NMSUB, LLGMX, LAGMAX, 1, N-1, 1, HEAD, ERR(2),
     +   L1, LNM1)
C
      CALL EISII(NMSUB, LLAG, LAG, -IABS(LAGMAX), IABS(LAGMAX), 1, HEAD,
     +   ERR(8), LLGMXM, LLGMXP)
C
      CALL EISGE(NMSUB, LNF, NF, 1, 1, HEAD, ERR(9), LNF)
C
      IF (ERR(1) .OR. ERR(2) .OR. ERR(9)) GO TO 50
C
      CALL EISGE(NMSUB, LLDS, LDSTAK, LDSMIN, 9, HEAD, ERR(14), LLDS)
C
      DO 40 I=1,15
         IF (ERR(I)) GO TO 50
   40 CONTINUE
C
      RETURN
C
   50 CONTINUE
      IERR = 1
      RETURN
C
      END
