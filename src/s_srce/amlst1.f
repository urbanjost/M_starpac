*AMLST1
      SUBROUTINE AMLST1 (IAMHD, PAR, NPAR, MSPECT, NFAC, VCVL, LVCVL,
     +  SCALE, LSCALE, STPT, LSTPT, IPARMN, IPARMX, LBLTYP, T975, IFIXD)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE PRINTS THE PARAMETERS FOR THE ARIMA ROUTINES.
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
      REAL
     +   T975
      INTEGER
     +   IAMHD,IPARMN,IPARMX,LBLTYP,LSCALE,LSTPT,LVCVL,NFAC,NPAR
C
C  ARRAY ARGUMENTS
      REAL
     +   PAR(*),SCALE(*),STPT(*),VCVL(*)
      INTEGER
     +   IFIXD(*),MSPECT(NFAC,4)
C
C  LOCAL SCALARS
      REAL
     +   FPLM,PLL,PUL,RATIO,SDPAR
      INTEGER
     +   IPRT,J,K,L,LL,LPAR,ORDER
C
C  LOCAL ARRAYS
      CHARACTER
     +   FIXED(3)*1
C
C  EXTERNAL FUNCTIONS
      REAL
     +   R1MACH
      EXTERNAL R1MACH
C
C  EXTERNAL SUBROUTINES
      EXTERNAL FIXPRT,IPRINT
C
C  INTRINSIC FUNCTIONS
      INTRINSIC SQRT
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     CHARACTER*1 FIXED(3)
C        THE CHARACTERS USED TO LABEL THE PARAMETERS FIXED OR NOT.
C     REAL FPLM
C        THE FLOATING POINT LARGEST MAGNITUDE.
C     INTEGER IAMHD
C        THE INDICATOR VALUE USED TO DESIGNATE THE TYPE OF LIST
C        TO BE GENERATED
C        IF IAMHD=1, THE LIST IS FOR THE INITIAL SUMMARY OF THE
C                    ESTIMATION ROUTINES.
C        IF IAMHD=2, THE LIST IS FOR THE INITIAL REPORT OF THE
C                    FORECASTING ROUTINES.
C        IF IAMHD=3, THE LIST IS FOR THE FINAL REPORT OF THE
C                    ESTIMATION ROUTINES.
C     INTEGER IFIXD(NPAR)
C        THE INDICATOR VALUES USED TO DESIGNATE WHETHER THE
C        PARAMETERS ARE TO BE OPTIMIZED OR ARE TO BE HELD FIXED.
C        IF IFIXED(I).NE.0, THEN PAR(I) WILL BE HELD FIXED.
C        IF IFIXED(I).EQ.0, THEN PAR(I) WILL BE OPTIMIZED.
C     INTEGER IPARMN
C        THE SMALLEST PARAMETER INDEX INCLUDED IN THIS TERM.
C     INTEGER IPARMX
C        THE LARGEST PARAMETER INDEX INCLUDED IN THIS TERM.
C     INTEGER IPRT
C        THE UNIT NUMBER FOR PRINTED OUTPUT.
C     INTEGER LVCVL
C        THE DIMENSION OF VECTOR VCVL.
C     INTEGER J
C        AN INDEX VARIABLE.
C     INTEGER L
C        AN INDEX VARIABLE.
C     INTEGER LBLTYP
C        THE TYPE OF LABLE TO BE PRINTED, WHERE
C        1 INDICATES THE TERM IS AUTOREGRESSIVE AND
C        2 INDICATES THE TERM IS MOVING AVERAGE
C     INTEGER LL
C        AN INDEX VARIABLE.
C     INTEGER LPAR
C        AN INDEX VARIABLE.
C     INTEGER LSCALE
C        THE DIMENSION OF VECTOR SCALE.
C     INTEGER LSTPT
C        THE DIMENSION OF VECTOR STPT.
C     INTEGER MSPECT(NFAC,4)
C        THE ARRAY CONTAINING THE VALUES OF P, D, Q, AND S FOR EACH FACT
C     INTEGER NFAC
C        THE NUMBER OF FACTORS IN THE MODEL
C     INTEGER NPAR
C        THE NUMBER OF PARAMETERS IN THE MODEL.
C     INTEGER ORDER
C        THE ORDER OF B FOR THE PARAMETER BEING PRINTED
C     REAL PAR(NPAR)
C        THE CURRENT ESTIMATES OF THE PARAMETERS.
C     REAL PLL
C        THE LOWER CONFIDENCE LIMIT FOR A GIVEN PARAMETER.
C     REAL PUL
C        THE UPPER CONFIDENCE LIMIT FOR A GIVEN PARAMETER.
C     REAL RATIO
C        THE RATIO OF A GIVEN PARAMETER VALUE TO ITS STANDARD ERROR.
C     REAL SCALE(LSCALE)
C        THE TYPICAL SIZE OF THE PARAMETERS.
C     REAL SDPAR
C        THE STANDARD DEVIATION OF A GIVEN PARAMETER VALUE.
C     REAL STPT(LSTPT)
C        THE STEP SIZE ARRAY.
C     REAL T975
C        THE VALUE OF THE 97.5 PERCENT POINT FUNCTION FOR THE
C        T DISTRIBUTION.
C     REAL VCVL(LVCVL)
C        THE LOWER HALF OF THE VARIANCE-COVARIANCE MATRIX, STORED
C        ROW WISE.
C
C
      FPLM = R1MACH(2)
C
      CALL IPRINT(IPRT)
C
C     PRINT NEXT SET OF TERMS
C
      LPAR = 0
      DO 1 J=1,IPARMX
         IF (IFIXD(J).EQ.0) LPAR = LPAR + 1
    1 CONTINUE

      DO 40 J=1,NFAC
        IF ((MSPECT(J,LBLTYP).EQ.0) .AND. (LBLTYP.NE.2)) GO TO 40
        IF (LBLTYP.NE.2) IPARMX = IPARMX + MSPECT(J,LBLTYP)
        IF (LBLTYP.EQ.2) IPARMX = IPARMX + 1
        ORDER = 0
        DO 30 L = IPARMN, IPARMX
          ORDER = ORDER + MSPECT(J,4)
          IF (IAMHD.EQ.2) GO TO 25
          CALL FIXPRT(IFIXD(L), FIXED)
          IF (LBLTYP.EQ.1) WRITE(IPRT, 1000) L, J, ORDER,
     +      (FIXED(K),K=1,3), PAR(L)
          IF (LBLTYP.EQ.2) WRITE(IPRT, 1004) L,
     +      (FIXED(K),K=1,3), PAR(L)
          IF (LBLTYP.EQ.3) WRITE(IPRT, 1005) L, J, ORDER,
     +      (FIXED(K),K=1,3), PAR(L)
            IF (IAMHD.EQ.3) GO TO 10
C
            IF (IFIXD(L).EQ.0) GO TO 5
            WRITE (IPRT, 1007)
            GO TO 10
C
    5       CONTINUE
            IF (SCALE(1).LE.0.0E0) WRITE (IPRT, 1001) STPT(L)
            IF (SCALE(1).GT.0.0E0) WRITE (IPRT, 1002) SCALE(L), STPT(L)
   10     CONTINUE
          IF (IAMHD .EQ. 1) GO TO 30
C
          IF (IFIXD(L).EQ.0) GO TO 20
          WRITE(IPRT, 1006)
          GO TO 30
C
   20     CONTINUE
          LPAR = LPAR + 1
          RATIO = FPLM
          LL = LPAR*(LPAR-1)/2 + LPAR
          IF (VCVL(LL).GT.0.0E0) RATIO = PAR(L)/SQRT(VCVL(LL))
          SDPAR = SQRT(VCVL(LL))
          PLL = PAR(L) - T975*SDPAR
          PUL = PAR(L) + T975*SDPAR
          WRITE(IPRT, 1003) SDPAR, RATIO, PLL, PUL
          GO TO 30
   25     CONTINUE
          IF (LBLTYP.EQ.1) WRITE(IPRT, 1010) L, J, ORDER, PAR(L)
          IF (LBLTYP.EQ.2) WRITE(IPRT, 1014) L, PAR(L)
          IF (LBLTYP.EQ.3) WRITE(IPRT, 1015) L, J, ORDER, PAR(L)
   30   CONTINUE
        IPARMN = IPARMX + 1
   40 CONTINUE
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT(1X, I5, 2X, 'AR (FACTOR', I2, ')',4X,I5,6X,3A1,E17.8)
 1001 FORMAT ('+', 65X, 7HDEFAULT, E17.8)
 1002 FORMAT ('+', 55X, 2E17.8)
 1003 FORMAT ('+', 55X, 4(2X, E15.8))
 1004 FORMAT(1X, I5, 13X, 'MU', 4X, '  ---' ,6X,3A1,E17.8)
 1005 FORMAT(1X, I5, 2X, 'MA (FACTOR', I2, ')',4X,I5,6X,3A1,E17.8)
 1006 FORMAT('+', 55X, 4(14X, '---'))
 1007 FORMAT('+', 69X, '---', 14X, '---')
 1010 FORMAT(1X, I5, 2X, 'AR (FACTOR', I2, ')',4X,I5,E17.8)
 1014 FORMAT(1X, I5, 13X, 'MU', 4X, '  ---' ,E17.8)
 1015 FORMAT(1X, I5, 2X, 'MA (FACTOR', I2, ')',4X,I5,E17.8)
      END
