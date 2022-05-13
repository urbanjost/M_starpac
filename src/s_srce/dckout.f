*DCKOUT
      SUBROUTINE DCKOUT(XM, IXM, N, M, NROW, NETA, NTAU, NPAR, MSG,
     +   LMSG, PAR, SCALE, LSCALE, HDR, PAGE, WIDE, ISUBHD, PRTFXD,
     +   IFIXD)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE PRINTS THE RESULTS OF THE DERIVATIVE CHECKING
C     SUBROUTINE
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
     +   ISUBHD,IXM,LMSG,LSCALE,M,N,NETA,NPAR,NROW,NTAU
      LOGICAL
     +   PAGE,PRTFXD,WIDE
C
C  ARRAY ARGUMENTS
      REAL
     +   PAR(NPAR),SCALE(LSCALE),XM(IXM,M)
      INTEGER
     +   IFIXD(NPAR),MSG(LMSG)
C
C  SUBROUTINE ARGUMENTS
      EXTERNAL HDR
C
C  LOCAL SCALARS
      INTEGER
     +   I,IMAX,IMIN,INDEX,IPRT,J,K,NPERL
      CHARACTER
     +   BLANK*1
C
C  LOCAL ARRAYS
      LOGICAL
     +   FTNOTE(6)
      CHARACTER
     +   FIXED(3)*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL FIXPRT,IPRINT
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MIN
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     CHARACTER*1 BLANK
C        THE CHARACTER BLANK.
C     CHARACTER*1 FIXED(3)
C        THE CHARACTERS USED TO LABEL THE PARAMETERS FIXED OR NOT.
C     LOGICAL FTNOTE(6)
C        THE ARRAY WHICH CONTROLS PRINTING OF FOOTNOTES.
C     EXTERNAL HDR
C        THE NAME OF THE ROUTINE WHICH PRODUCES THE HEADING
C     INTEGER I
C        AN INDEX VARIABLE
C     INTEGER IFIXD(NPAR)
C        THE INDICATOR VALUES USED TO DESIGNATE WHETHER THE
C        PARAMETERS ARE TO BE OPTIMIZED OR ARE TO BE HELD FIXED.  IF
C        IFIXD(I).NE.0, THEN PAR(I) WILL BE OPTIMIZED.  IF
C        IFIXD(I).EQ.0, THEN PAR(I) WILL BE HELD FIXED.
C     INTEGER IMAX, IMIN
C        THE LARGEST AND SMALLEST INDEX VALUE TO BE PRINTED ON EACH
C        LINE.
C     INTEGER INDEX
C        THE INDEX VALUE TO BE PRINTED.
C     INTEGER IPRT
C        THE UNIT NUMBER FOR PRINTED OUTPUT.
C     INTEGER ISUBHD
C        AN INDICATOR VALUE SPECIFYING SUBHEADINGS TO BE PRINTED BY
C        ROUTINE HDR.
C     INTEGER IXM
C        THE FIRST DIMENSION OF THE INDEPENDENT VARIABLE ARRAY.
C     INTEGER J
C        AN INDEX VARIABLE.
C     INTEGER K
C        AN INDEX VARIABLE.
C     INTEGER LMSG
C        THE LENGTH OF THE VECTOR MSG.
C     INTEGER LSCALE
C        THE LENGTH OF VECTOR SCALE.
C     INTEGER M
C        THE NUMBER OF INDEPENDENT VARIABLES.
C     INTEGER MSG(LMSG)
C        AN ARRAY USED TO STORE MESSAGE PARAMETERS.
C     INTEGER NETA
C        THE NUMBER OF RELIABLE DIGITS IN THE MODEL.
C     INTEGER NPAR
C        THE NUMBER OF UNKNOWN PARAMETERS IN THE MODEL.
C     INTEGER NPERL
C        THE NUMBER OF VALUES TO BE PRINTED PER LINE.
C     INTEGER NROW
C        THE NUMBER OF THE ROW OF THE INDEPENDENT VARIABLE ARRAY AT
C        WHICH THE DERIVATIVE IS TO BE CHECKED.
C     INTEGER NTAU
C        THE NUMBER OF DIGITS OF AGREEMENT REQUIRED BETWEEN THE
C        APPROXIMATED DERIVATIVES AND THE USER-SUPPLIED DERIVATIVES.
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
C       THE VARIABLE USED TO INDICATE WHETHER THE HEADING SHOULD
C        FULL WIDTH (TRUE) OR NOT (FALSE).
C     REAL XM(IXM,M)
C        THE INDEPENDENT VARIABLE.
C
      DATA BLANK /' '/
C
      CALL IPRINT(IPRT)
C
C     INITIALIZE ARRAY FIXED
C
      DO 10 K=1,3
         FIXED(K) = BLANK
   10 CONTINUE
C
      CALL HDR(PAGE, WIDE, ISUBHD)
C
C     SET UP FOR FOOTNOTES
C
      DO 20 I=1,6
         FTNOTE(I) = .FALSE.
   20 CONTINUE
C
      IF (MSG(1).LE.0) GO TO 40
C
      DO 30 I=1,NPAR
         IF ((MSG(I+1).EQ.0) .OR. (MSG(I+1).EQ.2)) GO TO 30
         K = MSG(I+1) - 2
         IF (K.EQ.-1) K = 5
         FTNOTE(1) = .TRUE.
         FTNOTE(K+1) = .TRUE.
   30 CONTINUE
C
C     PRINT REPORT
C
   40 CONTINUE
C
      WRITE (IPRT,1000)
      IF (FTNOTE(1)) WRITE (IPRT,1040)
      IF (PRTFXD) WRITE (IPRT,1160)
      IF (.NOT.PRTFXD) WRITE (IPRT,1170)
C
      IF (SCALE(1).LE.0.0E0) GO TO 60
C
      DO 50 I=1,NPAR
         IF (PRTFXD) CALL FIXPRT(IFIXD(I), FIXED)
         K = MSG(I+1) - 2
         IF (K.EQ.-1) K = 5
         IF (K.EQ.-2) WRITE (IPRT,1010) I, (FIXED(J),J=1,3), PAR(I),
     +      SCALE(I)
         IF (K.EQ.0) WRITE (IPRT,1020) I, (FIXED(J),J=1,3), PAR(I),
     +      SCALE(I)
         IF (K.GE.1) WRITE (IPRT,1030) I, (FIXED(J),J=1,3), PAR(I),
     +      SCALE(I), K
   50 CONTINUE
      GO TO 80
C
   60 CONTINUE
C
      DO 70 I=1,NPAR
         IF (PRTFXD) CALL FIXPRT(IFIXD(I), FIXED)
         K = MSG(I+1) - 2
         IF (K.EQ.-1) K = 5
         IF (K.EQ.-2) WRITE (IPRT,1180) I, (FIXED(J),J=1,3), PAR(I)
         IF (K.EQ.0) WRITE (IPRT,1190) I, (FIXED(J),J=1,3), PAR(I)
         IF (K.GE.1) WRITE (IPRT,1200) I, (FIXED(J),J=1,3), PAR(I), K
   70 CONTINUE
C
   80 CONTINUE
C
C     PRINT FOOTNOTES
C
      IF (.NOT.FTNOTE(1)) GO TO 90
C
      WRITE (IPRT,1060)
      IF (FTNOTE(2)) WRITE (IPRT,1070)
      IF (FTNOTE(3)) WRITE (IPRT,1080)
      IF (FTNOTE(4)) WRITE (IPRT,1090)
      IF (FTNOTE(5)) WRITE (IPRT,1100)
      IF (FTNOTE(6)) WRITE (IPRT,1050)
C
   90 CONTINUE
C
      WRITE (IPRT,1110) NETA
      WRITE (IPRT,1120) NTAU
C
C     PRINT OUT ROW OF INDEPENDENT VARIABLE WHICH WAS CHECKED.
C
      WRITE (IPRT,1130) NROW
      NPERL = 7
C
      DO 100 I=1,M,NPERL
         IMIN = I
         IMAX = MIN(I+NPERL-1,M)
         WRITE (IPRT,1140) (INDEX,INDEX=IMIN,IMAX)
         WRITE (IPRT,1150) (XM(NROW,INDEX),INDEX=IMIN,IMAX)
  100 CONTINUE
      WRITE (IPRT,1210) N
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (//)
 1010 FORMAT (1X, I3, 5X, 3A1, 2G17.8, 10X, 2HOK)
 1020 FORMAT (1X, I3, 5X, 3A1, 2G17.8, 7X, 9HINCORRECT)
 1030 FORMAT (1X, I3, 5X, 3A1, 2G17.8, 5X, 14HQUESTIONABLE (, I1, ')')
 1040 FORMAT (62X, 1H*)
 1050 FORMAT (/48H  (5) USER-SUPPLIED AND APPROXIMATED DERIVATIVES,
     +   14H DISAGREE, BUT/5X, 37H APPROXIMATED DERIVATIVE IS QUESTIONA,
     +   11HBLE BECAUSE, 6H RATIO/5X, 30H OF RELATIVE CURVATURE TO RELA,
     +   17HTIVE SLOPE IS TOO, 6H HIGH.)
 1060 FORMAT (/53H * NUMBERS IN PARENTHESES REFER TO THE FOLLOWING NOTE,
     +   2HS.)
 1070 FORMAT (/48H  (1) USER-SUPPLIED AND APPROXIMATED DERIVATIVES,
     +   11H AGREE, BUT/5X, 40H BOTH ARE ZERO.  RECHECK AT ANOTHER ROW.)
 1080 FORMAT (/48H  (2) USER-SUPPLIED AND APPROXIMATED DERIVATIVES,
     +   15H MAY AGREE, BUT/5X, 36H USER-SUPPLIED DERIVATIVE IS IDENTIC,
     +   9HALLY ZERO, 17H AND APPROXIMATED/5X, 21H DERIVATIVE IS ONLY A,
     +   18HPPROXIMATELY ZERO., 25H  RECHECK AT ANOTHER ROW.)
 1090 FORMAT (/48H  (3) USER-SUPPLIED AND APPROXIMATED DERIVATIVES,
     +   14H DISAGREE, BUT/5X, 37H USER-SUPPLIED DERIVATIVE IS IDENTICA,
     +   9HLLY ZERO., 12H  RECHECK AT/5X, 13H ANOTHER ROW.)
 1100 FORMAT (/48H  (4) USER-SUPPLIED AND APPROXIMATED DERIVATIVES,
     +   14H DISAGREE, BUT/5X, 37H APPROXIMATED DERIVATIVE IS QUESTIONA,
     +   11HBLE BECAUSE, 13H EITHER RATIO/5X, 22H OF RELATIVE CURVATURE,
     +   25H TO RELATIVE SLOPE IS TOO, 9H HIGH, OR/5X, 13H SCALE(K) IS ,
     +   6HWRONG.)
 1110 FORMAT (/43H NUMBER OF RELIABLE DIGITS IN MODEL RESULTS, 25X,
     +   6H(NETA), 1X, I5)
 1120 FORMAT (/40H NUMBER OF DIGITS IN DERIVATIVE CHECKING, 9H AGREEMEN,
     +   11HT TOLERANCE, 8X, 6H(NTAU), 1X, I5)
 1130 FORMAT (/45H ROW NUMBER AT WHICH DERIVATIVES WERE CHECKED, 23X,
     +   6H(NROW), 1X, I5/42H   -VALUES OF THE INDEPENDENT VARIABLES AT,
     +   9H THIS ROW)
 1140 FORMAT (10X, 5HINDEX, I5, 6I15)
 1150 FORMAT (10X, 5HVALUE, 7(1X, G14.7)/)
 1160 FORMAT (52X, 10HDERIVATIVE/7X, 24HPARAMETER STARTING VALUE, 6X,
     +   5HSCALE, 10X, 10HASSESSMENT/1X, 5HINDEX, 2X, 5HFIXED, 6X,
     +   5H(PAR), 12X, 7H(SCALE)/)
 1170 FORMAT (17X, 9HPARAMETER, 26X, 10HDERIVATIVE/15X, 12HSTARTING VAL,
     +   2HUE, 8X, 5HSCALE, 10X, 10HASSESSMENT/1X, 5HINDEX, 13X,
     +   5H(PAR), 12X, 7H(SCALE)/)
 1180 FORMAT (1X, I3, 5X, 3A1, G17.8, 7X, 7HDEFAULT, 13X, 2HOK)
 1190 FORMAT (1X, I3, 5X, 3A1, G17.8, 7X, 7HDEFAULT, 10X, 9HINCORRECT)
 1200 FORMAT (1X, I3, 5X, 3A1, G17.8, 7X, 7HDEFAULT, 8X, 11HQUESTIONABL,
     +   3HE (, I1, ')')
 1210 FORMAT (/23H NUMBER OF OBSERVATIONS, 48X, 3H(N), 1X, I5)
      END
