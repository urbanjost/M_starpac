*PGMS
      SUBROUTINE PGMS (YFFT, N, NFFT, LYFFT, IEXTND, NF, PER, LPER,
     +   FREQ, LFREQ, NPRT)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS IS THE USER CALLABLE ROUTINE FOR COMPUTING
C     THE (RAW) PERIODOGRAM OF A SERIES (LONG CALL).
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
     +   IEXTND,LFREQ,LPER,LYFFT,N,NF,NFFT,NPRT
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   FREQ(*),PER(*),YFFT(*)
C
C  SCALARS IN COMMON
      INTEGER
     +   IERR
C
C  LOCAL SCALARS
      INTEGER
     +   IPRT,NFFT2
      LOGICAL
     +   ERR01,ERR02,ERR03,ERR04,ERR05,HEAD
C
C  LOCAL ARRAYS
      CHARACTER
     +   LLFREQ(8)*1,LLPER(8)*1,LLYFFT(8)*1,LN(8)*1,NMSUB(6)*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL EISGE,ENFFT,IPRINT,PGMMN
C
C  COMMON BLOCKS
      COMMON /ERRCHK/IERR
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     LOGICAL ERR01, ERR02, ERR03, ERR04, ERR05
C        VALUE(S) INDICATING WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT
C        (FALSE).
C     DOUBLE PRECISION FREQ(LFREQ)
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
C     INTEGER IEXTND
C        THE INDICATOR VARIABLE USED TO DESIGNATE WHETHER ZERO
C        (IEXTND .EQ. 0) OR THE SERIES MEAN (IEXTND .NE. 0) IS TO BE
C        USED TO EXTEND THE SERIES.
C     INTEGER IPRT
C        THE LOGICAL UNIT USED FOR PRINTED OUTPUT.
C     INTEGER LFREQ
C        THE LENGTH OF THE ARRAY FREQ.
C     CHARACTER*1 LLFREQ(8), LLPER(8), LLYFFT(8), LN(8)
C        THE ARRAY(S) CONTAINING THE NAME(S) OF THE PARAMETER(S) CHECKED
C        FOR ERRORS.
C     INTEGER LPER
C        THE LENGTH OF THE ARRAY PER.
C     INTEGER LYFFT
C        THE LENGTH OF THE VECTOR YFFT.
C     INTEGER N
C        THE INTEGER NUMBER OF OBSERVATIONS.
C     INTEGER NF
C        THE NUMBER OF FREQUENCIES AT WHICH THE PERIODGRAM IS
C        TO BE COMPUTED.
C     INTEGER NFFT
C        THE EFFECTIVE LENGTH OF THE SERIES TO BE TRANSFORMED.
C     INTEGER NFFT2
C        THE EFFECTIVE SERIES LENGTH ACTUALLY USED.
C     INTEGER NPRT
C        THE VARIABLE CONTROLING PRINTED OUTPUT, WHERE
C        IF NPRT .LE. -2, THE OUTPUT CONSISTS OF A PAGE PLOT OF THE
C                         PERIODOGRAM ON A LOG-LINEAR SCALE,
C        IF NPRT .EQ. -1, THE OUTPUT CONSISTS OF A PAGE PLOT OF THE
C                         PERIODOGRAM IN DECIBELS ON A LINEAR SCALE,
C        IF NPRT .EQ.  0, THE OUTPUT IS SUPPRESSED,
C        IF NPRT .EQ.  1, THE OUTPUT CONSISTS OF A VERTICAL PLOT OF THE
C                         PERIODOGRAM IN DECIBELS ON A LINEAR SCALE.
C        IF NPRT .GE.  2, THE OUTPUT CONSISTS OF A VERTICAL PLOT OF THE
C                         PERIODOGRAM ON A LOG-LINEAR SCALE,
C     CHARACTER*1 NMSUB(6)
C        THE ARRAY CONTAINING THE NAME OF THIS SUBROUTINE.
C     DOUBLE PRECISION PER(LPER)
C        THE ARRAY IN WHICH THE PERIODOGRAM IS STORED.
C     DOUBLE PRECISION YFFT(LYFFT)
C        THE ARRAY CONTAINING THE OBSERVED TIME SERIES.
C
C     SET UP NAME ARRAYS
C
      DATA
     +  NMSUB(1),  NMSUB(2),  NMSUB(3),  NMSUB(4),  NMSUB(5),  NMSUB(6)
     + /     'P',       'G',       'M',       'S',       ' ',       ' '/
      DATA
     + LLFREQ(1), LLFREQ(2), LLFREQ(3), LLFREQ(4), LLFREQ(5),
     +  LLFREQ(6), LLFREQ(7), LLFREQ(8)
     +  /'L','F','R','E','Q',' ',' ',' '/
      DATA
     + LLPER(1), LLPER(2), LLPER(3), LLPER(4), LLPER(5),
     +  LLPER(6), LLPER(7), LLPER(8) /'L','P','E','R',' ',' ',' ',' '/
      DATA
     + LLYFFT(1), LLYFFT(2), LLYFFT(3), LLYFFT(4), LLYFFT(5),
     +  LLYFFT(6), LLYFFT(7), LLYFFT(8)
     +  /'L','Y','F','F','T',' ',' ',' '/
      DATA
     + LN(1), LN(2), LN(3), LN(4), LN(5), LN(6), LN(7), LN(8)
     + /'N',' ',' ',' ',' ',' ',' ',' '/
C
C     SET UP FOR ERROR CHECKING
C
      IERR = 0
      HEAD = .TRUE.
C
C     CALL ERROR CHECKING ROUTINES
C
      CALL EISGE(NMSUB, LN, N, 17, 1, HEAD, ERR01, LN)
      IF (ERR01) GO TO 5
C
      CALL ENFFT(NMSUB, NFFT, 2, N, LYFFT, NFFT2, HEAD, ERR02)
      NF = NFFT2/2
C
      CALL EISGE(NMSUB, LLYFFT, LYFFT, NFFT2, 9, HEAD, ERR03, LLYFFT)
C
      CALL EISGE(NMSUB, LLPER, LPER, NF, 9, HEAD, ERR04, LLPER)
C
      CALL EISGE(NMSUB, LLFREQ, LFREQ, NF, 9, HEAD, ERR05, LLFREQ)
C
      IF (ERR02 .OR. ERR03 .OR. ERR04 .OR. ERR05) GO TO 5
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
      CALL PGMMN (YFFT, N, NFFT2, IEXTND, NF, PER, LPER, YFFT, FREQ,
     +   LFREQ, NPRT, NMSUB)
C
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (/42H THE CORRECT FORM OF THE CALL STATEMENT IS//
     +  '       CALL PGMS (YFFT, N, NFFT, LYFFT,'/
     +  '      +           IEXTND, NF, PER, LPER, FREQ, LFREQ, NPRT)')
      END
