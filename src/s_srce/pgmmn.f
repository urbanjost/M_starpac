*PGMMN
      SUBROUTINE PGMMN (YFFT, N, NFFT, IEXTND, NF, PER, LPER, YAXIS,
     +   FREQ, LFREQ, NPRT, NMSUB)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS THE MAIN ROUTINE FOR COMPUTING THE RAW PERIODOGRAM.
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
     +   IEXTND,LFREQ,LPER,N,NF,NFFT,NPRT
C
C  ARRAY ARGUMENTS
      REAL
     +   FREQ(LFREQ),PER(LPER),YAXIS(LFREQ),YFFT(NFFT)
      CHARACTER
     +   NMSUB(6)*1
C
C  LOCAL SCALARS
      REAL
     +   YEXTND
      INTEGER
     +   I,N1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL AMEAN,PGMEST,PGORD,PGOUT,SETFRQ
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     REAL FREQ(LFREQ)
C        THE ARRAY IN WHICH THE FREQUENCIES CORRESPONDING TO THE
C        INTEGRATED SPECTRUM VALUES ARE STORED.
C     INTEGER I
C        AN INDEX VARIABLE.
C     INTEGER IEXTND
C        THE INDICATOR VARIABLE USED TO DESIGNATE WHETHER ZERO
C        (IEXTND .EQ. 0) OR THE SERIES MEAN (IEXTND .NE. 0) IS TO BE
C        USED TO EXTEND THE SERIES.
C     INTEGER LFREQ
C        THE LENGTH OF THE ARRAY FREQ.
C     INTEGER LPER
C        THE LENGTH OF THE ARRAY PER.
C     INTEGER N
C        THE INTEGER NUMBER OF OBSERVATIONS IN EACH SERIES
C     INTEGER NF
C        THE NUMBER OF FREQUENCIES AT WHICH THE PERIODGRAM IS
C        TO BE COMPUTED.
C     INTEGER NFFT
C        THE EFFECTIVE LENGTH OF THE SERIES TO BE TRANSFORMED.
C     CHARACTER*1 NMSUB(6)
C        THE NAME OF THE CALLING SUBROUTINE
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
C     REAL PER(LPER)
C        THE ARRAY IN WHICH THE PERIODOGRAM IS STORED.
C     REAL YAXIS(LFREQ)
C        THE ARRAY IN WHICH THE Y AXIS VALUES TO BE PLOTTED ARE STORED.
C     REAL YEXTND
C        THE VALUE USED TO EXTEND THE SERIES.
C     REAL YFFT(NFFT)
C        THE ARRAY CONTAINING THE OBSERVED TIME SERIES.
C
      YEXTND = 0.0E0
      IF (IEXTND .NE. 0) CALL AMEAN (YFFT, N, YEXTND)
C
C     EXTEND THE PERIODOGRAM ARRAY BY ITS MEAN  OR ZERO TO THE
C     EXTENDED LENGTH NFFT.
C
      N1 = N + 1
C
      DO 40 I = N1, NFFT
         YFFT(I) = YEXTND
   40 CONTINUE
C
C     COMPUTE THE PERIODOGRAM.
C
      CALL PGMEST (YFFT, NFFT, NF, 1.0E0, PER, LPER)
C
C     SET FREQUENCIES FOR PERIODOGRAM VALUES
C
      CALL SETFRQ (FREQ, NF, 1, 0.0E0, 0.5E0, 1.0E0)
C
      IF (NPRT .EQ. 0) RETURN
C
C     SET Y CO-ORDINATES FOR PERIODOGRAM PLOT.
C
      CALL PGORD (PER, NF, YAXIS, NPRT)
C
C     PLOT PERIODOGRAM IF OUTPUT NOT SUPPRESSED
C
      CALL PGOUT (YAXIS, FREQ, NF, NPRT, NMSUB)
C
      RETURN
C
      END
