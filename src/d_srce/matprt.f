*MATPRT
      SUBROUTINE MATPRT (X, Y, NC, IPRT, MODE, CODE, IRDIM)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS SUBROUTINE TAKES A SQUARE MATRIX AND PRINTS EITHER ITS
C     LOWER TRIANGULAR PART OR THE FULL MATRIX WITH OR WITHOUT DOUBLE
C     PRINTING.
C
C     WRITTEN BY - LINDA L. MITCHELL
C                  STATISTICAL ENGINEERING LAB/BOULDER
C                  NATIONAL BUREAU OF STANDARDS
C
C     CREATION DATE  -  APRIL 2, 1981
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   CODE,IPRT,IRDIM,MODE,NC
C
C  ARRAY ARGUMENTS
      DOUBLE PRECISION
     +   X(IRDIM,NC),Y(IRDIM,NC)
C
C  LOCAL SCALARS
      DOUBLE PRECISION
     +   TEMP
      INTEGER
     +   I,IWIDTH,J,K,KM,KN,L,NF
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MIN,SQRT
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     INTEGER CODE
C                  IF 1 -SINGLE PRINTED LINE, X ONLY (Y IS DUMMY ARG)
C                     2 -DOUBLE PRINTED LINE, BOTH X AND Y
C     INTEGER I
C                  ROW NUMBER
C     INTEGER IPRT
C                  THE OUTPUT UNIT NUMBER
C     INTEGER IRDIM
C                  ROW INDEX OF X
C     INTEGER IWIDTH
C                  THE WIDTH OF THE OUTPUT DEVICE.
C     INTEGER J
C                  FIRST COLUMN IN THE SET TO BE PRINTED
C     INTEGER K
C                  COLUMN NUMBER IN THE POSSIBLE SET OF NF
C     INTEGER KM
C                  LAST COLUMN IN THE SET
C                   LIMITED TO VALUES OF J-1 PLUS A NUMBER BETWEEN 1 AND
C                   NF (INCLUSIVE)
C     INTEGER KN
C                  LAST COLUMN TO PRINT WHEN PRINTING LOWER TRIANGLE
C     INTEGER L
C                  FIRST ROW TO PRINT FOR THIS SET
C     INTEGER MODE
C                  IF 0, LOWER TRIANGULAR PART PRINTED
C                     1, FULL MATRIX PRINTED
C                     2, LOWER TRIANGULAR PART IS PRINTED WITH
C                        SQUARE ROOTS OF THE DIAGONAL
C     INTEGER NC
C                  ROW AND COLUMN DIMENSION OF X
C     INTEGER NF
C                  THE NUMBER OF COLUMNS THAT CAN BE PRINTED, GIVEN
C                  THE WIDTH IWIDTH OF THE OUTPUT DEVICE.
C     DOUBLE PRECISION TEMP
C                  A TEMPORARY LOCATION
C     DOUBLE PRECISION X(IRDIM,NC)
C                  NC BY NC INPUT MATRIX
C     DOUBLE PRECISION Y(IRDIM,NC)
C                  MATRIX TO BE PRINTED ON THE SECOND LEVEL IF CODE=2
C
      IWIDTH = 132
      NF = MIN(7, (IWIDTH - 7)/17)
      L = 1
      DO 20 J=1,NC, NF
         KN = MIN(NC, J+NF-1)
         WRITE(IPRT,1000) (K,K=J,KN)
         WRITE(IPRT,1030)
         IF ((MODE.EQ.00) .OR. (MODE.EQ.2)) L = J
         DO 10 I=L,NC
            TEMP = X(I,I)
            KM = KN
            IF ((MODE.EQ.0) .OR. (MODE.EQ.2))
     +         KM = J + MIN(I-L, NF-1)
            IF ((MODE.EQ.2) .AND. ((I.GE.J) .AND. (I.LE.KM)))
     +         X(I,I) = SQRT(X(I,I))
            WRITE(IPRT,1010) I, (X(I,K),K=J,KM)
            IF (CODE.EQ.2) WRITE(IPRT,1020) (Y(I,K),K=J,KM)
            IF (CODE.EQ.2) WRITE(IPRT,1030)
            X(I,I) = TEMP
   10    CONTINUE
   20 CONTINUE
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (/4X, 7HCOLUMN , 7(I9, 8X))
 1010 FORMAT (4X, I6, 1X, 7(3X, G14.8))
 1020 FORMAT (9X, 7(3X, G14.8))
 1030 FORMAT (4X)
      END
