*SUMOT
      SUBROUTINE SUMOT(STS, N, NNZW, WTS)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS SUBROUTINE REPORTS THE RESULTS OF A STAT FAMILY
C     COMPUTATION OF 53 SELECTED STATISTICS.  THERE MAY OR
C     MAY NOT BE WEIGHTS.
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   N,NNZW
      LOGICAL
     +   WTS
C
C  ARRAY ARGUMENTS
      REAL
     +   STS(53)
C
C  LOCAL SCALARS
      INTEGER
     +   I,IPRT
C
C  LOCAL ARRAYS
      INTEGER
     +   ITEMP(10)
C
C  EXTERNAL SUBROUTINES
      EXTERNAL IPRINT,VERSP
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     INTEGER I
C        LOOP INDEX.
C     INTEGER IPRT
C        THE OUTPUT UNIT.
C     INTEGER ITEMP(10)
C        A TEMPORARY VECTOR OF INTEGER STORAGE.
C     INTEGER N
C        INPUT PARAMETER.  THE LENGTH OF THE ORIGINAL DATA VECTOR.
C     INTEGER NNZW
C        INPUT PARAMETER.  THE NUMBER OF VALUES IN THE ORIGINAL DATA
C        VECTOR WITH WEIGHTS GREATER THAN 0.0E0.
C     REAL STS(53)
C        INPUT PARAMETER.  THE VECTOR OF 53 STATISTICS COMPUTED.
C        ROW STATISTIC                    ROW STATISTIC
C         1  LENGTH OF VECTOR             TESTS FOR NONRANDOMNESS
C         2  NUMBER OF NONZERO WEIGHTS    23  NUMBER OF RUNS UP AND DOWN
C        MEASURES OF LOCATION             24  EXPECTED NUMBER OF RUNS
C         3  UNWEIGHTED MEAN              25  S.D. OF NUMBER OF RUNS
C         4  WEIGHTED MEAN                26  MEAN SQR. SUCCESSIVE DIFF.
C         5  MEDIAN                       27  MEAN SQR. SUCC. DIFF./VAR.
C         6  MID-RANGE                    DEVIATIONS FROM WTD MEAN
C         7  25 P.C. UNWTD. TRIMMED MEAN  28  NUMBER OF + SIGNS
C         8  25 P.C. WTD. TRIMMED MEAN    29  NUMBER OF - SIGNS
C        MEASURES OF DISPERSION           30  NUMBER OF RUNS
C         9  STANDARD DEVIATION (S.D.)    31  EXPECTED NUMBER OF RUNS
C        10  S.D. OF MEAN                 32  S.D. OF RUNS
C        11  RANGE                        33  DIFF./S.D. OF RUNS
C        12  MEAN VARIATION               OTHER STATISTICS
C        13  VARIANCE (VAR.)              34  MINIMUM
C        14  COEFFICIENT OF VARIATION     35  MAXIMUM
C        CONFIDENCE INTERVALS             36  BETA 1
C        15  LOWER CONFIDENCE LIMIT, MEAN 37  BETA 2
C        16  UPPER CONFIDENCE LIMIT, MEAN 38  WTD. SUM OF VALUES
C        17  LOWER CONFIDENCE LIMIT, S.D. 39  WTD. SUM OF SQUARES
C        18  UPPER CONFIDENCE LIMIT, S.D. 40  WTD. SUM OF SQRD. DEVS.
C        LINEAR TREND STATISTICS          41  STUDENTS T
C        19  SLOPE                        42  WTD. SUM OF ABS. VALUES
C        20  S.D. OF SLOPE                43  WTD. AVG. ABS. VALUES
C        21  SLOPE/S.D. OF SLOPE = T      44-53 FREQ. DISTRIBUTION
C        22  PROB ( X .GT. ABS(OBS. T))
C     LOGICAL WTS
C        INPUT PARAMETER.  A FLAG TO INDICATE WHETHER OR NOT THERE ARE
C        WEIGHTS.
C
C     BEGIN PRINTOUT
C
      CALL IPRINT(IPRT)
C
C     PRINT HEADING
C
      CALL VERSP(.TRUE.)
C
C     PRINT NUMBERS OF OBSERVATIONS, RAW AND NONZERO WEIGHTED.
C
      IF (.NOT.WTS) WRITE (IPRT,1000)
      IF (WTS) WRITE (IPRT,1010)
      IF (NNZW.NE.N) GO TO 10
      WRITE (IPRT,1020) NNZW
      GO TO 20
   10 WRITE (IPRT,1030) NNZW, N
      WRITE (IPRT,1040)
C
C     PRINT FREQUENCY DISTRIBUTIONS
C
   20 DO 30 I=1,10
         ITEMP(I) = STS(I+43)
   30 CONTINUE
      WRITE (IPRT,1050) (ITEMP(I),I=1,10)
C
C     PRINT MEASURES OF LOCATION AND DISPERSION
C
      WRITE (IPRT,1060)
      IF (STS(4).NE.0.0E0)
     +   WRITE (IPRT,1070) (STS(I+2),STS(I+8),I=1,6)
      IF (STS(4).EQ.0.0E0)
     +   WRITE (IPRT,1080) (STS(I+2),STS(I+8),I=1,5), STS(8)
C
C     PRINT CONFIDENCE INTERVALS
C
      WRITE (IPRT,1090) (STS(I),I=15,18)
C
C     PRINT LINEAR TREND AND OTHER STATISTICS, AND PRINT HEADING FOR
C     TESTS FOR NONRANDOMNESS
C
      WRITE (IPRT,1100)
     +   (STS(I),STS(I+15),I=19,22), (STS(I),I=38,41)
      ITEMP(1) = STS(23)
      ITEMP(2) = STS(28)
      ITEMP(3) = STS(29)
      ITEMP(4) = STS(30)
C
C     PRINT TESTS FOR NONRANDOMNESS
C
      WRITE (IPRT,1110) ITEMP(1), STS(42), STS(24), STS(43),
     +   (STS(I),I=25,27), (ITEMP(I),I=2,4), (STS(I),I=31,33)
C
C     PRINT FOOTNOTE
C
      WRITE (IPRT,1120)
      RETURN
C
 1000 FORMAT('+STATISTICAL ANALYSIS')
 1010 FORMAT('+WEIGHTED STATISTICAL ANALYSIS')
 1020 FORMAT(//5X, 4HN = , I5)
 1030 FORMAT(//5X, 4HN = , I5, 32H (NO. OF NON-ZERO WTS)    LENGTH,
     +   2H =, I5)
 1040 FORMAT(/5X, 45HALL COMPUTATIONS ARE BASED ON OBSERVATIONS WI,
     +   19HTH NON-ZERO WEIGHTS)
 1050 FORMAT(//5X, 28HFREQUENCY DISTRIBUTION (1-6), 7X, 10I6)
 1060 FORMAT(//5X, 26HMEASURES OF LOCATION (2-2), 34X, 10HMEASURES O,
     +   18HF DISPERSION (2-6))
 1070 FORMAT(/10X, 26HUNWEIGHTED MEAN          =, 1PE15.7, 20X,
     +   26HWTD STANDARD DEVIATION   =, E15.7/10X, 17HWEIGHTED MEAN    ,
     +   9H        =, E15.7, 20X, 26HWEIGHTED S.D. OF MEAN    =,
     +   E15.7/10X, 26HMEDIAN                   =, E15.7, 20X, 6HRANGE ,
     +   20H                   =, E15.7/10X, 23HMID-RANGE              ,
     +   3H  =, E15.7, 20X, 26HMEAN DEVIATION           =, E15.7/10X,
     +   26H25 PCT UNWTD TRIMMED MEAN=, E15.7, 20X, 16HVARIANCE        ,
     +   10H         =, E15.7/10X, 26H25 PCT WTD TRIMMED MEAN  =,
     +   E15.7, 20X, 26HCOEF. OF. VAR. (PERCENT) =, E15.7)
 1080 FORMAT(/10X, 26HUNWEIGHTED MEAN          =, 1PE15.7, 20X,
     +   26HWTD STANDARD DEVIATION   =, E15.7/10X, 17HWEIGHTED MEAN    ,
     +   9H        =, E15.7, 20X, 26HWEIGHTED S.D. OF MEAN    =,
     +   E15.7/10X, 26HMEDIAN                   =, E15.7, 20X, 6HRANGE ,
     +   20H                   =, E15.7/10X, 23HMID-RANGE              ,
     +   3H  =, E15.7, 20X, 26HMEAN DEVIATION           =, E15.7/10X,
     +   26H25 PCT UNWTD TRIMMED MEAN=, E15.7, 20X, 16HVARIANCE        ,
     +   10H         =, E15.7/10X, 26H25 PCT WTD TRIMMED MEAN  =,
     +   E15.7, 20X, 26HCOEFFICIENT OF VARIATION =, 13H    UNDEFINED/
     +   98X, 14H(MEAN IS ZERO))
 1090 FORMAT(///20X, 46HA TWO-SIDED 95 PCT CONFIDENCE INTERVAL FOR MEA,
     +   4HN IS, 1PE14.7, 4H TO , E14.7, 6H (2-2)/20X, 13HA TWO-SIDED 9,
     +   37H5 PCT CONFIDENCE INTERVAL FOR S.D. IS, E14.7, 4H TO ,
     +   E14.7, 6H (2-7))
 1100 FORMAT(///5X, 30HLINEAR TREND STATISTICS (5-1) , 30X, 6HOTHER ,
     +   10HSTATISTICS//10X, 5HSLOPE, 20X, 1H=, 1PE15.7, 20X, 7HMINIMUM,
     +   18X, 1H=, E15.7/10X, 13HS.D. OF SLOPE, 12X, 1H=, E15.7, 20X,
     +   7HMAXIMUM, 18X, 1H=, E15.7/10X, 26HSLOPE/S.D. OF SLOPE = T  =,
     +   E15.7, 20X, 8HBETA ONE, 17X, 1H=, E15.7/10X, 14HPROB EXCEEDING,
     +   21H ABS VALUE OF OBS T =, 0PF6.3, 20X, 8HBETA TWO, 17X, 1H=,
     + 1PE15.7/71X, 17HWTD SUM OF VALUES, 8X, 1H=, E15.7/71X, 7HWTD SUM,
     +   11H OF SQUARES, 7X, 1H=, E15.7/5X, 24HTESTS FOR NON-RANDOMNESS,
     +   42X, 22HWTD SUM OF DEV SQUARED, 4H   =, E15.7/71X, 9HSTUDENTS ,
     +   'T', 15X, 1H=, E15.7)
 1110 FORMAT(10X, 26HNO. OF RUNS UP AND DOWN  =, I5, 30X, 9HWTD SUM A,
     +   17HBSOLUTE VALUES  =, 1PE15.7/
     +   10X, 26HEXPECTED NO. OF RUNS     =,
     +   0PF7.1, 28X, 26HWTD AVE ABSOLUTE VALUES  =, 1PE15.7/
     +   10X, 26HS.D. OF NO. OF RUNS      =, 0PF8.2/
     +   10X, 26HMEAN SQ SUCCESSIVE DIFF  =, 1X, 1PE16.7/
     +   10X, 26HMEAN SQ SUCC DIFF/VAR    =, 0PF9.3///
     +   10X, 24HDEVIATIONS FROM WTD MEAN//
     +   15X, 21HNO. OF + SIGNS      =, I5/
     +   15X, 21HNO. OF - SIGNS      =, I5/
     +   15X, 21HNO. OF RUNS         =, I5/
     +   15X, 21HEXPECTED NO. OF RUNS=, F7.1/
     +   15X, 12HS.D. OF RUNS, 8X, 1H=, F8.2/
     +   15X, 21HDIFF./S.D. OF RUNS  =, F9.3)
 1120 FORMAT(///49H NOTE - ITEMS IN PARENTHESES REFER TO PAGE NUMBER,
     +   36H IN NBS HANDBOOK 91 (NATRELLA, 1966))
      END
