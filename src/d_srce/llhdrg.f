*LLHDRG
      SUBROUTINE LLHDRG(PAGE, WIDE, ISUBHD)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS SUBROUTINE PRINTS THE PAGE HEADINGS FOR THE UNRESTRICTED
C     LINEAR LEAST SQUARES ROUTINES.
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
     +   ISUBHD
      LOGICAL
     +   PAGE,WIDE
C
C  LOCAL SCALARS
      INTEGER
     +   IPRT
C
C  EXTERNAL SUBROUTINES
      EXTERNAL IPRINT,VERSP
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     INTEGER IPRT
C        THE UNIT NUMBER FOR PRINTED OUTPUT.
C     INTEGER ISUBHD
C        AN INDICATOR VALUE SPECIFYING SUBHEADINGS TO BE PRINTED.
C     LOGICAL PAGE
C        THE VARIABLE USED TO INDICATE WHETHER A GIVEN SECTION OF
C        THE OUTPUT IS TO BEGIN ON A NEW PAGE (TRUE) OR NOT (FALSE).
C     LOGICAL WIDE
C        THE VARIABLE USED TO INDICATE WHETHER THE HEADING SHOULD
C        BE FULL WIDTH (TRUE) OR NOT (FALSE).
C
C
      CALL IPRINT(IPRT)
      IF (PAGE) WRITE (IPRT,1020)
      CALL VERSP(WIDE)
      IF (PAGE) WRITE (IPRT,1000)
      IF (.NOT.PAGE) WRITE (IPRT,1010)
      PAGE = .TRUE.
C
      IF (ISUBHD.EQ.0) RETURN
C
      GO TO (10), ISUBHD
C
   10 WRITE (IPRT,1030)
C
      RETURN
C
C     FORMAT STATEMENTS FOR PAGE HEADINGS
C
 1000 FORMAT (32H+LINEAR LEAST SQUARES ESTIMATION,
     +  ' WITH USER-SPECIFIED MODEL, CONTINUED')
 1010 FORMAT ('+', 63('*')/
     +   1X, 34H*  LINEAR LEAST SQUARES ESTIMATION,
     +   ' WITH USER-SPECIFIED MODEL  *'/ 1X, 63('*'))
 1020 FORMAT ('1')
 1030 FORMAT (//30H SUMMARY OF INITIAL CONDITIONS/ 1X, 30('-'))
      END
