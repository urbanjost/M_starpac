*EISGE
      SUBROUTINE EISGE(NMSUB, NMVAR1, NVAL, NMIN, MSGTYP, HEAD, ERROR,
     +   NMVAR2)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE CHECKS WHETHER THE VALUE   NVAL   IS GREATER THAN
C     OR EQUAL TO   NMIN   AND PRINTS A DIAGNOSTIC IF IT IS NOT.
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
     +   MSGTYP,NMIN,NVAL
      LOGICAL
     +   ERROR,HEAD
C
C  ARRAY ARGUMENTS
      CHARACTER
     +   NMSUB(6)*1,NMVAR1(8)*1,NMVAR2(8)*1
C
C  LOCAL SCALARS
      INTEGER
     +   I,IPRT
C
C  EXTERNAL SUBROUTINES
      EXTERNAL EHDR,IPRINT
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     LOGICAL ERROR
C        VALUE(S) INDICATING WHETHER AN ERROR WAS DETECTED (TRUE) OR NOT
C        (FALSE).
C     LOGICAL HEAD
C        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
C        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
C        OF HEAD WILL BE CHANGED TO FALSE.
C     INTEGER I
C        AN INDEX ARGUMENT.
C     INTEGER IPRT
C        THE UNIT NUMBER FOR PRINTED OUTPUT.
C     INTEGER MSGTYP
C        AN ARGUMENT USED TO INDICATE THE TYPE OF MESSAGE TO BE
C        PRINTED, WHERE IF ERROR IS TRUE AND
C        MSGTYP = 1 THE INPUT VALUE WAS TOO SMALL BASED
C                   ON LIMITS IMPOSED BY STARPAC
C        MSGTYP = 2 THE INPUT VALUE WAS TOO SMALL BASED ON OTHER INPUT
C                   ARGUMENTS.
C        MSGTYP = 3 THE INPUT VALUE WAS TOO SMALL BASED ON OTHER INPUT
C                   ARGUMENTS, WHERE THE VALUE INDICATES THE FIRST
C                   DIMENSION OF A DIMENSIONED ARRAY
C                   N.B.  IT IS ASSUMED THAT THE DIMENSION NAME IS THE
C                         ARRAY NAME PRECEDED BY THE LETTER I.  IF THE
C                         ARRAY NAME IS 6 LETTERS, THE DIMENSION NAME
C                         SHOULD OMIT THE LAST LETTER.  THE DIMENSION
C                         NAME WILL BE PRINTED USING (NMVAR(I),I=1,6),
C                         AND THE ARRAY NAME USING (NMVAR(I),I=2,7).
C        MSGTYP = 4 THE INPUT VALUE WAS TOO SMALL BASED ON OTHER INPUT
C                   ARGUMENTS, WHERE THE VALUE INDICATES THE SECOND
C                   DIMENSION OF A DIMENSIONED ARRAY
C                   N.B.  IT IS ASSUMED THAT THE DIMENSION NAME IS THE
C                         ARRAY NAME PRECEDED BY THE LETTER J.  IF THE
C                         ARRAY NAME IS 6 LETTERS, THE DIMENSION NAME
C                         SHOULD OMIT THE LAST LETTER.  THE DIMENSION
C                         NAME WILL BE PRINTED USING (NMVAR(I),I=1,6),
C                         AND THE ARRAY NAME USING (NMVAR(I),I=2,7).
C        MSGTYP = 5 THE ARGUMENT BEING CHECKED IS LDSTAK.
C                   NO LONGER USED.
C        MSGTYP = 6 THE ARGUMENT INDICATES THE FIRST DIMENSION OF
C                   AN ARRAY BEING CHECKED AGAINST THE NUMBER OF
C                   UNFIXED PARAMETERS.
C        MSGTYP = 7 THE INPUT VALUE WAS TOO SMALL BASED ON OTHER INPUT
C                   ARGUMENTS, WHERE THE VALUE INDICATES THE
C                   DIMENSION OF A VECTOR.
C                   N.B.  IT IS ASSUMED THAT THE DIMENSION NAME IS THE
C                         ARRAY NAME PRECEDED BY THE LETTER L.  IF THE
C                         ARRAY NAME IS 6 LETTERS, THE DIMENSION NAME
C                         SHOULD OMIT THE LAST LETTER.  THE DIMENSION
C                         NAME WILL BE PRINTED USING (NMVAR(I),I=1,6),
C                         AND THE ARRAY NAME USING (NMVAR(I),I=2,7).
C        MSGTYP = 8 THE INPUT VALUE WAS TOO SMALL BASED ON OTHER INPUT
C                   ARGUMENTS, WHERE THE VALUE INDICATES THE
C                   DIMENSION OF THE VECTORS ACOV AND NLPPA.
C        MSGTYP = 9 THE INPUT VALUE WAS TOO SMALL BASED ON LIMITS
C                   IMPOSED BY STARPAC, WHERE THE VALUE INDICATES THE
C                   DIMENSION OF A VECTOR.
C                   N.B.  IT IS ASSUMED THAT THE DIMENSION NAME IS THE
C                         ARRAY NAME PRECEDED BY THE LETTER L.  IF THE
C                         ARRAY NAME IS 6 LETTERS, THE DIMENSION NAME
C                         SHOULD OMIT THE LAST LETTER.  THE DIMENSION
C                         NAME WILL BE PRINTED USING (NMVAR(I),I=1,6),
C                         AND THE ARRAY NAME USING (NMVAR(I),I=2,7).
C     INTEGER NMIN
C        THE MINIMUM ACCEPTABLE VALUE FOR THE ARGUMENT BEING TESTED.
C     CHARACTER*1 NMSUB(6)
C        THE CHARACTERS OF THE CALLING SUBROUTINES NAME.
C     CHARACTER*1 NMVAR1(8)
C        THE CHARACTERS OF THE NAME OF THE ARGUMENT BEING CHECKED.
C     CHARACTER*1 NMVAR2(8)
C        THE CHARACTERS OF THE NAME OF THE ARGUMENT BEING CHECKED
C        AGAINST.
C     INTEGER NVAL
C        THE INPUT VALUE OF THE ARGUMENT BEING CHECKED.
C
      ERROR = .FALSE.
C
      IF (NVAL .GE. NMIN) RETURN
C
      ERROR = .TRUE.
C
      CALL IPRINT (IPRT)
C
      CALL EHDR(NMSUB, HEAD)
C
      WRITE (IPRT, 1000) (NMVAR1(I), I=1,6), NVAL
C
      GO TO (20, 30, 40, 50, 60, 70, 80, 90, 100), MSGTYP
C
C     PRINT MESSAGE FOR VALUE TOO SMALL BASED ON LIMITS IMPOSED
C     BY STARPAC.
C
   20 WRITE (IPRT, 1010) (NMVAR1(I), I=1,6), NMIN
      RETURN
C
C     PRINT MESSAGE FOR VALUE TOO SMALL BASED ON OTHER INPUT
C     ARGUMENTS.
C
   30 WRITE (IPRT, 1020) (NMVAR1(I), I=1,6), (NMVAR2(I), I=1,8)
      RETURN
C
C     PRINT MESSAGE FOR VALUE TOO SMALL, WHERE VALUE INDICATED THE
C     FIRST DIMENSION OF A DIMENSIONED ARRAY.
C
   40 WRITE (IPRT, 1030) (NMVAR1(I), I=2,7), (NMVAR1(I), I=1,6),
     +   (NMVAR2(I), I=1,8)
      RETURN
C
C     PRINT MESSAGE FOR VALUE TOO SMALL, WHERE VALUE INDICATED THE
C     SECOND DIMENSION OF A DIMENSIONED ARRAY.
C
   50 WRITE (IPRT, 1040) (NMVAR1(I), I=2,7), (NMVAR1(I), I=1,6),
     +   (NMVAR2(I), I=1,8)
      RETURN
C
C     PRINT MESSAGE FOR VALUE TOO SMALL, WHEN ARGUMENT IS LDSTAK.
C
   60 WRITE(IPRT, 1050) NMIN
      RETURN
C
C     PRINT MESSAGE FOR VALUE TOO SMALL, WHERE VALUE INDICATED THE
C     FIRST DIMENSION OF A DIMENSIONED ARRAY CHECK AGAINST THE NUMBER OF
C     UNFIXED PARAMETERS.
C
   70 WRITE (IPRT, 1060) (NMVAR1(I), I=2,7), (NMVAR1(I), I=1,6)
      RETURN
C
C     PRINT MESSAGE FOR VALUE TOO SMALL, WHERE VALUE INDICATED THE
C     DIMENSION OF A VECTOR.
C
   80 WRITE (IPRT, 1070) (NMVAR1(I), I=2,7), (NMVAR1(I), I=1,6),
     +   (NMVAR2(I), I=1,8)
      RETURN
C
C     PRINT MESSAGE FOR VALUE TOO SMALL, WHERE VALUE INDICATED THE
C     DIMENSION OF THE VECTORS ACOV AND NLPPA.
C
   90 WRITE (IPRT, 1080) (NMVAR1(I), I=1,6), (NMVAR2(I), I=1,8)
      RETURN
C
C     PRINT MESSAGE FOR VALUE TOO SMALL, WHERE VALUE INDICATED THE
C     DIMENSION OF A VECTOR.
C
  100 WRITE (IPRT, 1090) (NMVAR1(I), I=2,7), (NMVAR1(I), I=1,6),
     +   NMIN
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (/20H THE INPUT VALUE OF , 6A1, 4H IS , I5, '.')
 1010 FORMAT(
     +   27H THE VALUE OF THE ARGUMENT , 6A1,
     +   34H MUST BE GREATER THAN OR EQUAL TO , I5, '.')
 1020 FORMAT(
     +   27H THE VALUE OF THE ARGUMENT , 6A1,
     +   34H MUST BE GREATER THAN OR EQUAL TO , 8A1, '.')
 1030 FORMAT(
     +   24H THE FIRST DIMENSION OF , 6A1,
     +   30H, AS INDICATED BY THE ARGUMENT/
     +    1X, 6A1, 35H, MUST BE GREATER THAN OR EQUAL TO , 8A1, '.')
 1040 FORMAT(
     +   25H THE SECOND DIMENSION OF , 6A1,
     +   30H, AS INDICATED BY THE ARGUMENT/
     +    1X, 6A1, 35H, MUST BE GREATER THAN OR EQUAL TO , 8A1, '.')
 1050 FORMAT(
     +   55H THE DIMENSION OF THE DOUBLE PRECISION VECTOR DSTAK, AS,
     +   13H INDICATED BY/
     +   54H THE ARGUMENT LDSTAK, MUST BE GREATER THAN OR EQUAL TO,
     +   I5, '.')
 1060 FORMAT(
     +   24H THE FIRST DIMENSION OF , 6A1,
     +   30H, AS INDICATED BY THE ARGUMENT/
     +    1X, 6A1, 34H, MUST BE GREATER THAN OR EQUAL TO,
     +   34H THE NUMBER OF UNFIXED PARAMETERS.)
 1070 FORMAT(
     +   15H THE LENGTH OF , 6A1,
     +   30H, AS INDICATED BY THE ARGUMENT/
     +    1X, 6A1, 35H, MUST BE GREATER THAN OR EQUAL TO , 8A1, '.')
 1080 FORMAT(
     +   29H THE LENGTH OF ACOV AND NLPPA,
     +   30H, AS INDICATED BY THE ARGUMENT/
     +    1X, 6A1, 35H, MUST BE GREATER THAN OR EQUAL TO , 8A1, '.')
 1090 FORMAT(
     +   15H THE LENGTH OF , 6A1,
     +   30H, AS INDICATED BY THE ARGUMENT/
     +    1X, 6A1, 35H, MUST BE GREATER THAN OR EQUAL TO , I6, '.')
C
      END
