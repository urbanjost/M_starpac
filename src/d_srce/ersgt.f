*ERSGT
      SUBROUTINE ERSGT(NMSUB, NMVAR, VAL, VALMN, MSGTYP, HEAD, ERROR,
     +   NMMIN)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THE ROUTINE CHECKS WHETHER THE VALUE   VAL   IS GREATER THAN
C        VALMN   , AND PRINTS A DIAGNOSTIC IF IT IS NOT.
C
C     WRITTEN BY  -  JANET R. DONALDSON
C                    STATISTICAL ENGINEERING DIVISION
C                    NATIONAL BUREAU OF STANDARDS, BOULDER, COLORADO
C
C     CREATION DATE  -  JUNE 10, 1982
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      DOUBLE PRECISION
     +   VAL,VALMN
      INTEGER
     +   MSGTYP
      LOGICAL
     +   ERROR,HEAD
C
C  ARRAY ARGUMENTS
      CHARACTER
     +   NMMIN(8)*1,NMSUB(6)*1,NMVAR(8)*1
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
C        PRINTED, WHERE IF ERROR IS .TRUE. AND
C        MSGTYP = 1 THE INPUT VALUE WAS TOO SMALL BASED
C                   ON LIMITS IMPOSED BY STARPAC
C        MSGTYP = 2 THE INPUT VALUE WAS TOO SMALL BASED ON
C                   LIMITS BASED ON OTHER INPUT ARGUMENTS.
C     CHARACTER*1 NMMIN(8)
C        THE NAME OF THE ARGUMENT SPECIFYING THE MINIMUM.
C     CHARACTER*1 NMSUB(6)
C        THE CHARACTERS OF THE CALLING SUBROUTINES NAME.
C     CHARACTER*1 NMVAR(8)
C        THE CHARACTERS OF THE ARGUMENTS NAME.
C     DOUBLE PRECISION VAL
C        THE INPUT VALUE OF THE ARGUMENT BEING CHECKED.
C     DOUBLE PRECISION VALMN
C        THE MINIMUM VALUE THE ARGUMENT CAN VALIDLY HAVE.
C
      ERROR = .FALSE.
C
      IF (VAL .GT. VALMN) RETURN
C
      ERROR = .TRUE.
      CALL IPRINT (IPRT)
      CALL EHDR(NMSUB, HEAD)
C
      WRITE (IPRT, 1000) (NMVAR(I),I=1,6), VAL
C
      GO TO (10, 20), MSGTYP
C
   10 WRITE (IPRT, 1010) (NMVAR(I),I=1,6), VALMN
      RETURN
C
   20 WRITE (IPRT, 1020) (NMVAR(I),I=1,6), (NMMIN(I),I=1,8)
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (/20H THE INPUT VALUE OF , 6A1, 4H IS , G15.8, '.')
 1010 FORMAT(
     +   27H THE VALUE OF THE ARGUMENT , 6A1,
     +   22H MUST BE GREATER THAN , G21.14, '.')
 1020 FORMAT(
     +   27H THE VALUE OF THE ARGUMENT , 6A1,
     +   22H MUST BE GREATER THAN , 8A1, '.')
C
      END
