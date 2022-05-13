*ERAGTP
      SUBROUTINE ERAGTP (NMSUB, NMVAR, YMMN, NVMX, HEAD, MSGTYP, NV,
     +   NMMIN)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE PRINTS THE ERROR MESSAGES FOR ERAGT AND ERAGTM.
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
      REAL
     +   YMMN
      INTEGER
     +   MSGTYP,NV,NVMX
      LOGICAL
     +   HEAD
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
C     LOGICAL HEAD
C        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
C        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
C        OF HEAD WILL BE CHANGED TO FALSE.
C     INTEGER I
C        AN INDEX ARGUMENT.
C     INTEGER IPRT
C        THE UNIT NUMBER FOR PRINTED OUTPUT.
C     INTEGER MSGTYP
C        THE INDICATOR ARGUMENT FOR THE TYPE OF MESSAGE.
C        IF (MSGTYP.GE.3) THE MESSAGE PRINTED WILL USE NMMIN
C        OTHERWISE IT WILL USE YMMN.
C        IF (MSGTYP = 1 OR 3) NO VIOLATIONS ARE ALLOWED.
C        IF (MSGTYP = 2 OR 4) THE NUMBER OF VIOLATIONS MUST
C                             BE LESS THAN   NVMX   .
C     CHARACTER*1 NMMIN(8)
C        THE NAME OF THE ARGUMENT SPECIFYING THE MINIMUM.
C     CHARACTER*1 NMSUB(6)
C        THE CHARACTERS OF THE CALLING ROUTINES NAME.
C     CHARACTER*1 NMVAR(8)
C        THE CHARACTERS OF THE PARAMETERS NAME.
C     INTEGER NV
C        THE NUMBER OF VIOLATIONS FOUND.
C     INTEGER NVMX
C        THE LARGEST NUMBER OF VIOLATIONS ALLOWED.
C     REAL YMMN
C        THE MINIMUM ACCEPTABLE VALUE.
C
      CALL IPRINT(IPRT)
      CALL EHDR(NMSUB, HEAD)
C
      IF (MSGTYP.LE.2)
     +   WRITE (IPRT, 1000) (NMVAR(I),I=1,6), YMMN, NV
      IF (MSGTYP.GE.3)
     +   WRITE (IPRT, 1005) (NMVAR(I),I=1,6), (NMMIN(I),I=1,8), NV
C
      GO TO (10, 20, 30, 40), MSGTYP
C
   10 WRITE(IPRT, 1010) (NMVAR(I),I=1,6), YMMN
      RETURN
C
   20 WRITE(IPRT, 1020) (NMVAR(I),I=1,6), YMMN, NVMX
      RETURN
C
   30 WRITE(IPRT, 1030) (NMVAR(I),I=1,6), (NMMIN(I),I=1,8)
      RETURN
C
   40 WRITE(IPRT, 1040) (NMVAR(I),I=1,6), (NMMIN(I),I=1,8), NVMX
      RETURN
C
C     FORMAT STATEMENTS
C
 1000 FORMAT (/
     +   31H THE NUMBER OF VALUES IN ARRAY , 6A1,
     +   23H LESS THAN OR EQUAL TO , 1PE14.7, 4H IS , I6, '.')
 1005 FORMAT (/
     +   31H THE NUMBER OF VALUES IN ARRAY , 6A1,
     +   23H LESS THAN OR EQUAL TO , 8A1, 4H IS , I6, '.')
 1010 FORMAT(
     +   25H THE VALUES IN THE ARRAY , 6A1,
     +   26H MUST ALL BE GREATER THAN , 1PE14.7, '.')
 1020 FORMAT(
     +   35H THE NUMBER OF VALUES IN THE ARRAY , 6A1,
     +   23H LESS THAN OR EQUAL TO , 8A1/
     +   19H MUST BE LESS THAN , I5, '.')
 1030 FORMAT(
     +   25H THE VALUES IN THE ARRAY , 6A1,
     +   26H MUST ALL BE GREATER THAN , 1PE14.7, '.')
 1040 FORMAT(
     +   35H THE NUMBER OF VALUES IN THE ARRAY , 6A1,
     +   23H LESS THAN OR EQUAL TO , 8A1/
     +   19H MUST BE LESS THAN , I5, '.')
C
      END
