*ENFFT
      SUBROUTINE ENFFT(NMSUB, NFFT, NDIV, N, LYFFT, NFFT2, HEAD, ERROR)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS ROUTINE CHECKS WHETHER THE VALUE NFFT IS SUCH THAT NFFT-2 IS
C     DIVISIBLE BY NDIV AND HAS NO PRIME FACTORS GREATER THAN 23, AND
C     THE PRODUCT OF THE SQUARE FREE PRIME FACTORS OF NFFT - 2 DO NOT
C     EXCEED 209, I.E., THE VALUE OF NFFT MEETS THE REQUIREMENTS OF
C     THE EXTENDED LENGTH OF THE SERIES REQUIRED FOR ANY ROUTINE
C     USING THE SINGLETON FFT PROVIDING THE PROPER VALUE OF NDIV
C     IS CHOSEN.
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
     +   LYFFT,N,NDIV,NFFT,NFFT2
      LOGICAL
     +   ERROR,HEAD
C
C  ARRAY ARGUMENTS
      CHARACTER
     +   NMSUB(6)*1
C
C  LOCAL SCALARS
      INTEGER
     +   IPRT,NFFT1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL EHDR,IPRINT,SETESL
C
C  INTRINSIC FUNCTIONS
      INTRINSIC MAX
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
C     INTEGER IPRT
C        THE UNIT NUMBER FOR PRINTED OUTPUT.
C     INTEGER LYFFT
C        THE LENGTH OF THE VECTOR CONTAINING THE SERIES TO BE EXTENDED.
C     INTEGER N
C        THE ACTUAL NUMBER OF OBSERVATIONS IN THE SERIES.
C     CHARACTER*1 NMSUB(6)
C        THE CHARACTERS OF THE CALLING SUBROUTINES NAME.
C     INTEGER NFFT
C        THE USER SUPPLIED EXTENDED SERIES LENGTH.
C     INTEGER NFFT1
C        THE MAXIMUM OF NFFT AND N+2.
C     INTEGER NFFT2
C        THE SMALLEST EXTENDED SERIES LENGTH WHICH EQUALS OR
C        EXCEEDS NFFT AND WHICH MEETS THE REQUIREMENTS OF
C        SINGLETONS FFT CODE.
C
      ERROR = .FALSE.
      CALL IPRINT (IPRT)
C
      IF (NFFT .GE. N+2) GO TO 20
C
C     PRINT WARNING
C
      CALL EHDR(NMSUB, HEAD)
C
      WRITE (IPRT, 1050) N
C
   20 CONTINUE
      NFFT1 = MAX(NFFT, N+2)
      CALL SETESL(NFFT1-2, NDIV, NFFT2)
C
      IF (NFFT .EQ. NFFT2) RETURN
C
C     PRINT WARNING
C
      CALL EHDR(NMSUB, HEAD)
C
      WRITE (IPRT, 1020) NFFT, NFFT2
C
      IF (NFFT .GT. LYFFT) GO TO 40
C
      WRITE (IPRT, 1030) NFFT2
      RETURN
C
   40 CONTINUE
C
      ERROR = .TRUE.
C
      WRITE (IPRT, 1040) NFFT2, LYFFT
      RETURN
C
C     FORMAT STATEMENTS
C
 1020 FORMAT (/
     +   40H THE INPUT VALUE OF THE PARAMETER NFFT (, I5,
     +   15H) DOES NOT MEET/
     +   51H THE REQUIREMENTS OF SINGLETONS FFT CODE.  THE NEXT,
     +   13H LARGER VALUE/
     +   15H WHICH DOES IS , I5, '.')
 1030 FORMAT (/
     +   11H THE VALUE , I5, 37H WILL BE USED FOR THE EXTENDED SERIES,
     +    8H LENGTH.)
 1040 FORMAT (/
     +   20H HOWEVER, THE VALUE , I5, 27H EXCEEDS THE LENGTH LYFFT (,
     +   I5, 8H) OF THE/
     +   58H VECTOR YFFT, AND THEREFORE CANNOT BE USED AS THE EXTENDED/
     +   43H SERIES LENGTH WITHOUT REDIMENSIONING YFFT.)
 1050 FORMAT (/
     +   56H THE EXTENDED SERIES LENGTH (NFFT) MUST EQUAL OR EXCEED,/
     +   45H THE NUMBER OF OBSERVATIONS IN THE SERIES (N=, I5,
     +    9H) PLUS 2.)
C
      END
