*STATER
      SUBROUTINE STATER(NMSUB, WT, N, LDSTAK, WTS, NNZW, STACK, IERR)
C
C     LATEST REVISION  -  03/15/90  (JRD)
C
C     THIS SUBROUTINE CHECKS INPUT PARAMETERS TO THE USER
C     CALLABLE MEMBERS OF THE STAT FAMILY OF ROUTINES
C     FOR ERRORS AND REPORTS ANY THAT IT FINDS, BESIDES
C     RETURNING A FLAG INDICATING THAT ERRORS HAVE BEEN
C     FOUND.
C
C     WRITTEN BY - JOHN E. KOONTZ
C                  STATISTICAL ENGINEERING DIVISION
C                  NATIONAL BUREAU OF STANDARDS
C
C     CREATION DATE  -  MAY 17, 1982
C
C
C  VARIABLE DECLARATIONS
C
C  SCALAR ARGUMENTS
      INTEGER
     +   IERR,LDSTAK,N,NNZW
      LOGICAL
     +   STACK,WTS
C
C  ARRAY ARGUMENTS
      REAL
     +   WT(*)
      CHARACTER
     +   NMSUB(6)*1
C
C  LOCAL SCALARS
      INTEGER
     +   LDSMIN,NZW
      LOGICAL
     +   HEAD,IER1,IER2,IER3
C
C  LOCAL ARRAYS
      CHARACTER
     +   LLDS(8)*1,LN(8)*1,LTHREE(8)*1,LWT(8)*1
C
C  EXTERNAL SUBROUTINES
      EXTERNAL EISGE,ERVWT,LDSCMP
C
C     VARIABLE DEFINITIONS (ALPHABETICALLY)
C
C     LOGICAL HEAD
C        A FLAG INDICATING WHETHER THE HEADING SHOULD BE PRINTED
C        (TRUE) OR NOT (FALSE).  IF A HEADING IS PRINTED, THE VALUE
C        OF HEAD WILL BE CHANGED TO FALSE.
C     INTEGER IERR
C        OUTPUT PARAMETER.  A FLAG INDICATING WHETHER OR
C        NOT AN ERROR HAS BEEN FOUND.  0 = OK, 1 = ERROR.
C     LOGICAL IER1
C        TRUE IF N .LT. 3
C     LOGICAL IER2
C        TRUE IF LDSTAK .LT. (N + 13)/2.0E0
C     LOGICAL IER3
C        TRUE IF SOME WT .LT. 0.0E0 OR NNZW .LT. 3
C     INTEGER LDSMIN
C        MINIMUM LENGTH OF FRAMEWORK AREA IN DOUBLE
C        PRECISION ELEMENTS.
C     INTEGER LDSTAK
C        INPUT PARAMETER.  THE NUMBER OF LOCATIONS PROVIDED IN
C        THE FRAMEWORK AREA.
C     CHARACTER*1 LLDS(8), LN(8), LTHREE(8), LWT(8)
C        THE ARRAY(S) CONTAINING THE NAME(S) FO THE VARIALBE(S) CHECKED
C        FOR ERRORS
C     INTEGER N
C        INPUT PARAMETER.  THE NUMBER OF ELEMENTS IN Y AND WT.
C     CHARACTER*1 NMSUB(6)
C        THE NAME OF THE CALLING SUBROUTINE
C     INTEGER NNZW
C        OUTPUT PARAMETER.  IF WTS, THEN SET EQUAL TO THE
C        NUMBER OF VALUES IN WT WHICH ARE POSITIVE.  ELSE,
C        UNDEFINED.
C     INTEGER NZW
C        THE NUMBER OF ZERO WEIGHTS.
C     LOGICAL STACK
C        A FLAG INDICATING WHETHER THIS ROUTINE USES THE STACK (TRUE)
C        OR NOT (FALSE).
C     REAL WT(N)
C        INPUT PARAMETER.  THE VECTOR OF WEIGHTS CORRESPONDING
C        TO THE VECTOR Y.
C     LOGICAL WTS
C        INPUT PARAMETER.  A FLAG INDICATING WHETHER OR NOT
C        THERE IS REALLY A VECTOR WT (TRUE), OR ONLY A DUMMY PARAMETER
C        (FALSE).
C
C     INITIALIZE NAME VECTORS
C
      DATA LLDS(1), LLDS(2), LLDS(3), LLDS(4), LLDS(5), LLDS(6),
     +  LLDS(7), LLDS(8) / 'L', 'D', 'S', 'T', 'A', 'K', ' ', ' '/
      DATA     LN(1),    LN(2),    LN(3),    LN(4),    LN(5),    LN(6),
     +  LN(7), LN(8)  / 'N', ' ', ' ', ' ', ' ', ' ', ' ', ' '/
      DATA LTHREE(1),LTHREE(2),LTHREE(3),LTHREE(4),LTHREE(5),LTHREE(6),
     +  LTHREE(7), LTHREE(8) / 'T', 'H', 'R', 'E', 'E', ' ', ' ', ' '/
      DATA    LWT(1),   LWT(2),   LWT(3),   LWT(4),   LWT(5),   LWT(6),
     +  LWT(7), LWT(8) / 'W', 'T', ' ', ' ', ' ', ' ', ' ', ' '/
C
C     INITIALIZE ERROR FLAGS
C
      IER1 = .FALSE.
      IER2 = .FALSE.
      IER3 = .FALSE.
C
      IERR = 0
C
      HEAD = .TRUE.
C
C     CHECK TO SEE THAT THERE ARE AT LEAST THREE DATA POINTS.
C
      CALL EISGE(NMSUB, LN, N, 3, 2, HEAD, IER1, LTHREE)
C
C     CHECK TO SEE THAT AN AMOUNT OF WORK AREA EQUAL
C     IN LENGTH TO THE REQUIREMENTS OF THE PERMUTATION
C     VECTOR WILL BE AVAILABLE.
C
      IF (STACK) THEN
         CALL LDSCMP(1, 0, N, 0, 0, 0, 'S', 0, LDSMIN)
         CALL EISGE(NMSUB, LLDS, LDSTAK, LDSMIN, 9, HEAD, IER2, LLDS)
      END IF
C
C     IF THERE ARE WEIGHTS
C     CHECK TO SEE THAT AT LEAST THREE DATA ITEMS HAVE NONZERO WEIGHTS.
C
      NNZW = N
      IF (WTS) THEN
         CALL ERVWT(NMSUB, LWT, WT, N, 3, HEAD, NNZW, NZW, 1, IER3,
     +              LTHREE)
      END IF
C
C     SEE IF ANY ERRORS WERE FOUND.
C
      IF (IER1 .OR. IER2 .OR. IER3) IERR = 1
      RETURN
      END
