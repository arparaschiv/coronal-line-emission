C
C PURPOSE: FINDS NEXT WORD IN TEXT FROM INDEX K0
C
C INPUTS:
C     text   input string
c     k0     beginning search position
C OUTPUTS:
c     k1     index position of beginning of word
c     k2     index position of end of word
C COMMON:
C
C COMMENTS: October 6, 1999, P. JUDGE
C 
      SUBROUTINE GETWRD(TEXT,K0,K1,K2)
C
C  FINDS NEXT WORD IN TEXT FROM INDEX K0. NEXT WORD IS TEXT(K1:K2)
C  THE NEXT WORD STARTS AT THE FIRST ALPHANUMERIC CHARACTER AT K0
C  OR AFTER. IT ENDS WITH THE LAST ALPHANUMERIC CHARACTER IN A ROW
C  FROM THE START
C
      INCLUDE 'PREC'
      INTEGER MSEPAR
      PARAMETER (MSEPAR=7)
      CHARACTER*(*) TEXT
      CHARACTER SEPAR(MSEPAR)
      INTEGER K0,K1,K2,I,J
      DATA SEPAR/' ','(',')','=','*','/',','/
C
      K1=0
      DO 400 I=K0,LEN(TEXT)
        IF(K1.EQ.0) THEN
          DO 100 J=1,MSEPAR
            IF(TEXT(I:I).EQ.SEPAR(J)) GOTO 200
  100     CONTINUE
          K1=I
C
C  NOT START OF WORD
C
  200     CONTINUE
        ELSE
          DO 300 J=1,MSEPAR
            IF(TEXT(I:I).EQ.SEPAR(J)) GOTO 500
  300     CONTINUE
        ENDIF
  400 CONTINUE
C
C  NO NEW WORD. RETURN K1=K2=0
C
      K1=0
      K2=0
      GOTO 999
C
C  NEW WORD IN TEXT(K1:I-1)
C
  500 CONTINUE
      K2=I-1
C
  999 CONTINUE
      RETURN
      END
C
C***********************************************************************
C
