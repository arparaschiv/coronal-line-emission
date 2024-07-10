C
C PURPOSE: LEFT JUSTIFIES THE STRING TEXT
C
C INPUTS:
C
C OUTPUTS:
C
C COMMON:
C
C COMMENTS: OCTOBER 6, 1999, P. JUDGE
C 
      SUBROUTINE LJUST(TEXT)
C
C  LEFT JUSTIFIES THE STRING TEXT
C
      INCLUDE 'PREC'
      CHARACTER*(*) TEXT
C
      L=LEN(TEXT)
      DO 100 J=1,L
        IF(TEXT(J:J).NE.' ') GOTO 200
  100 CONTINUE
  200 CONTINUE
      DO 300 I=1,L
        IF(J.LE.L) THEN
          TEXT(I:I)=TEXT(J:J)
        ELSE
          TEXT(I:I)=' '
        ENDIF
        J=J+1
  300 CONTINUE
C
      RETURN
      END
C
C************************************************************************
C
