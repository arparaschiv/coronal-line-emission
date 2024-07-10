FUNCTION XLINPOT,Y,Z,ycen,zcen,bcen
;
;       VECTOR POTENTIAL AX FROM A LINE CURRENT CENTERED
;
XLINPOT=0.
N=n_elements(ycen)
FOR J=0,N-1 DO BEGIN 
   YD=Y-YCEN(J)
   ZD=Z-ZCEN(J)
   R=SQRT(YD*YD+ZD*ZD)
   XLINPOT= XLINPOT -BCEN(J)*ALOG(R)
ENDFOR
RETURN,XLINPOT
END

FUNCTION PLASMOID,P
;
; compute plasmoid from
; inputs to cle
;

ax=p
ax.data=0  

n=6
yc=[0.85,1.05, 1.09, 1.13, 1.17, 1.21]
zc=[0.,0., 0., 0., 0., 0.]
;b=[1., 1., .3, .3, .1]
b=[20.,1., 1., 1., 1., 1.]
Y=GET_MAP_XP(AX)
Z=GET_MAP_YP(AX)

AX.data=XLINPOT(Y,Z,yc,zc,b)
RETURN,AX

END

