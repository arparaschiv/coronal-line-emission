;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Document name: redp.pro
; Created by:    Philip Judge, September 12, 2000
;
; Last Modified: Tue Nov 29 14:28:00 2005 by judge (Philip Judge) on niwot.local
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
PRO redp,q,u,x,y,nmx,x1,y1,q1,u1
aa = size(q)
x1 = x
y1 = y
u1 = u
q1 = q
IF(aa(1) GT nmx) THEN BEGIN 
   sx = nmx
   sy = round(float(sx)/aa(1)*aa(2))
   q1=congrid(q,sx,sy)
   u1=congrid(u,sx,sy)
   xx = findgen(nmx)/(nmx-1)
   yy = findgen(sy)/(sy-1)
   x1 = min(x)+xx*(max(x)-min(x))
   y1 = min(y)+yy*(max(y)-min(y))
ENDIF
return
END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of 'redp.pro'.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
