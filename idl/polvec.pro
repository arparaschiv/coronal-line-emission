;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Document name: polvec.pro
; Created by:    Philip Judge, December 8, 2005
;
; Last Modified: Tue Apr 25 11:05:03 2006 by judge (judge) on macniwot.local
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
PRO polvec,s,npts = npts,col = col,xran = xran,yran = yran
IF(n_elements(npts) EQ 0) THEN npts = 30
leng = 1.1
delta = 0.
x = get_map_xp(s(0))
y = get_map_yp(s(0))
qq = s(1).data
uu = s(2).data
;
r = sqrt(x*x+y*y)
disk = where(r LE 1.0,k)
IF(k NE 0) THEN BEGIN
   qq(disk) = !values.f_nan
   uu(disk) = !values.f_nan
ENDIF

   
redp,qq,uu,x,y,npts,x1,y1,q1,u1
IF(n_elements(xran) NE 0) THEN BEGIN 
   k = where (x1 GT xran(0) AND x1 LT xran(1))
   x1 = x1(k)
   u1 = u1(k,*)
   q1 = q1(k,*)
ENDIF
;   
IF(n_elements(yran) NE 0) THEN BEGIN 
   k = where ( y1 GT yran(0) AND y1 LT yran(1))
   y1 = y1(k)
   u1 = u1(*,k)
   q1 = q1(*,k)
ENDIF
;
pline,q1,u1,leng=leng,x1+delta,y1+delta,/over,/clip,/nohead,col = col
return
END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of 'polvec.pro'.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
