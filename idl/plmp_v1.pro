;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Document name: plmp.pro
; Created by:    Philip Judge, December 8, 2005
;
; Last Modified: Thu Jul 13 17:10:47 2006 by judge (judge) on macniwot.local
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
PRO plmp,map,dran = dran,_extra = _extra,top = top,bottom = bottom,left = left,$
         xspace = xspace
;
; get data ranges  
;
IF(n_elements(dran) EQ 0) THEN BEGIN 
      dd = median(sqrt(map.data*map.data))
      dran = [0,dd]*3
      titl1 = map.id
      iv = strpos(titl1,'V')
      IF(iv GE 0) THEN dran = [-dd,dd]*3
      ii = strpos(titl1,'Stokes I')
      IF(ii GE 0) THEN dran = dran*5./3
ENDIF
;
; plot the map
;
noerase = !p.multi(0) NE 0
plot_map,map,dran = dran,title = title,_extra = _extra,noerase = noerase
;
;
; color bar
;
IF(n_elements(xspace) EQ 0) THEN xspace = 1
x = 1.04
y = 0.58
;
;x=0.76
;y=0.05
;
IF(n_elements(left) NE 0) THEN x = 0.2
IF(n_elements(top) NE 0) THEN y = 0.7
dx = 0.05
dy = 0.2
;
csize = !p.charsize
IF(csize EQ 0) THEN csize = 1
;m = !p.multi
;n = m(1)
width = (!x.window(1)-!x.window(0))
x = !x.window(0) + x*width
dx = dx*width
x1 = x+dx
x0 = x-dx/3
x2 = x1+dx*csize*xspace
;n = m(2)
height = (!y.window(1)-!y.window(0))
y = !y.window(0) + y*height
dy = dy*height
y1 = y+dy
y0 = y-dy/20
y2 = y1+dy/20
print,'x0,x2,y0,y2',x0,x2,y0,y2
polyfill,[x0,x2,x2,x0,x0],[y0,y0,y2,y2,y0],color=0,/norm
dmin = min(dran) &  dmax = max(dran)
color_bar,x,x1,y,y1,max = dmax,min = dmin,top = top,$
   title = '',/right,charsize = charsize,/norm,color = -1 ;!d.n_colors-1
return
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of 'plmp.pro'.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
