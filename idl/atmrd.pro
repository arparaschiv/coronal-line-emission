;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Document name: atmrd.pro
; Created by:    Philip Judge, August 31, 2000
;
; Last Modified: Mon Jan 30 18:36:42 2006 by judge (judge) on macniwot.local
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
; Read in atmospheric model
;
PRO atmrd, atm, file = file
IF(n_elements(file) EQ 0) THEN file = 'ATMOS'
openr,lu,/get,file,/f77_unformatted
message,/inf,'reading '+file+' ...'
;

gxmax = 0. &  gxmin = 0. & ngx = 10 
gymax = 0. &  gymin = 0. & ngy = 12
gzmax = 0. &  gzmin = 0. & ngz = 12
cname = '      '
forrd,lu,cname
forrd,lu,gxmin,gxmax,ngx
forrd,lu,gymin,gymax,ngy
forrd,lu,gzmin,gzmax,ngz
IF(ngx+ngy+ngz EQ 0) THEN BEGIN 
   message,/inf,'Integers not written correctly'
   read,'Enter NGX, NGY, NGZ',ngx,ngy,ngz
ENDIF
;
x = gxmin+ (gxmax-gxmin) *findgen(ngx)/(ngx-1)
y = gymin+ (gymax-gymin) *findgen(ngy)/(ngy-1)
z = gzmin+ (gzmax-gzmin) *findgen(ngz)/(ngz-1)
bx = fltarr(ngx,ngy,ngz)
by = bx
bz = bx
nne = bx
te = bx
v = bx
vt = bx
ix = 0l &  iy = 0l &  iz = 0l 
rbx = 0. &  rby = 0. &  rbz = 0. &  rnne = 0. &  rte = 0. &  rv = 0. &  rvt = 0.
WHILE NOT eof(lu) DO BEGIN 
   forrd,lu,ix,iy,iz
   forrd,lu,rbx,rby,rbz,rnne,rte,rv,rvt
   ix = ix-1
   iy = iy-1
   iz = iz-1
   bx(ix,iy,iz) = rbx
   by(ix,iy,iz) = rby
   bz(ix,iy,iz) = rbz
   nne(ix,iy,iz) = rnne
   te(ix,iy,iz) = rte
   v(ix,iy,iz) = rv
   vt(ix,iy,iz) = rvt
ENDWHILE
free_lun,lu
;
;
readme = (['variables are ',$
                    'file:     file containing atmospheric parameters',$
                    'model:    name of fortran routine used to determine coronal parameters',$
                    'x,y,z:    positions in solar radii',$
                    '          (x is along LOS, y is E-W, z N-S in observer''s frame)',$
                    'bx,by,bz: x,y,z-components of magnetic field in G',$
                    'nne:      electron density in /cm3',$
                    'te:       electron temperature in K',$
                    'v:        LOS velocity in km/s (+ve is a red-shift)',$
                    'vt:       "turbulent" velocity in km/s'])

atm = {file:file,model:cname,x:x,y:y,z:z,bx:bx,by:by,bz:bz,nne:nne,te:te,v:v,vt:vt,$
       readme:readme}

;hdata=ptr_new(nne)
;message,' use IDL> slicer3,hdata',/inf
return
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of 'atmrd.pro'.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
