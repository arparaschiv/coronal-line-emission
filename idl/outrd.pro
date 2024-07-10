;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Document name: outrd.pro
; Created by:    Philip Judge, November 19, 2002
;
; Last Modified: Wed Mar 15 17:37:18 2006 by judge (judge) on macniwot.local
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
;
PRO outrd,s,kr = kr,file = file,no_rot=no_rot,silent = silent
;
;+   
; reads output of the fortran program cle.
; kr is the transition index (the atom contains several radiative
;    transitions, kr is the    
;
IF(n_elements(file) EQ 0) THEN file = 'OUT'
IF(n_elements(kr) EQ 0) THEN kr = 0
IF(n_elements(no_rot) EQ 0) THEN no_rot = 0
;
notes = transpose(['To read data, do e.g. ',$
                     '  IDL> outrd,s,file=''OUT''',$
                               ' ',$
                     'Output from cle (contact: judge@ucar.edu)',$
                     '  s is a map structure, where ',$
                     '  s(0) is the intensity I integrated over the a line',$
                     '  s(1) is Q or Px integrated over the line',$
                     '  s(2) is U or Py integrated over the line',$
                     '  s(3) is V integrated over the line (using weights of -1,+1',$
                     '          for wavelengths < and > line center at rest)',$
                     '  s(4) is Vmag=V but computed using the magnetograph formula',$
                               ' ',$
                     '-> maps are manipulated using mapping software on SSW, e.g.',$
                     '  IDL> chkarg,''plot_map''',$
                     '  IDL> plot_map,s(3) ; plots a map of Stokes V'])
full_notes = transpose(['  s.q are delta-wavelengths across line profile [Angstrom]',$
                             '  s.full are stokes vectors full(q,x,y)',$
                             '  units are erg/cm^2/s/Angstrom'])
;
IF(n_params(0) EQ 0) THEN BEGIN 
   print,notes
   print,''
   print,'If IWLINE was NE 0 on running cle, full Stokes vectors are read and stored. Then'
   print, full_notes
   return
ENDIF
;
openr,lu,/get,file
;
IF(n_elements(silent) EQ 0) THEN message,/inf,'reading '+file+' ...'
str = ''
readto,'QNORM',lu,s
QNORM = float(getwrd(s,1,1))
readto,'IWLINE',lu,s
IWLINE = fix(getwrd(s,1,1))
readto,'CRTN',lu,s
CRTN = getwrd(s,1,1)
readto,'1 GRID',lu
readto,'DIRECTION ',lu
readf,lu,str
s1 = getwrd(str,1,1,delim=':')    
gxmin = float(getwrd(s1,0,0))
gxmax = float(getwrd(s1,1,1))
ngx = long(getwrd(s1,2,2))
readf,lu,str
s1 = getwrd(str,1,1,delim=':')    
gymin = float(getwrd(s1,0,0))
gymax = float(getwrd(s1,1,1))
ngy = long(getwrd(s1,2,2))
readf,lu,str
s1 = getwrd(str,1,1,delim=':')    
gzmin = float(getwrd(s1,0,0))
gzmax = float(getwrd(s1,1,1))
ngz = long(getwrd(s1,2,2))
gx = gxmin+(gxmax -gxmin)/(ngx-1)*findgen(ngx)
gy = gymin+(gymax -gymin)/(ngy-1)*findgen(ngy)
gz = gzmin+(gzmax -gzmin)/(ngz-1)*findgen(ngz)
;
cornotes = ''
IF(crtn EQ 'FFL') THEN BEGIN 
   readto,'1 CORONAL ROUTINE',lu
   cornotes = 'FFL: '
   readf,lu,str
   cornotes = 'FFL: '+str
   readf,lu,str
   cornotes = cornotes+' ' + str
   readf,lu,str
   cornotes = cornotes+' ' + str
   cornotes = strcompress(cornotes)
ENDIF
;
readto,'1 STOKES',lu
readto,' TOTAL NUMBER OF LINES ',lu,s
nline = long(getwrd(s,/last))
alamb = dblarr(nline)
nq = intarr(nline)
q = fltarr(nline,100)
nqmax = 0
FOR i = 0,nline-1 DO BEGIN
   readf,lu,str
   alamb(i) = double(getwrd(str,/last))
   readf,lu,str
   nq(i) = fix(getwrd(str,/last))
   nqmax = nqmax >  nq(i)
   dum = fltarr(nq(i))
   readf,lu,dum
   q(i,0:nq(i)-1) = dum
ENDFOR
;
q = q(*,0:nqmax-1)

emerg = fltarr(nline,5,ngy,ngz)
IF(iwline GT 0) THEN BEGIN 
   full =  fltarr(nline,5,nqmax,ngy,ngz)
   readto,' FULL EMERGENT STOKES PROFILES',lu
ENDIF ELSE BEGIN 
   readto,' FREQUENCY-INTEGRATED',lu
ENDELSE
;
count=0
FOR i = 0,ngy-1 DO BEGIN
   FOR j = 0,ngz-1 DO BEGIN
      FOR m = 0,nline-1 DO BEGIN 
          IF(iwline GT 0) THEN BEGIN 
             fdum = fltarr(nq(m))
             FOR ll = 0,4 DO BEGIN 
                readf,lu,fdum; ,form = '(10e11.3)'; read detailed profiles
                full(m,ll,0:nq(m)-1,i,j) = fdum
             endfor
          ENDIF
          dum = fltarr(5)
          readf,lu,dum
          emerg(m,*,i,j) = dum
          q1 = emerg(m,1,i,j)
          u1 = emerg(m,2,i,j)
	  xx=[q1,u1]
	  if(not no_rot) then xx = qu2vect(q1,u1)  
          emerg(m,1,i,j) = xx(0)
	  emerg(m,2,i,j) = xx(1)
          count+=1
;          print,i,j,m,count,dum,form='(4(i4,1x),5(1x,e9.2))'
       ENDFOR
   ENDFOR
ENDFOR
free_lun,lu
;
s = size(emerg)
IF(s(3) EQ 1) THEN BEGIN 
   i = reform(emerg(kr,0,*,*),1,s(4))
   px = reform(emerg(kr,1,*,*),1,s(4))
   py = reform(emerg(kr,2,*,*),1,s(4))
   v = reform(emerg(kr,3,*,*),1,s(4))
   vmag = reform(emerg(kr,4,*,*),1,s(4))
ENDIF ELSE BEGIN 
   i = reform(emerg(kr,0,*,*))
   px = reform(emerg(kr,1,*,*))
   py = reform(emerg(kr,2,*,*))
   v = reform(emerg(kr,3,*,*))
   vmag = reform(emerg(kr,4,*,*))
ENDELSE

p = sqrt(px*px+py*py)
;
; output to maps
;
xc = (gymax+gymin)/2. 
dx =  (gymax-gymin)/(ngy)
yc = (gzmax+gzmin)/2. 
dy =  (gzmax-gzmin)/(ngz)
;
x = gymin + findgen(ngy)*dx
y = gzmin + findgen(ngz)*dy
;
; units
;
punit = 955.                     ;arcsec
unit = 'arcseconds'
punit = 1.   ;solar radii
unit = 'solar radii'
xc = xc*punit
dx = dx*punit
yc = yc*punit
dy = dy*punit
x = x*punit
y = y*punit
;
im = make_map(i,xc = xc,yc = yc,id = 'Stokes I',dx = dx,dy = dy,yun = unit,xun = unit)
add_prop,im,lambda = alamb(kr)
add_prop,im,dunits = 'erg/cm^2/s'
add_prop,im,notes = notes(0,2:*)
add_prop,im,cornotes = cornotes
IF(iwline GT 0) THEN BEGIN 
   add_prop,im,q =    reform(q(kr,0:nq(kr)-1))
   fu = reform(full(kr,0,*,*,*))
   add_prop,im,full = fu
   add_prop,im,full_notes = full_notes
ENDIF
;
;
s = replicate(im,5)
sx='Px'
sy='Py'
if(no_rot) then begin
  sx='Q (ref=z-axis)'
  sy='U (ref=z-axis)'
endif
;
s(0) = im
s(1).data = px &  s(1).id = 'Stokes '+sx
s(2).data = py &  s(2).id = 'Stokes '+sy
s(3).data = v & s(3).id = 'Stokes V'
s(4).data = vmag & s(4).id = 'Stokes V (magnetograph formula)'

IF(n_elements(silent) EQ 0) THEN BEGIN 
   print,''
   print,s(0).notes
   print,''
ENDIF
IF(IWLINE GT 0) THEN BEGIN
   FOR i = 1,4 DO s(i).full = reform(full(kr,i,*,*,*))
   IF(n_elements(silent) EQ 0) THEN $
      message,/inf,'Full stokes profiles were also outputted, stored e.g., in s.q, s.full'
   IF(n_elements(silent) EQ 0) THEN print,s(0).full_notes
ENDIF

;
;
return
END
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of 'outrd.pro'.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
