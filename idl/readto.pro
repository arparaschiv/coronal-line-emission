;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Document name: readto.pro
; Created by:    Philip Judge, November 21, 2002
;
; Last Modified: Thu Nov 21 10:38:34 2002 by judge (Philip Judge) on judgepc.hao.ucar.edu
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
PRO readto,str,lu,s
l = strlen(str)
s = ''
while strmid(s,0,l) NE str DO BEGIN 
   readf,lu,s
ENDWHILE

return
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of 'readto.pro'.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
