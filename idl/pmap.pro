;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Document name: pmap.pro
; Created by:    Philip Judge, December 6, 2005
;
; Last Modified: Tue Dec  6 17:44:05 2005 by judge (Philip Judge) on niwot.local
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
FUNCTION pmap,s,relative = relative
;   
; makes a linear polarization magnitude map
; from the map structure s read by outrd,s
;   
p = s(0)
p.id = 'P '
p.data = sqrt(s(1).data*s(1).data+ s(2).data*s(2).data)
IF(n_elements(relative) NE 0) THEN BEGIN 
   p.data = p.data/s(0).data
   p.id = 'P/I '
endif
return,p
end   

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of 'pmap.pro'.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
