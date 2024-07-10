;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Document name: qu2vect.pro
; Created by:    Philip Judge, September 6, 2000
;
; Last Modified: Tue Apr 17 18:12:02 2007 by judge (judge) on edlen.local
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
FUNCTION qu2vect,q,u
; +  
; rotates q and u vectors from CLE, computed there using a ref direction 
; parallel to the z-axis [N-S], to a polarization vector in the frame of 
; (y,z)  [ E-W, N-S] in the plane of the sky. The input numbers q,u are
; replaced with these components.   
;
; the CLE routine corona.f contains the definition of the 
;   reference direction.
; -  
IF(q EQ 0. AND u eq 0.) THEN return,[0.,0.]
;
alpha = -0.5*atan(u,q)
p = sqrt(q*q+u*u)
alpha = -alpha+!dpi/2
return,p*[cos(alpha),sin(alpha)]
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of 'qu2vect.pro'.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
