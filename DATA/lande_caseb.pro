;--------------------------------------------------------------
; This returns the Lande factor in the case b) 
;--------------------------------------------------------------
function lande_caseb, J, N, S, Lambda
	 return, (Lambda^2/(2*N*(N+1))*(J*(J+1)+N*(N+1)-S*(S+1))+J*(J+1)-N*(N+1)+S*(S+1)) / (J*(J+1))
end

;--------------------------------------------------------------
; This returns the effective Lande factor in the case b) 
;--------------------------------------------------------------
function effec_lande_caseb, J, N, S, Lambdau, Lambdal
	 lu = lande_caseb(J, N, S, Lambdau)
	 ll = lande_caseb(J, N, S, Lambdal)
	 	 
end

