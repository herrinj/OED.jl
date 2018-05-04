
import jInv.ForwardShare.getSensMatSize

function  getSensMatSize(pFor::Union{UnconstrainedOEDParam, EqconstrainedOEDParam, IneqconstrainedOEDParam})
	 (m,n,nrhs) = getProblemSize(pFor.Design)
	 return n*nrhs,getNumberOfDesignParameters(pFor.Design)
end


import jInv.ForwardShare.getSensMatVec
function getSensMatVec(x::Vector,p::Vector,pFor::Union{UnconstrainedOEDParam,EqconstrainedOEDParam})
	M   = getM(p,pFor.Design,false)
	D   = getD(p,pFor.Design)

	(m,n,nrhs) = getProblemSize(pFor.Design)
	nc   = getNumberOfConstraints(pFor)
	mv   = zeros(n,nrhs)
	for k=1:nrhs
		# get derivatives of design
		dpMx   = getdpMx(p,pFor.Xr[:,k],pFor.Design)
		dpMTMx = getdpMTx(p,M*pFor.Xr[:,k],pFor.Design)
		dpMTx  = getdpMTx(p,D[:,k],pFor.Design)
	 	dpD    = getdpD(p,k,pFor.Design)

		dwHx = M'*(dpMx*x) + dpMTMx*x
	    dwb  = M'*(dpD*x)  + dpMTx*x
		mv[:,k] = (pFor.Ainv\[dwb-dwHx;zeros(nc)])[1:n]
	end
    return vec(mv)  
end

import jInv.ForwardShare.getSensTMatVec
function getSensTMatVec(x::Vector,p::Vector,pFor::Union{UnconstrainedOEDParam,EqconstrainedOEDParam})

	M    = getM(p,pFor.Design,false)
	D    = getD(p,pFor.Design)
    n    = size(M,2)
    nrhs = size(D,2)
	nc   = getNumberOfConstraints(pFor)

    Hinvx = (pFor.Ainv\[reshape(x,n,nrhs);zeros(nc,nrhs)])[1:n,:]
    
	mv = zeros(getNumberOfDesignParameters(pFor.Design))
	for k=1:nrhs
		# get derivatives of design
		dpMx   = getdpMx(p,pFor.Xr[:,k],pFor.Design)
		dpMTMx = getdpMTx(p,M*pFor.Xr[:,k],pFor.Design)
		dpMTx  = getdpMTx(p,D[:,k],pFor.Design)
	 	dpD    = getdpD(p,k,pFor.Design)

    	dwHx = dpMx'*(M*Hinvx[:,k]) + dpMTMx'*Hinvx[:,k]
		dwb  = dpD'*(M*Hinvx[:,k]) + dpMTx'*Hinvx[:,k]
		mv  += dwb-dwHx
	end
    return mv
end