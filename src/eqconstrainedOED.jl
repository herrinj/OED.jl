export EqconstrainedOEDParam



import jInv.ForwardShare.ForwardProbType

"""
type EqconstrainedOEDParam <: ForwardProbType

Parameter for OED with equality constrained lower-level problem

	f(p) = argmin_f 0.5*| M(p)*f - d(p)|^2 + 0.5*|L*f|^2 
	       subject to Ce*f = ce

STATUS: this is under development

Fields:

	Design  - design type (OEDTypeAParam)
	L       - Regularization Operator
	Ce      - equality constraint matrix
	ce      - equality constraints
	Xr      - field for reconstructions
	Ainv    - field for linear solvers
"""
type EqconstrainedOEDParam <: ForwardProbType
   Design::AbstractOEDParam   
   L::AbstractArray{Float64}
   Ce::AbstractArray{Float64,2}
   ce::AbstractArray{Float64}
   Xr::AbstractArray{Float64,2}
   Ainv
end

function getKKT(H::SparseMatrixCSC{Float64}, Ce::AbstractArray)
	ne = size(Ce,1)
	if ne==0
		return H
	else
		return [H Ce';Ce spzeros(ne,ne)]
	end
end

function getKKT(H::Array{Float64}, Ce::AbstractArray)
	ne = size(Ce,1)
	if ne==0
		return H
	else
		return [H Ce';Ce 0]
	end
end

getNumberOfConstraints(pFor::EqconstrainedOEDParam) = size(pFor.Ce,1)

function Base.display(pFor::EqconstrainedOEDParam)
	(m,n,nrhs) = getProblemSize(pFor.Design)
	println("--- OED.EqconstrainedOEDParam <: ForwardProbType ---")
	println("Type of Design:                    $(typeof(pFor.Design))")
	println("Size of regularization operator:   $(size(pFor.L))")
	println("Number of constraints:             $(size(pFor.Ce,1))")
end

import jInv.ForwardShare.getData
function getData(p::Vector{Float64},pFor::EqconstrainedOEDParam,doClear=false)
    
	# get design
	M  = getM(p,pFor.Design,true)
	D  = getD(p,pFor.Design)
	Ce = pFor.Ce
	ce = pFor.ce
	L  = pFor.L
	n  = size(M,2)
	ne = size(Ce,1)

	# build and factorize KKT system
	H   = M'*M + L'*L
	KKT = getKKT(H,Ce)
    pFor.Ainv = lufact(KKT)
	
    xl  = pFor.Ainv\[M'*D;ce]
    Xr  = xl[1:n,:]
    
	pFor.Xr = Xr
    return Xr,pFor
end



