export UnconstrainedOEDParam, getNumberOfConstraints

import jInv.ForwardShare.ForwardProbType

"""
type UnconstrainedOEDParam <: ForwardProbType

Parameter for OED with unconstrained lower level problem

	f(p) = argmin_f 0.5*| M(p)*f - d(p)|^2 + 0.5*|L*f|^2 

STATUS: this is under development

Fields:

	Design  - design type (OEDTypeAParam)
	L       - regularization Operator
	Xr      - field for reconstructions
	Ainv    - field for linear solvers
"""
type UnconstrainedOEDParam <: ForwardProbType
   Design::AbstractOEDParam   
   L 
   Xr
   Ainv
end

getNumberOfConstraints(pFor::UnconstrainedOEDParam) = 0

function Base.display(pFor::UnconstrainedOEDParam)
	(m,n,nrhs) = getProblemSize(pFor.Design)
	println("--- OED.UnconstrainedOEDParam <: ForwardProbType ---")
	println("Type of Design:                    $(typeof(pFor.Design))")
	println("Size of linear system:             $m x $n")
	println("Number of right hand sides:        $nrhs")
	println("Size of regularization operator:   $(size(pFor.L))")
end

import jInv.ForwardShare.getData
function getData(p::Vector{Float64},pFor::UnconstrainedOEDParam,doClear::Bool=false)
	# get design
	M  = getM(p,pFor.Design, true)
	D  = getD(p,pFor.Design)
	L  = pFor.L

	# build and factorize normal equations
	H   = M'*M + L'*L
	pFor.Ainv = cholfact(H)
	
	Xr = pFor.Ainv\(M'*D)

    pFor.Xr = Xr
    return Xr,pFor
end

