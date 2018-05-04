export IneqconstrainedOEDParam, getIneqconstrainedOEDParam


import jInv.ForwardShare.ForwardProbType

"""
type IneqconstrainedOEDParam <: ForwardProbType

Parameter for OED with equality constrained lower-level problem

	f(p) = argmin_f 0.5*| M(p)*f - d(p)|^2 + 0.5*|L*f|^2
	       subject to Ce*f = ce and Ci*f >= ci

STATUS: this is under development

Fields:

	Design         - design type (OEDTypeAParam)
	L              - Regularization Operator
	Ce             - equality constraint matrix
	ce             - equality constraints
	Ci             - inequality constraints
	ci             - equality constraints
	Xr             - field for reconstructions
	Yr             - field for slacks
	Le             - field for Lagrange multipliers of equality constraints
	Li             - field for Lagrange multipliers of inequality constraints
	Ainv           - field for factorizations
	maxIterIPQP    - maxIter for interior point method
	outIPQP        - flag for IP output (0:no output, 1:final status, 2:complete history)
	ftolIPQP       - tolerance on optimality condition
	mutolIPQP      - tolerance on centrality
	solveKKT       - KKT solver (options: solveFullKKT, solveAugKKT, solveNormalEqKKT,...)
"""
type IneqconstrainedOEDParam <: ForwardProbType
   Design::AbstractOEDParam
   L::AbstractArray{Float64}
   Ce::AbstractArray{Float64}
   ce::Array{Float64}
   Ci::AbstractArray{Float64}
   ci::Array{Float64}
   Xr::Array{Float64}
   Yr::Array{Float64}
   Le::Array{Float64}
   Li::Array{Float64}
   Ainv
   maxIterIPQP::Int64
   outIPQP::Int64
   ftolIPQP::Real
   mutolIPQP::Real
   solveKKT::Function
end

function Base.display(pFor::IneqconstrainedOEDParam)
	(m,n,nrhs) = getProblemSize(pFor.Design)

	println("--- OED.IneqconstrainedOEDParam <: ForwardProbType --- ")
	println("Design:                            $(typeof(pFor.Design))")
	println("Size of linear system              $m x $n")
	println("Number of right hand sides         $nrhs")
	println("Size of regularization operator:   $(size(pFor.L))")
	println("Number of equality constraints:    $(size(pFor.Ce,1))")
	println("Number of inequality constraints:  $(size(pFor.Ci,1))")
	println("max. iter of IPQP solver           $(pFor.maxIterIPQP)")
	println("verbosity of IPQP solver           $(pFor.outIPQP)")
	println("ftol of IPQP solver                $(pFor.ftolIPQP)")
	println("centrality tol of IPQP solver      $(pFor.mutolIPQP)")
	println("linear solver in IPQP              $(pFor.solveKKT)")
end

getNumberOfConstraints(pFor::IneqconstrainedOEDParam) = size(pFor.Ce,1)+size(pFor.Ci,1)

"""
function getIneqconstrainedOEDParam(Design,L,Ce,ce,Ci,ci,kwargs...)

	constructs and returns an IneqconstrainedOEDParam

	Required Inputs:

		Design         - design parameter
		L              - regularization operator
		Ce,ce          - equality constraints
		Ci,ci          - inequality constraints

	Optional keyword arguments:

		Xr             - starting guess for inverse solutions (default=[])
		Yr             - starting guess for slacks (default=[])
		Le             - starting guess for Lagrange multipliers assoc. with equality
		Li             - starting guess for Lagrange multipliers assoc. with inequality
		maxIterIPQP    - maximum number of iteration for ipqp(...), default=20
		outIPQP        - verbosity of ipqp(...), default=0
		ftolIPQP       - ftol of ipqp(...), default=1e-13
		mutolIPQP      - mutol of ipqp(...), default=1e-13
		solveKKT       - linear solver used in ipqp(...), default=OED.solveNormalEqKKT
"""
function getIneqconstrainedOEDParam(Design::AbstractOEDParam,L::AbstractArray{Float64},Ce,ce,Ci,ci;Xr=[],Yr=[],Le=[],Li=[], maxIterIPQP=20,outIPQP=0,ftolIPQP=1e-13,mutolIPQP=1e-13,solveKKT=OED.solveNormalEqKKT)

	(m,n,nrhs) = getProblemSize(Design)

	Ce = (isempty(Ce)) ? zeros(0,n) : Ce
	ce = (isempty(ce)) ? zeros(0,nrhs) : ce
	if size(Ce,1) != size(ce,1);
		error("size(Ce,1) = $(size(Ce,1)) != $(size(ce,1)) size(ce,1) ")
    end

	return IneqconstrainedOEDParam(Design,L,Ce,ce,Ci,ci,Xr,Yr,Le,Li,[],maxIterIPQP,outIPQP,ftolIPQP,mutolIPQP,solveKKT)
end

import jInv.ForwardShare.getData
function getData(p::Vector{Float64},pFor::IneqconstrainedOEDParam,doClear=false)
    (m,n,nrhs) = getProblemSize(pFor.Design)

	# get design
	M  = getM(p,pFor.Design, true)
	D  = getD(p,pFor.Design)
	Ce = pFor.Ce
	ce = pFor.ce
	Ci = pFor.Ci
	ci = pFor.ci
	L  = pFor.L
	n  = size(M,2)
	ne = size(Ce,1)
	ni = size(Ci,1)

	# build QP: min_f 0.5*f'*Q*f + C'*f
	Q   = M'*M + L'*L
	C   = -M'*D

	# get initial values
	Xr = (isempty(pFor.Xr)) ? -(Q\C)                   : pFor.Xr
	Yr = (isempty(pFor.Yr)) ? Ci*Xr - ci               : pFor.Yr
	Le = (isempty(pFor.Le)) ? zeros(ne,nrhs)           : pFor.Le
	Li = (isempty(pFor.Li)) ? .1*ones(ni,nrhs)         : pFor.Li
	# make sure slack variable and lagrange multiplier are sufficiently large
	Yr = max.(.1,Yr)
	Li = max.(.1,Li)

    pFor.Xr = Xr; pFor.Yr = Yr; pFor.Le=Le; pFor.Li=Li

	Ainv = Array{Any}(nrhs)
	for k=1:nrhs
		xr,yr,lek,lik,his,Ainv[k] = ipqp(Q,vec(C[:,k]),Ce,vec(ce[:,k]), Ci,vec(ci[:,k]),			 				Xr[:,k],Yr[:,k],Le[:,k],Li[:,k],out=pFor.outIPQP,maxIter=pFor.maxIterIPQP,ftol=pFor.ftolIPQP)
		pFor.Xr[:,k]  = xr
		pFor.Yr[:,k]  = yr
		pFor.Le[:,k]  = lek
		pFor.Li[:,k]  = lik
	end
	pFor.Ainv = Ainv
    return copy(pFor.Xr),pFor
end

import jInv.ForwardShare.getSensMatVec
function getSensMatVec(x::Vector{Float64},p::Vector{Float64},pFor::IneqconstrainedOEDParam)
	(m,n,nrhs) = getProblemSize(pFor.Design)

	M   = getM(p,pFor.Design,false)
	D   = getD(p,pFor.Design)

	Ce  = pFor.Ce
    ce  = pFor.ce
	Ci  = pFor.Ci
    ci  = pFor.ci
    L   = pFor.L
	Xr  = pFor.Xr
	Yr  = pFor.Yr
	Le  = pFor.Le
	Li  = pFor.Li
	mi   = size(Ci,1)
    me   = size(Ce,1)

	# build QP: min_f 0.5*f'*Q*f + C'*f
	rebuildQ = false
	for k=1:nrhs
		if  (size(pFor.Ainv[k],1)!=n)
			println("Sens needs to rebuild")
			rebuildQ = true
			break
		end
	end
	Q   = rebuildQ ?  M'*M + L'*L : spzeros(n,n)
	C   = -M'*D

	mv   = zeros(n,nrhs)
	for k=1:nrhs
		# get derivatives of design
		dpMx  = getdpMx(p,pFor.Xr[:,k],pFor.Design)
		dpMTMx = getdpMTx(p,M*pFor.Xr[:,k],pFor.Design)
		dpMTx = getdpMTx(p,D[:,k],pFor.Design)
	 	dpD   = getdpD(p,k,pFor.Design)

		# derivative of right hand side
		dpHx = M'*(dpMx*x) + dpMTMx*x
		dpb  = M'*(dpD*x)  + dpMTx*x

		tmp,     = pFor.solveKKT(Q,Ce,Ci,Yr[:,k],Li[:,k],dpb-dpHx,zeros(me),zeros(mi),zeros(mi),KKTinv=pFor.Ainv[k],doClear=false)
		(dx,dy,dle,dli) = OED.splitVecIPQP(tmp,n,me,mi)
		mv[:,k]  = dx
	end
    return vec(-mv)
end

import jInv.ForwardShare.getSensTMatVec
function getSensTMatVec(x::Vector{Float64},p::Vector{Float64},pFor::IneqconstrainedOEDParam)
	(m,n,nrhs) = getProblemSize(pFor.Design)

	M    = getM(p,pFor.Design,false)
	D    = getD(p,pFor.Design)

	Ce  = pFor.Ce
    ce  = pFor.ce
	Ci  = pFor.Ci
    ci  = pFor.ci
    L   = pFor.L
	Xr  = pFor.Xr
	Yr  = pFor.Yr
	Le  = pFor.Le
	Li  = pFor.Li
	mi   = size(Ci,1)
    me   = size(Ce,1)

	# build QP: min_f 0.5*f'*Q*f + C'*f
	C   = -M'*D
	rebuildQ = false
	for k=1:nrhs
		if   (size(pFor.Ainv[k],1)!=n)
			println("Sens T needs to rebuild")
			rebuildQ = true
			break
		end
	end
	Q   = rebuildQ ?  M'*M + L'*L : spzeros(n,n)

    x = reshape(x,n,nrhs)
	mv = zeros(getNumberOfDesignParameters(pFor.Design))
	for k=1:nrhs
    	tmp,      = pFor.solveKKT(Q,Ce,Ci,Yr[:,k],Li[:,k],-x[:,k],zeros(me),zeros(mi),zeros(mi),KKTinv=pFor.Ainv[k],doClear=false)
		(dx,dy,dle,dli) = OED.splitVecIPQP(tmp,n,me,mi)

		# get derivatives of design
		dpMx  = getdpMx(p,pFor.Xr[:,k],pFor.Design)
		dpMTMx = getdpMTx(p,M*pFor.Xr[:,k],pFor.Design)
		dpMTx = getdpMTx(p,D[:,k],pFor.Design)
	 	dpD   = getdpD(p,k,pFor.Design)

    	dwHx = dpMx'*(M*dx) + dpMTMx'*dx
		dwb  = dpD'*(M*dx) + dpMTx'*dx
		mv  += dwb-dwHx
	end
    return mv
end
