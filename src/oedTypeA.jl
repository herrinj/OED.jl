export OEDTypeAParam, OEDParam, getProblemSize


"""
type OED.OEDTypeAParam

data structure for OED Type A problem where design is given by

M(p) = sdiag(p)*A   and D(p) = sdiag(p)*D	

Fields:
	A - forward operator
	D - data	
"""
type OEDTypeAParam <: AbstractOEDParam
	A::AbstractArray{Float64} 
	D::Array{Float64,2}
	M::AbstractArray{Float64,2}
end

function Base.display(param::OEDTypeAParam)
	println("--- OED.OEDTypeAParam ---")
	println("size(A) = $(size(param.A))")
	println("size(D) = $(size(param.D))")
end


getProblemSize(param::OEDTypeAParam) = (size(param.A,1),size(param.A,2),size(param.D,2))
getNumberOfDesignParameters(pFor::OEDTypeAParam) = size(pFor.A,1)

function getM(p::Vector{Float64},param::OEDTypeAParam,doClear::Bool=true)
	if doClear || isempty(param.M) 
	  param.M =  sdiag(p)*param.A
	end
	return param.M
end
getdpMx(p::Vector{Float64},x::Vector{Float64},param::OEDTypeAParam)  = opDiagonal(param.A*x)
getdpMTx(p::Vector{Float64},x::Vector{Float64},param::OEDTypeAParam) = (LinearOperator(param.A)'*opDiagonal(x))
getD(p::Vector{Float64},param::OEDTypeAParam)       = sdiag(p)*param.D
getdpD(p::Vector{Float64},k::Int,param::OEDTypeAParam)     = opDiagonal(param.D[:,k]); # sdiag(param.D[:,k])



