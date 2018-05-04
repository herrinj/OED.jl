export AbstractOEDParam, getProblemSize, getM, getdpMx, getdpMTx, getD, getdpD, getNumberOfDesignParameters

"""
	abstract OED.AbstractOEDParam

	generic type for OED problems. For an example see OEDTypeAParam
"""
abstract type AbstractOEDParam end


"""
m,n,nrhs = getProblemSize(pFor::AbstractOEDParam)

"""
getProblemSize(pFor::AbstractOEDParam) = error("nyi for pFor of type $(typeof(pFor))")

"""
function OED.getM

constructs matrix M(p) for OED problems

Input:

p     - design parameters
pFor  <: AbstractOEDParam (describing the construction)

Output:

M     - forward operator
"""
getM(p::Vector{Float64},pFor::AbstractOEDParam)       = error("nyi for pFor of type $(typeof(pFor))")

"""
function OED.getdpMx

constructs gradient matrix  d_p (M(p)*x)

Input:

p     - design parameters
x     - projection
pFor  <: AbstractOEDParam (describing the construction)

Output:

M     - forward operator
"""
getdpMx(p::Vector{Float64},x::Vector{Float64},pFor::AbstractOEDParam)  = error("nyi for pFor of type $(typeof(pFor))")

"""
function OED.getdpMTx

constructs gradient matrix  d_p (M(p)'*x)

Input:

p     - design parameters
x     - projection
pFor  <: AbstractOEDParam (describing the construction)

Output:

M     - forward operator
"""
getdpMTx(p::Vector{Float64},x::Vector{Float64},pFor::AbstractOEDParam) = error("nyi for pFor of type $(typeof(pFor))")

"""
function OED.getD

get data for current design

Input:

p     - design parameters
pFor  <: AbstractOEDParam (describing the construction)

Output:

M     - forward operator
"""
getD(p::Vector{Float64},pFor::AbstractOEDParam)       = error("nyi for pFor of type $(typeof(pFor))")

"""
function OED.getdpD

gets gradient matrix d_p (D(p)[:,k])

Input:

p     - design parameters
k     - index of data
pFor  <: AbstractOEDParam (describing the construction)

Output:

M     - forward operator
"""
getdpD(p::Vector{Float64},k::Int,pFor::AbstractOEDParam)     = error("nyi for pFor of type $(typeof(pFor))")

"""
function OED.getNumberOfDesignParameters(pFor)
"""
getNumberOfDesignParameters(pFor::AbstractOEDParam)     = error("nyi for pFor of type $(typeof(pFor))")
