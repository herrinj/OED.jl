export l1Reg

"""
l1Reg(m,mref,M)

l1-regularizer: R = sum(abs(m-mref))
"""
function l1Reg(m::Vector{Float64},mref::Vector{Float64},M::AbstractMesh)
	dm   = m .- mref
    dR   = ones(size(m))
    Rc   = sum(dm)
    return Rc,dR,1e-4*speye(length(m))
end
