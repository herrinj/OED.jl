export rotation2D

"""
getQ(xc)

builds basis for fast computation of rigid2D
"""
function getQ(xc::Array{Float64})
	n  = round(Int,length(xc)/2)
	xc = reshape(xc,n,2)
	return  kron(speye(2), [xc ones(n)])
end

"""
rotation2D(w,x,doDerivative,Q)

returns rotated points x and derivatives.
"""
function rotation2D(w::Float64,x::Array{Float64};doDerivative::Bool=false,Q::AbstractArray{Float64,2}=getQ(x),c::Vector{Float64}=zeros(2))

R  = [ cos(w) -sin(w);sin(w) cos(w)]
g  = (eye(2)-R)*vec(c)
f  = [R[1,1];R[1,2];g[1];R[2,1];R[2,2];g[2]]
y  = Q*f

if doDerivative
	dR = [-sin(w[1]) -cos(w[1]);cos(w[1]) -sin(w[1])]
	dg = -dR*vec(c)
	df = [dR[1,1];dR[1,2];dg[1];dR[2,1];dR[2,2];dg[2]]
	dy = Q*df
	return y,dy
else
	return y
end
end
