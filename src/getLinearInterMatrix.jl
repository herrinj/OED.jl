export getLinearInterMatrix

function getLinearInterMatrix(M::RegularMesh,x::AbstractArray{Float64})
	return getLinearInterMatrix(M.domain,M.n,x)
end

"""
function A = getLinearInterpolationMatrix(domain,m,x)

builds interpolation matrix A, s.t. A(x)*I(:) = linearInter(I,domain,x)


Input:
 domain - description of computational domain
 m     - discretization size of image
 x     - interpolation points

Output:
 A     - linear interpolation matrix

see also linearInter

"""
function getLinearInterMatrix(domain::Array{Float64},n::Array{Int},xin::AbstractArray{Float64})

	h   = vec(domain[2:2:end]-domain[1:2:end])./vec(n)
	dim = round(Int,length(domain)/2)
	np  = round(Int,length(xin)/dim)
	x   = reshape(copy(xin),np,dim)

	# Convert x and y to the coordinate system 1:m, 1:n
	for i=1:dim
	    x[:,i] = (x[:,i]-domain[2*i-1])/h[i] + .5
	end

	Valid(j) = (0.<x[:,j]) .& (x[:,j].<(n[j]+1))

	if dim==1
		valid = find( Valid(1) )
	elseif dim===2
		valid = find( Valid(1) .& Valid(2) )
	elseif dim==3
		valid = find( Valid(1) .& Valid(2) .& Valid(3) )
	end
	if isempty(valid)
	    return spzeros(prod(n),np)
	end;
	P     = round.(Int,floor.(x)); x = x-P;    # split x into integer/remainder
	p(j)  = P[valid,j]
	xi(j) = x[valid,j]

	if dim==1
		error("nyi")
    elseif dim==2
        A    = getMatrix2D(valid, p(1)  , p(2)  , np, n, (1-xi(1)).*(1-xi(2)));
        A   += getMatrix2D(valid, p(1)+1, p(2)  , np, n, (xi(1))  .*(1-xi(2)));
        A   += getMatrix2D(valid, p(1)  , p(2)+1, np, n, (1-xi(1)).*(xi(2)));
        A   += getMatrix2D(valid, p(1)+1, p(2)+1, np, n, (xi(1))  .*(xi(2)));

    elseif dim== 3
		error("nyi")
    end
return A
end

function getMatrix2D(I::Vector{Int},p1::Vector{Int},p2::Vector{Int},np::Int,n::Array{Int},weight::Vector{Float64})
    valid =  find((p1.>=1) .& (p1.<=n[1]) .& (p2.>=1) .& (p2.<=n[2]));
    return sparse(I[valid],p1[valid]+(p2[valid]-1)*n[1],weight[valid],np,prod(n));
 end
