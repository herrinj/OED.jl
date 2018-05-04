export getLinearInterMatrixTranspose

function findValid(js::Array{Int,2},n::Array{Int})
	nv  =  size(js,1)
	val = BitArray(nv)
	for k=1:nv
	 val[k] =  (js[k,1]>0) && (js[k,1] <=n[1]) && (js[k,2] > 0) && (js[k,2]<=n[2])
	end
	return val
end


"""
function [A,dATc] = getLinearInterMatrixTranspose(domain,m,x,T)

builds transpose of linear interpolation matrix A,
                           s.t. A(x)'*T(:) = linearInter(T,domain,x)

Input:
 domain    - description of computational domain
 m        - discretization size of image
 x        - interpolation points
 T        - interpolation coefficients
 varargin - optional parameters

Output:
 A     - linear interpolation matrix
 dAT   - derivative of (A'*T) w.r.t. x

"""
function getLinearInterMatrixTranspose(domain::Array{Float64},n::Array{Int},x::Array{Float64};doDerivative::Bool=false,Img::Vector{Float64}=zeros(0))

h   = vec(domain[2:2:end]-domain[1:2:end])./vec(n)
dim = round(Int,length(domain)/2)
np  = round(Int,length(x)/dim)
x   = reshape(copy(x),np,dim)

# Convert x and y to the coordinate system 1:m, 1:n
for i=1:dim
    x[:,i] = (x[:,i]-domain[2*i-1])/h[i] + .5
end

Valid(j) = trues(size((x[:,j])));#; (0<x(:,j) & x(:,j)<m(j)+1);      # determine indices of valid points

if dim==1
	valid = find( Valid(1) )
elseif dim===2
	valid = find( Valid(1) .& Valid(2) )
elseif dim==3
	valid = find( Valid(1) .& Valid(2) .& Valid(3) )
end

if isempty(valid)
    return spzeros(np,prod(n))
end;

P     = round.(Int,floor.(x)); x = x-P;    # split x into integer/remainder
p(j)  = P[valid,j]
xi(j) = x[:,j]

if dim==1
	error("nyi")
elseif dim==2
	ind  = valid

	# get indices of neighbouring points
	jk    = [p(1)   p(2)  ]
	j1k   = [p(1)+1 p(2)  ]
	jk1   = [p(1)   p(2)+1]
	j1k1  = [p(1)+1 p(2)+1]

	# compute weights for interpolation
	s1 = (1-xi(1)).*(1-xi(2))
	s2 = (xi(1)).*(1-xi(2))
	s3 = (1-xi(1)).*(xi(2))
	s4 =  (xi(1)).*(xi(2))
	ii = [ind  ind  ind  ind]
	jj = [jk ; j1k; jk1; j1k1]

	# points that are outside the domain
	# valid = (js::Array{Int,2}) -> find((js[:,1].>0)&(js[:,1].<=n[1])&(js[:,2].>0)&(js[:,2].<=n[2]))
	# valid = (js::Array{Int,2}) -> (js[:,1].>0)&(js[:,1].<=n[1])&(js[:,2].>0)&(js[:,2].<=n[2])
    valid = (ij::Array{Int,2}) -> findValid(ij,n)
    s2ind = (jj::Array{Int,2}) -> jj[:,1]+(jj[:,2]-1)*n[1]

        jkt = valid(jk); ijkt = ind[jkt]
		A1 = sparse(s2ind(jk[jkt,:])     ,ijkt   ,s1[ijkt]  ,prod(n),np);
        A2 = sparse(s2ind(j1k[valid(j1k),:])   ,ind[valid(j1k)]  ,s2[ind[valid(j1k)]] ,prod(n),np);
        A3 = sparse(s2ind(jk1[valid(jk1),:])   ,ind[valid(jk1)]  ,s3[ind[valid(jk1)]] ,prod(n),np);
        A4 = sparse(s2ind(j1k1[valid(j1k1),:]) ,ind[valid(j1k1)] ,s4[ind[valid(j1k1)]],prod(n),np);

        A =  A1 + A2 + A3 + A4;
        if doDerivative
			if isempty(Img); error("need to provide direction"); end

            dv1dxi = zeros(np); dv1deta = zeros(np)
            dv1dxi[ind] = -(1-xi(2)); dv1deta[ind] = -(1-xi(1))

            dv2dxi = zeros(np); dv2deta = zeros(np);
            dv2dxi[ind] = (1-xi(2)); dv2deta[ind] = -(xi(1));

            dv3dxi = zeros(np); dv3deta = zeros(np);
            dv3dxi[ind] = (-xi(2)); dv3deta[ind] = 1-xi(1);

            dv4dxi = zeros(np); dv4deta = zeros(np);
            dv4dxi[ind] = (xi(2)); dv4deta[ind] = xi(1);
            # built sdiags before creating function handle
            dv1dxi = sdiag(dv1dxi);
            dv2dxi = sdiag(dv2dxi);
            dv3dxi = sdiag(dv3dxi);
            dv4dxi = sdiag(dv4dxi);


			A1 = sparse(s2ind(jk[valid(jk),:])     ,ind[valid(jk)]   ,Img[ind[valid(jk)]]  ,prod(n),np);
	        A2 = sparse(s2ind(j1k[valid(j1k),:])   ,ind[valid(j1k)]  ,Img[ind[valid(j1k)]] ,prod(n),np);
	        A3 = sparse(s2ind(jk1[valid(jk1),:])   ,ind[valid(jk1)]  ,Img[ind[valid(jk1)]] ,prod(n),np);
	        A4 = sparse(s2ind(j1k1[valid(j1k1),:]) ,ind[valid(j1k1)] ,Img[ind[valid(j1k1)]],prod(n),np);

			Dxi =  A1*dv1dxi + A2*dv2dxi +  A3*dv3dxi + A4*dv4dxi

            # built sdiags before creating function handle
            dv1deta = sdiag(dv1deta);
            dv2deta = sdiag(dv2deta);
            dv3deta = sdiag(dv3deta);
            dv4deta = sdiag(dv4deta);
            Deta =  A1*dv1deta + A2*dv2deta + A3*dv3deta + A4*dv4deta
            dAT =  [Dxi/h[1]  Deta/h[2]];
		end
   elseif dim==3
	error("NYI")
	end
if !doDerivative
	return A
else
	return A, dAT;
end
end
