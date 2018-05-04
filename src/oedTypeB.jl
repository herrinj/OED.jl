export rotation2D, OEDTypeBParam, getOEDTypeBParam
export getTransmissionMatrix

"""
type OEDTypeBParam <: AbstractOEDParam

Fields:

	r          - number of rays
	l          - number of projections (design parameters)
	Ftrue      - true models
	Mimg       - image domain
	Mgrid      - rotating grid
	c          - center of rotation
	noiseLevel -
	Q          - temp storage for rotation2D
	xc         - cell-centers of Mgrid
	T          - transmission matrix
	M          - forward operator
	thetaOld   - previous theta (to check if M must be rebuilt)
"""
type OEDTypeBParam <: AbstractOEDParam
	r::Int
	l::Int
	Ftrue::Array{Float64,2}
	Mimg::RegularMesh
	Mgrid::RegularMesh
	c::Vector{Float64}
	noiseLevel::Float64
	xc::Array{Float64}
	T::AbstractArray{Float64}
	Q::Array{Float64}
	M::Array{Float64}
	thetaOld::Vector{Float64}
end
function Base.display(param::OEDTypeBParam)
	println("--- OED.OEDTypeBParam ---")
	println("number of rays   = $(param.r)")
	println("number of proj.  = $(param.l)")
	println("number of images = $(size(param.Ftrue,2))")
	println("image size       = $(param.Mimg.n)")
	println("image domain     = $(param.Mimg.domain)")
	println("grid size        = $(param.Mgrid.n)")
	println("gird domain      = $(param.Mgrid.domain)")
	println("noise level      = $(param.noiseLevel)")
	println("rotation center  = $(param.c)")

	print("temp variables   = ")
	(!isempty(param.xc)) && print("xc $(size(param.xc)), ")
	(!isempty(param.T)) && print("T $(size(param.T)), ")
	(!isempty(param.Q)) && print("Q $(size(param.Q)).")
	(!isempty(param.M)) && print("M $(size(param.M)).")
	print("\n")
end

export getOEDTypeBParam
"""
function getOEDTypeBParam

constructs OEDTypeBParam

Required Input:

	r          - number of rays
	l          - number of projections
	Ftrue      - true models
	Mimg       - image domain
	Mgrid      - rotating grid

Optional Inputs:

	c          - center of rotation
	noiseLevel -

"""
function getOEDTypeBParam(r::Int,l::Int,Ftrue,Mimg::RegularMesh,Mgrid::RegularMesh;
	c=(Mgrid.domain[2:2:end]+Mgrid.domain[1:2:end])/2, noiseLevel=1e-4)
	return OEDTypeBParam(r,l,Ftrue,Mimg,Mgrid,c,noiseLevel,[],[],[],[],[])
end

getProblemSize(param::OEDTypeBParam) = (param.r*param.l,param.Mimg.nc,size(param.Ftrue,2))
getNumberOfDesignParameters(pFor::OEDTypeBParam) = pFor.l

function getGrid(pFor::OEDTypeBParam)
	if isempty(pFor.xc)
		pFor.xc = getCellCenteredGrid(pFor.Mgrid)
	end
	return pFor.xc
end


function rotation2D(w::Float64,pFor::OEDTypeBParam;doDerivative::Bool=false)
	if isempty(pFor.Q) || size(pFor.Q,1) != length(pFor.xc)
		xc     = getGrid(pFor)
		pFor.Q = getQ(xc)
	end
	return rotation2D(w,pFor.xc,Q=pFor.Q,doDerivative=doDerivative,c=pFor.c)
end

function getTransmissionMatrix(pFor)
	if isempty(pFor.T)
		Mgrid = pFor.Mgrid
		r     = pFor.r;

		start = floor(Int,(Mgrid.n[1]-r)/2);
		v = ones(1,Mgrid.n[1]);
		T = zeros(r,prod(Mgrid.n));
		T[:, start*Mgrid.n[1]+1:start*Mgrid.n[1]+r*Mgrid.n[1]] = kron(eye(r),v);
		pFor.T = T;
	end
	return pFor.T
end


function getM(theta::Vector{Float64}, pFor::OEDTypeBParam,doClear::Bool=true)
	# This function gets the projection matrix by rotation and summation
	if doClear || isempty(pFor.thetaOld) || norm(theta-pFor.thetaOld,Inf)/norm(theta,Inf) > 1e-14
		Mgrid = pFor.Mgrid
		Mimg  = pFor.Mimg
		r     = pFor.r;
		l     = pFor.l

		# Get projection angles in radians
		psi = theta*(pi/180);
		if length(psi)!=l; error("getM theta must be length $l."); end;
		xc = getCellCenteredGrid(Mgrid);

		M     = zeros(r*l,Mimg.nc);

		T = getTransmissionMatrix(pFor)
		for i = 1:l
		  yc  = rotation2D(psi[i],pFor)
		  R   = getLinearInterMatrix(Mimg,yc)

		  M[(i-1)*r+1:i*r,:] = T*R
		end
		pFor.thetaOld = theta; pFor.M= M
	end
	return pFor.M
end
function getdpMx(theta::Vector{Float64},x::Vector{Float64},pFor::OEDTypeBParam)
#  This function gets the Jacobian matrix times a vector f

	Mgrid = pFor.Mgrid
	Mimg  = pFor.Mimg

	# Get projection angles in radians and create grid
	psi = theta*(pi/180);
	xc  = getCellCenteredGrid(Mgrid);

	# Build transmission matrix T
	T  = getTransmissionMatrix(pFor)

	# Build Jacobian matrix
	Jac = Array{Array}(pFor.l)
	for i = 1:pFor.l
	  yc,yp = rotation2D(psi[i],pFor,doDerivative=true)
	  Rc,dR = linearInter(reshape(x,tuple(Mimg.n...)),Mimg.domain,yc,doDerivative=true)
	  Jac[i] = (pi/180)*T*(dR*yp)
	end
	dpM = SparseMatrixCSC(sparse(Jac[1]))
	for i=2:pFor.l
		dpM = blkdiag(dpM,SparseMatrixCSC(sparse(Jac[i])))
	end
	return dpM
end

function getdpMTx(theta::Vector{Float64}, x::Vector{Float64}, pFor::OEDTypeBParam)
# This function gets the Jacobian matrix (with respect to the transpose)
#   times a vector x

	Mgrid = pFor.Mgrid
	Mimg  = pFor.Mimg

	# Get projection angles in radians
	psi = theta*(pi/180)
	xc  = getGrid(pFor)

	T   = getTransmissionMatrix(pFor)
	TTx = T'*reshape(x,pFor.r,pFor.l);
	Jac = zeros(Mimg.nc,0);

	for i = 1:pFor.l
	  yc,yp  = rotation2D(psi[i],pFor,doDerivative=true);
	  AT,dAT = getLinearInterMatrixTranspose(Mimg.domain,Mimg.n,yc,Img=TTx[:,i],doDerivative=true)
	  Jac = [Jac (pi/180)*dAT*yp];
	end
	return Jac
end


function getD(theta::Vector{Float64},param::OEDTypeBParam)
	D  = getM(theta,param) *param.Ftrue
	D += param.noiseLevel*randn(size(D))
	return D
end

function getdpD(theta::Vector{Float64},k::Int,param::OEDTypeBParam)
	return getdpMx(theta,param.Ftrue[:,k],param)
end
