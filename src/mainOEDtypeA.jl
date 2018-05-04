using MAT
using jInv.Mesh
using jInv.LinearSolvers
using jInv.ForwardShare
using jInv.InverseSolve
using jInv.Utils
using ArgParse
using OED

function parse_commandline()
	s = ArgParseSettings()
	
	@add_arg_table s begin
	
	"--matfile"
		help = "matlab file storing problem description"
	"--exName"
		help = "name of example"
	"--model"
		help = "model function. Choose: 1-> parallelTomo, 2->identity"
		arg_type = Int
		default = 1
	"--IPmode"
		help = "inverse problem model. Choose: 0->unconstrained, 1->eqconstrained, 2->nonneq, 3->box "
		arg_type = Int
		default = 0
	"--alphaMin"
		help = "minimum regularization parameter for sparsity"
		default = "1e-10"
	"--alphaMax"
		help = "maximum regularization parameter for sparsity"
		default = "1e-1"
	"--nAlpha"
		help = "number of alphas"
		arg_type = Int
		default = 5
	"--maxIter"
		help = "maximum number of iteration"
		arg_type = Int
		default = 20
	"--noiseLevel"
		help = "noiseLevel"
		default = "1e-3"
	"--nex"
		help = "Number of examples. If -1 all examples will be used. Otherwise nEx examples will be randomly choosen"
		arg_type=Int
		default =-1
	end
	return parse_args(s)
end

function main()
	parsed_args = parse_commandline()
	
	matfile = matread(parsed_args["matfile"])
	A       = matfile["A"]
	
	if haskey(matfile,"L")
		L = matfile["L"]
		if size(L,2) != size(A,2)
			error("user supplied L has wrong size")
		end
	else
		L = speye(size(A,2))
	end
	Xtrue   = matfile["x_true"]
	B       = A*Xtrue
	theta   = matfile["theta"]
	nt      = length(theta)
	n       = Int(matfile["n"])
	nex     = parsed_args["nex"]
	
	Ce = matfile["Ce"]
	ce = Ce*Xtrue
	Ci = matfile["Ci"]
	ci = matfile["ci"]

	if nex!==-1
idx = 1:nex
		# idx = randperm(size(Xtrue,2))[1:nex]
		Xtrue = Xtrue[:,idx]
		B     = B[:,idx]
		ce    = ce[:,idx]
		ci    = ci[:,idx]
	end
		
	exName  = parsed_args["exName"]
	
	if parsed_args["model"] == 1 # parallel projections
		l      = round(Int64,size(A,1)/nt)
		modFun = v->parallelTomoMod(v,l)
		m0     = ones(nt)
		model  = :parallel
	elseif parsed_args["model"] == 2 # identity
		modFun = identityMod
		m0     = ones(size(A,1))
		model  = :identity
	end
	
	ipDict  = Dict(0=>:unconstr, 1=>:eqconstr, 2=>:nonnegconstr, 3=>:boxconstr)
	IPmode  = ipDict[parsed_args["IPmode"]]
	
	maxIter = parsed_args["maxIter"]	
	# add noise
	noiseLevel = parse(parsed_args["noiseLevel"])
	N  = randn(size(B));
	N /= norm(vec(N));
	N = noiseLevel*norm(vec(B))*N
	B += N	
	
	# divide problems among available workers
	ids = mod(0:size(Xtrue,2)-1,nworkers())+1
	OEDparam = Array{OEDTypeAParam}(nworkers())
	for k=1:nworkers()
		OEDparam[k] = OEDTypeAParam(A,B[:,ids.==k],zeros(0,0))
	end
	
	alphaMin = parse(parsed_args["alphaMin"])
	alphaMax = parse(parsed_args["alphaMax"])
	nAlpha   = parsed_args["nAlpha"]
	Alphas   = logspace(log10(alphaMin),log10(alphaMax),nAlpha)
	println("Alphas = $Alphas")
	
	runOEDTypeA(OEDparam,Xtrue,L,Alphas,modFun,IPmode,m0,Ce=Ce,Ci=Ci,ce=ce,ci=ci,maxIter=maxIter,exName=@sprintf "%s-%s-noise-%1.2e" exName model noiseLevel)	

end

main()

