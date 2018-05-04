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
	"--lambda"
		help = "regularization parameter for imaging problem"
		default = "1e-1"
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
	"--l"
		help = "Number of projection angles (design parameters)"
		arg_type=Int
		default=2
	"--OEDtype"
		help = "Type of OED problem: 2 -> OEDTypeB (default), 3->OEDTypeC"
		arg_type=Int
		default=2
	"--r"
		help = "Number of projections per angle, default: -1 --> pick all"
		arg_type=Int
		default=-1
	"--pad"
		help = "Padding of domain for imaging matrix"
		arg_type=Int
		default=4
	end
	return parse_args(s)
end

function main()
	parsed_args = parse_commandline()
	
	matfile = matread(parsed_args["matfile"])
	A       = matfile["A"]
	lambda  = parse(parsed_args["lambda"])
	if haskey(matfile,"L")
            L   = matfile["L"]
            if size(L,2)!=size(A,2); error("user supplied L has wrong size"); end
        else
	    L       = sqrt(lambda)*speye(size(A,2))
        end
	Xtrue   = matfile["x_true"]
	r       = parsed_args["r"]
	l       = parsed_args["l"]
	n       = Int(matfile["n"])
	nex     = parsed_args["nex"]
	OEDtype = parsed_args["OEDtype"]	
	
	Ce = matfile["Ce"]
	ce = Ce*Xtrue
	Ci = matfile["Ci"]
	ci = matfile["ci"]
	
	if nex!==-1
		idx = randperm(size(Xtrue,2))[1:nex]
		Xtrue = Xtrue[:,idx]
		ci   = ci[:,idx]
		ce   = ce[:,idx]
	end
	
	Mimg  = getRegularMesh([-1 1 -1 1.],[n,n])
	# get padded domain
	pad   = parsed_args["pad"]
	domain = copy(Mimg.domain);
	domain[1:2:end] -= pad*Mimg.h
	domain[2:2:end] += pad*Mimg.h
	Mgrid = getRegularMesh(domain,[n+2*pad,n+2*pad])
	
	if r==-1
		r = n+2*pad
	end
	
	print("Mimg.h=$(Mimg.h)\tvs.\tMgrid.h=$(Mgrid.h)")
	m0 = collect(linspace(1,179,l+2)[2:end-1])
	if OEDtype==2
		# divide problems among available workers
		ids = mod(0:size(Xtrue,2)-1,nworkers())+1
		OEDparam = Array{OEDTypeBParam}(nworkers())
		for k=1:nworkers()
			OEDparam[k] = getOEDTypeBParam(r,l,Xtrue[:,ids.==k],Mimg,Mgrid)
		end
	elseif OEDtype==3
		# divide problems among available workers
		ids = mod(0:size(Xtrue,2)-1,nworkers())+1
		OEDparam = Array{OEDTypeCParam}(nworkers())
		for k=1:nworkers()
			OEDparam[k] = getOEDTypeCParam(r,l,Xtrue[:,ids.==k],Mimg,Mgrid,L=copy(L))
		end
		m0 = [m0; lambda]
		L  = 0*L
	end
	display(OEDparam[1])
		
	exName  = parsed_args["exName"]

	ipDict  = Dict(0=>:unconstr, 1=>:eqconstr, 2=>:nonnegconstr, 3=>:boxconstr)
	IPmode  = ipDict[parsed_args["IPmode"]]
	
	maxIter = parsed_args["maxIter"]	
	m0 = [ 1.0 89.0; 10 100; 45. 135.; 44.9 45.1; 90. 90]'
	
	for k=1:size(m0,2)
		mOpt = runOEDTypeB(OEDparam,Xtrue,L,IPmode,vec(m0[:,k]),Ce=Ce,Ci=Ci,ce=ce,ci=ci,maxIter=maxIter,exName=@sprintf "%s-k-%d" exName k )	
		println("m0  : $m0")
		println("mOpt: $mOpt")
	end
end

main()

