export runOEDTypeB

"""

function runOEDTypeA

Input:

	param
	Xtrue
	L
	n
	Alphas
	modFun
	IPmode
	m0
"""
function runOEDTypeB(OEDparam::Array,Xtrue,L,IPmode,m0;Ce=[],ce=[],Ci=[],ci=[],pcgMaxIter=10,maxIter=20,maxStep=1.0,exName="")

	ids = mod.(0:size(Xtrue,2)-1,nworkers())+1

	println("\n\n=== setup forward problem ===")
	if IPmode==:unconstr
		pFor = Array{UnconstrainedOEDParam}(nworkers())
		for j=1:nworkers()
			pFor[j] = UnconstrainedOEDParam(OEDparam[j],L,zeros(0,0),[])
		end
	elseif IPmode==:eqconstr
		if isempty(Ce) || isempty(ce)
			error("equality constraints must be provided. To run unconstrained OED, use IPmode=:unconstr.")
		end
		pFor = Array{EqconstrainedOEDParam}(nworkers())
		for j=1:nworkers()
			pFor[j] = EqconstrainedOEDParam(OEDparam[j],L,Ce,ce[:,ids.==j],zeros(0,0),[])
		end
	elseif IPmode==:nonnegconstr
		I = speye(size(Xtrue,1))
		pFor = Array{IneqconstrainedOEDParam}(nworkers())
		for j=1:nworkers()
				pFor[j] = getIneqconstrainedOEDParam(OEDparam[j],L,Ce,ce[:,ids.==j],I,0*OEDparam[j].Ftrue,
				          ftolIPQP=1e-3,mutolIPQP=1e-3,outIPQP=0)
		end

	elseif IPmode==:boxconstr
		if isempty(Ci) || isempty(ci)
			error("inequality constraints must be provided. To run equality constrained OED, use IPmode=:eqconstr.")
		end
		pFor = Array{IneqconstrainedOEDParam}(nworkers())

		for j=1:nworkers()
		pFor[j] = 	getIneqconstrainedOEDParam(OEDparam[j],L,Ce,ce[:,ids.==j],Ci,ci[:,ids.==j],
		ftolIPQP=1e-3,mutolIPQP=1e-3,outIPQP=0)
		end
	else
		error("Invalid IPmode = $IPmode")
	end
	display(pFor[1])

	np = length(m0)

	if nworkers()>1
		pMis  = Array{RemoteChannel}(nworkers())
		for j=1:nworkers()
			pMis[j] = initRemoteChannel(getMisfitParam,workers()[j], pFor[j],ones(size(OEDparam[j].Ftrue)),Xtrue[:,ids.==j],SSDFun)
		end
	else
		pMis  = getMisfitParam(pFor[1],ones(size(Xtrue)),Xtrue,SSDFun);
	end

	# make sure the full size problem is in there
	Minv        = getRegularMesh([0  1 0 1],[12 12]);
	pInv        = getInverseParam(Minv,identityMod,smallnessReg,0.,zeros(np),0*ones(np),180*ones(np));
	pInv.pcgMaxIter = pcgMaxIter
	pInv.maxIter    = maxIter
	pInv.maxStep    = maxStep
	pInv.HesPrec.applyPrec = (a,b,x,c) -> x
	display(pInv)

	println("\n\n=== optimize projection angles (started=$(now())) ===")
	mIter = zeros(length(m0),maxIter+1)
	Diter = zeros(size(Xtrue,1),size(Xtrue,2),maxIter+1)
	if nworkers()==1
		storeInterm = (mc,Dc,iter,pInv,pMis) -> (mIter[:,iter+1] = mc; Diter[:,:,iter+1] = Dc);
	else
		storeInterm  = (mc,Dc,iter,pInv,pMis) -> (mIter[:,iter+1] = mc; for k=1:nworkers(); Diter[:,ids.==k,iter+1] = fetch(Dc[k]); end);
	end
	tic()
	mOpt,Xo,flag,His = projGN(copy(m0),pInv,pMis,dumpResults=storeInterm,solveGN=projGNexplicit);

	Xopt = 0*Xtrue
	if nworkers()>1
		Xopt = Xo
	else
		for k=1:nworkers()
			Xopt[:,ids.==k] = fetch(Xo[k])
		end
	end

	time = toc()
	println("=== finished in $time secs ===")
	if !isempty(exName)
		fileName = @sprintf "%s-%s.mat" exName IPmode
		logName = @sprintf "%s-%s.log" exName IPmode
		mv("jInv.out",logName,remove_destination=true)
		println("=== save results in $fileName===")
		His.Dc = [];
		results = Dict("mOpt"=>mOpt, "Xopt"=>Xopt,"flag"=>flag,"His"=>His,"time"=>time,"mIter"=>mIter,"Diter"=>Diter)
		matwrite(fileName,results)
	end
	return mOpt
end
