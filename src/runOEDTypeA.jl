export runOEDTypeA

"""
function runOEDTypeA

Input:

	param
	Xtrue
	L
	n
	theta
	Alphas
	model
	IPmode
"""
function runOEDTypeA(OEDparam::Array{OEDTypeAParam},Xtrue::Array{Float64},L::AbstractArray,Alphas,modFun::Function,IPmode::Symbol,m0::Vector{Float64};Ce=[],ce=[],Ci=[],ci=[],pcgMaxIter=10,maxIter=20,maxStep=1.0,exName="")

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
				pFor[j] = getIneqconstrainedOEDParam(OEDparam[j],L,Ce,ce[:,ids.==j],I,0*Xtrue[:,ids.==j],
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
	A  = copy(pFor[1].Design.A)
	D  = zeros(size(A,1),size(Xtrue,2))
	for k=1:nworkers()
		D[:,ids.==k] = pFor[k].Design.D
	end

	for k=1:length(Alphas)
		alpha = Alphas[k]
		println("\n\n==== processing alpha=($alpha) which is $(k) out of $(length(Alphas)) ====\n\n")
		println("\n\n=== setup first optimization to identify non-zeros ===")
		if nworkers()>1
			pMis  = Array{RemoteChannel}(nworkers())
			for j=1:nworkers()
				pFor[j].Design.A = A
				pFor[j].Design.D = D[:,ids.==j]
				pMis[j] = initRemoteChannel(getMisfitParam,workers()[j], pFor[j],ones(size(Xtrue[:,ids.==j])),Xtrue[:,ids.==j],SSDFun)
			end
		else
			pFor[1].Design.A = A
			pFor[1].Design.D = D
			pMis        = getMisfitParam(pFor[1],ones(size(Xtrue)),Xtrue,SSDFun);
		end

		# make sure the full size problem is in there
		Minv        = getRegularMesh([0  1 0 1],[12 12]);
		pInv        = getInverseParam(Minv,modFun,l1Reg,alpha,zeros(np),zeros(np),Inf*ones(np));
		pInv.pcgMaxIter = pcgMaxIter
		pInv.maxIter    = maxIter
		pInv.maxStep    = maxStep
		display(pInv)

		println("\n\n=== first inversion started $(now()) ===")
		tic()
		mOpt1,Xo,flag1,His1 = projGN(copy(m0),pInv,pMis);
		time1 = toc()
		wOpt1, = modFun(mOpt1)

		Xopt1 = 0*Xtrue
		if nworkers() == 1
			Xopt1 = Xo
		else
			for jj=1:nworkers()
				Xopt1[:,ids.==jj] = fetch(Xo[jj])
			end
		end

		println("model:   nnz=$(countnz(mOpt1)), sparsity=$(countnz(mOpt1)/length(mOpt1))")
		println("weights: nnz=$(countnz(wOpt1)), sparsity=$(countnz(wOpt1)/length(wOpt1))")

		println("\n\n=== setup second optimization to tune nonzero weights ===")
		if nworkers()>1
			pMis  = Array{RemoteChannel}(nworkers())
			for j=1:nworkers()
				pFor[j].Design.A = A[wOpt1.>0,:]
				pFor[j].Design.D = D[wOpt1.>0,ids.==j]
				pMis[j] = initRemoteChannel(getMisfitParam,workers()[j], pFor[j],ones(size(Xtrue[:,ids.==j])),Xtrue[:,ids.==j],SSDFun)
			end
		else
			pFor[1].Design.A = A[wOpt1.>0,:]
			pFor[1].Design.D = D[wOpt1.>0,:]
			pMis        = getMisfitParam(pFor[1],ones(size(Xtrue)),Xtrue,SSDFun);
		end

		nnzM = countnz(mOpt1)
		pInv.boundsLow  = zeros(nnzM)
		pInv.boundsHigh = Inf*ones(nnzM)
		pInv.mref       = zeros(nnzM)
		pInv.alpha      = 0.0
		pInv.HesPrec.applyPrec = (a,b,x,c) -> x
		display(pInv)

		println("\n\n=== second inversion started $(now()) ===")
		tic()
		mOpt2,Xo,flag2,His2 = projGN(mOpt1[mOpt1.>0],pInv,pMis);
		time2 = toc()
		wOpt2, = modFun(mOpt2)
		Xopt2 = 0*Xtrue
		if nworkers() == 1
			Xopt2 = Xo
		else
			for jj=1:nworkers()
				Xopt2[:,ids.==jj] = fetch(Xo[jj])
			end
		end

		if !isempty(exName)
			fileName = @sprintf "%s-%s-alpha-%2.2e-nnz-%d.mat" exName IPmode alpha countnz(mOpt2)
			logName = @sprintf "%s-%s-alpha-%2.2e-nnz-%d.log" exName IPmode alpha countnz(mOpt2)
			mv("jInv.out",logName,remove_destination=true)
			println("=== save results in $fileName===")
			His1.Dc = []; His2.Dc = []
			results = Dict("mOpt1"=>mOpt1, "Xopt1"=>Xopt1,"flag1"=>flag1,"His1"=>His1,"time1"=>time1,"wOpt1"=>wOpt1,
			               "mOpt2"=>mOpt2, "Xopt2"=>Xopt2,"flag2"=>flag2,"His2"=>His2,"time2"=>time2,"wOpt2"=>wOpt2)
			matwrite(fileName,results)
		end
		m0 = mOpt1;
	end
	return m0
end
