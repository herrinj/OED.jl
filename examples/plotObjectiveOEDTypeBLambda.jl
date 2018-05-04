using MAT
using jInv.Mesh
using jInv.LinearSolvers
using jInv.ForwardShare
using jInv.InverseSolve
using jInv.Utils
using OED

model   = 1;
IPmodes = [0;1;2;3];
lambdas = logspace(-4.,2,20);
nrep    = 4;

r       = 48
l       = 2
pad     = 4
nex     = 20
matfile = "data/OED_shepp.mat"
resfile = "results/plotObjectiveB/OED_shepp"
matfile = matread(matfile)
n       = Int(matfile["n"])
L       = .1*speye(n*n)
Xtrue   = matfile["x_true"]

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

Mimg   = getRegularMesh([-1 1 -1 1.],[n,n])
domain = copy(Mimg.domain);
domain[1:2:end] -= pad*Mimg.h
domain[2:2:end] += pad*Mimg.h
Mgrid = getRegularMesh(domain,[n+2*pad,n+2*pad])
wGrid = getRegularMesh([0 180 0 180],[32 32]);


# divide problems among available workers
ids = mod(0:size(Xtrue,2)-1,nworkers())+1
OEDparam = Array{OEDTypeBParam}(nworkers())
for k=1:nworkers()
	OEDparam[k] = getOEDTypeBParam(r,l,Xtrue[:,ids.==k],Mimg,Mgrid)
	OEDparam[k].noiseLevel = 1e0
end
ipDict  = Dict(0=>:unconstr, 1=>:eqconstr, 2=>:nonnegconstr, 3=>:boxconstr)

pFor = 1;
for k=3:3
	IPmode  = ipDict[IPmodes[k]]
	outfile = @sprintf "%s-%s.mat" resfile IPmode
	outf    = matread(outfile)
    F       = outf["F"]
	F[F.==0] = NaN
	w1,w2 = getCellCenteredAxes(wGrid)
	i,j     = ind2sub((length(w1),length(w2)), indmin(F))
	wOpt = [w1[i];w2[j]]

	println("wOpt = $wOpt")
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
				          ftolIPQP=1e-3,mutolIPQP=1e-3,outIPQP=1)
		end

	elseif IPmode==:boxconstr
		if isempty(Ci) || isempty(ci)
			error("inequality constraints must be provided. To run equality constrained OED, use IPmode=:eqconstr.")
		end
		pFor = Array{IneqconstrainedOEDParam}(nworkers())

		for j=1:nworkers()
		pFor[j] = 	getIneqconstrainedOEDParam(OEDparam,L,Ce,ce[:,ids.==j],Ci,ci[:,ids.==j],
		ftolIPQP=1e-3,mutolIPQP=1e-3,outIPQP=0)
		end
	else
		error("Invalid IPmode = $IPmode")
	end
	display(pFor[1])


	F    = zeros(length(lambdas),nrep)
	dF   = zeros(length(lambdas),nrep)
	Dc    = zeros(size(Xtrue,1),size(Xtrue,2),length(lambdas))

	for j=1:nrep
	for k=1:length(lambdas)
	    println("lambda=$(lambdas[k])")
			if nworkers()>1
				pMis  = Array{RemoteChannel}(nworkers())
				for j2=1:nworkers()
					pFor[j2].L = sqrt(lambdas[k])*speye(n*n)
					pMis[j2] = initRemoteChannel(getMisfitParam,workers()[j2], pFor[j2],ones(size(OEDparam[j2].Ftrue)),Xtrue[:,ids.==j2],SSDFun)
				end
		        Dcr,F[k,j],dFk,d2F,pMis,tMis = computeMisfit(wOpt,pMis);
				dF[k,j] = norm(dFk)
				for j2=1:nworkers()
					Dc[:,ids.==j2,k] = fetch(Dcr[j2])
				end
			else
				pFor[1].L = sqrt(lambdas[k])*speye(n*n)
				pMis  = getMisfitParam(pFor[1],ones(size(Xtrue)),Xtrue,SSDFun);
		        Dc[:,:,k],F[k,j],dFk,d2F,pMis,tMis = computeMisfit(wOpt,pMis);
				dF[k,j] = norm(dFk)
			end

	end
end
	outfile2 = @sprintf "%s-lambdas-%s.mat" resfile IPmode
	matwrite(outfile2,Dict("F"=>F,"lambdas"=>lambdas,"Dc"=>Dc,"dF"=>dF,"wOpt"=>wOpt))
end
