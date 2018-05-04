using MAT
using jInv.Mesh
using jInv.LinearSolvers
using jInv.ForwardShare
using jInv.InverseSolve
using jInv.Utils
using OED

model = 1;
IPmodes = [0;1;2;3];
lambda = .1;

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
ids = mod.(0:size(Xtrue,2)-1,nworkers())+1
OEDparam = Array{OEDTypeBParam}(nworkers())
for k=1:nworkers()
	OEDparam[k] = getOEDTypeBParam(r,l,Xtrue[:,ids.==k],Mimg,Mgrid)
	OEDparam[k].noiseLevel = 1e-5
end

ipDict  = Dict(0=>:unconstr, 1=>:eqconstr, 2=>:nonnegconstr, 3=>:boxconstr)

for k=1:length(IPmodes)
	IPmode  = ipDict[IPmodes[k]]
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

	if nworkers()>1
		pMis  = Array{RemoteChannel}(nworkers())
		for j=1:nworkers()
			pMis[j] = initRemoteChannel(getMisfitParam,workers()[j], pFor[j],ones(size(OEDparam[j].Ftrue)),Xtrue[:,ids.==j],SSDFun)
		end
	else
		pMis  = getMisfitParam(pFor[1],ones(size(Xtrue)),Xtrue,SSDFun);
	end
	w1,w2 = getCellCenteredAxes(wGrid)
	F     = zeros(length(w1),length(w2))
	for k1=1:length(w1)
	    println("k1=$k1")
	    for k2=1:k1
	        Dc,F[k1,k2],dF,d2F,pMis,tMis = computeMisfit([w1[k1],w2[k2]],pMis);
	    end
	end
	outfile = @sprintf "%s-%s.mat" resfile IPmode
	matwrite(outfile,Dict("F"=>F,"domain"=>wGrid.domain,"n"=>wGrid.n))
end
