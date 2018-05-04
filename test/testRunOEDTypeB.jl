using OED
using Base.Test
using MAT
using jInv.Mesh

r = 3;
l = 3;
nrhs = 10
Mimg  = getRegularMesh([-1 1 -1 1.],[4,4])
Mgrid = getRegularMesh([-1.2 2.4 -1.6 1.4],[8,8])

Ftrue = randn(Mimg.nc,nrhs).^2
for k=1:nrhs
	Ftrue[k,k] = 10;
end
OEDparam = getOEDTypeBParam(r,l,Ftrue,Mimg,Mgrid)
display(OEDparam)
L = speye(Mimg.nc)

IPmodes = (:unconstr,:nonnegconstr,:nonnegconstr,:boxconstr)

Ce = ones(1,Mimg.nc); ce = Ce*Ftrue;
Ci = [speye(Mimg.nc);-speye(Mimg.nc)];
ci = [minimum(Ftrue)*ones(Mimg.nc,nrhs); -maximum(Ftrue)*ones(Mimg.nc,nrhs)]
m0 = rand(l)
mOpt = []
for ip=IPmodes
mO = runOEDTypeB([OEDparam],Ftrue,L,ip,copy(m0);Ce=Ce,ce=ce,Ci=Ci,ci=ci,pcgMaxIter=10,maxIter=3,maxStep=10.0,exName="/tmp/OEDb")
push!(mOpt,mO)
end
