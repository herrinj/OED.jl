using OED
using Base.Test
using jInv.InverseSolve


n = 10
nrhs = 1
theta = 1:10
Alphas = logspace(-10,-9,3);

A = sdiag(1./linspace(1,n,n).^4)
N = 5;
x = rand(n,nrhs); x[4:7,:] = 1;
D = A*x  + 1e-8* randn(size(A,1),nrhs)/1e3
OEDparam = OEDTypeAParam(A,D,zeros(0,0))

L = speye(n)

modSettings = [(v->parallelTomoMod(v,N), rand(2));
			(identityMod,rand(size(A,1)))]

IPmodes = (:unconstr,:nonnegconstr,:boxconstr)

Ce = ones(1,n); ce = Ce*x;
Ci = speye(n);  ci = zeros(n,nrhs)
for ip=IPmodes
	for mod=modSettings
runOEDTypeA([OEDparam],x,L,Alphas,mod[1],ip,mod[2];Ce=Ce,ce=ce,Ci=Ci,ci=ci,pcgMaxIter=10,maxIter=2,maxStep=1.0,exName="")
	end
end