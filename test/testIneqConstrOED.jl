using jInv.Utils
using Base.Test
using jInv.ForwardShare
using jInv.Mesh
using OED

n = 100
nrhs = 2
# A = abs(sprandn(N,n,.2))
A = sdiag(1./linspace(1,n,n).^4)
A = A[1:2:end,:]
N = size(A,1)
x = rand(n,nrhs); x[4:7,:] = 1;
D = A*x  + 1e-1* randn(N,nrhs)/1e3
OEDparamA = OEDTypeAParam(A,D,zeros(0,0))

# setup OED type B param
r = 20;
l = 8;
Mimg  = getRegularMesh([-1 1 -1 1.],[10,10])
Mgrid = getRegularMesh([-1.2 2.4 -1.6 1.4],[20,20])
OEDparamB = getOEDTypeBParam(r,l,x,Mimg,Mgrid,noiseLevel=0.)

OEDs = (OEDparamA, OEDparamB)



println("--- solve inequality constrained problem ---")
bounds = [.3 .9]
L = .1*speye(n)
Ce = ones(1,n)
ce = Ce*x
Ci = [speye(n);-speye(n)]
ci = [bounds[1]*ones(n,nrhs); -bounds[2]*ones(n,nrhs)]
pForBox    = getIneqconstrainedOEDParam(OEDparamA, L,Ce,ce,Ci,ci,maxIterIPQP=40,outIPQP=-1)  
pForNonNeg = getIneqconstrainedOEDParam(OEDparamA, L,[],[],speye(n),0*ones(n,nrhs),maxIterIPQP=40,outIPQP=-1)  


pFors = (pForBox, pForNonNeg)
for j=1:length(OEDs)
    for k=1:length(pFors)
		pFors[k].Design = OEDs[j]
		display(pFors[k])
		
		np  =  getNumberOfDesignParameters(OEDs[j])
		p0      = rand(np)
		xr, = getData(p0,pFors[k])
		@test norm(pFors[k].Ce*xr-pFors[k].ce) < 1e-12
		@test all(pFors[k].Ci*xr.>= pFors[k].ci)
		println("rel.error: $(norm(xr-x)/norm(x))")
		
		println("--- test derivative of LS solution w.r.t. design ---")
		checkDerivative(rand(np),pFors[k])	
		
		println("--- adjoint test ---")
		(nw,nv) = getSensMatSize(pFors[k])
		v = randn(nv)
		w = randn(nw)
		t1 = dot(w,getSensMatVec(v,p0,pFors[k]))
		t2 = dot(v,getSensTMatVec(w,p0,pFors[k]))
		@test abs(t1-t2)/abs(t1) < 1e-10
	end
end