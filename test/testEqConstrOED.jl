using jInv.Utils
using Base.Test
using OED
using jInv.ForwardShare
using jInv.Mesh

n = 100
nrhs = 2
# A = abs(sprandn(N,n,.2))
A = sdiag(1./linspace(1,n,n).^4)
A = A[1:end-2,:]
N = size(A,1)
x = rand(n,nrhs); x[4:7,:] = 1;
D = A*x  + 0* randn(N,nrhs)/1e3
OEDparamA = OEDTypeAParam(A,D,zeros(0,0))

# setup OED type B param
r = 20;
l = 8;
Mimg  = getRegularMesh([-1 1 -1 1.],[10,10])
Mgrid = getRegularMesh([-1.2 2.4 -1.6 1.4],[20,20])
OEDparamB = getOEDTypeBParam(r,l,x,Mimg,Mgrid,noiseLevel=0.)

OEDs = (OEDparamA, OEDparamB)

for k=1:2
	display(OEDs[k])
	
	L = 1*speye(n)
	Ce = ones(1,n)
	ce = Ce*x
	pFor = EqconstrainedOEDParam(OEDs[k], L,Ce,ce, zeros(0,0), [])  
	
	println("--- solve LS problem ---")
	display(pFor)
	np  =  getNumberOfDesignParameters(OEDs[k])
	p0      = rand(np)
	xr, = getData(p0,pFor)
	@test norm(Ce*xr-ce) < 1e-12
	println("rel.error: $(norm(xr-x)/norm(x))")

	println("--- test derivative of LS solution w.r.t. design ---")
	pass, = checkDerivative(rand(np),pFor)
	@test pass

	println("--- adjoint test ---")
	nw,nv = getSensMatSize(pFor)
	v = randn(nv)
	w = randn(nw)
	t1 = dot(w,getSensMatVec(v,p0,pFor))
	t2 = dot(v,getSensTMatVec(w,p0,pFor))
	@test abs(t1-t2)/abs(t1) < 1e-10
end


