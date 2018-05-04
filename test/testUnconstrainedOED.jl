using OED
using jInv.Utils
using jInv.ForwardShare
using Base.Test
using jInv.Mesh

n = 100
nrhs = 5

A = sdiag(1./linspace(1,n,n).^4)
A = A[1:2:end,:]
N = size(A,1)
x = rand(n,nrhs); x[4:7,:] = 1;
D = A*x  + 1e-8* randn(N,nrhs)/1e3
OEDparamA = OEDTypeAParam(A,D,zeros(0,0))

# setup OED type B param
r = 20;
l = 2;
Mimg  = getRegularMesh([-1 1 -1 1.],[10,10])
Mgrid = getRegularMesh([-1.2 2.4 -1.6 1.4],[20,20])
OEDparamB = getOEDTypeBParam(r,l,x,Mimg,Mgrid,noiseLevel=0.)
OEDparamC = getOEDTypeCParam(r,l,x,Mimg,Mgrid,noiseLevel=0.)


L = 1e-1*speye(n)
pForA = UnconstrainedOEDParam(OEDparamA,L,[],[])
pForB = UnconstrainedOEDParam(OEDparamB,L,[],[])
pForC = UnconstrainedOEDParam(OEDparamC,0*L,[],[])

pFors = (pForA,pForB,pForC)

println("\n\t---test OED with unconstrained lower-level problem---")
for k=1:3
	display(pFors[k])
	np  =  getNumberOfDesignParameters(pFors[k].Design)
	p0      = rand(np)

	println("\t\tcheck reconstruction")
	xr, = getData(p0,pFors[k])
	println("\t\t\trel.error: $(norm(xr-x)/norm(x))")
	
	print("\t\tcheckDerivative...")
	
	pass, = checkDerivative(p0,pFors[k])
	@test pass
	print("passed!\n")
	
	print("\t\tadjoint test...")
	(nw,nv) = getSensMatSize(pFors[k])
	w = rand(nw)
	v = rand(nv)
	t1 = dot(w,getSensMatVec(v,p0,pFors[k]))
	t2 = dot(v,getSensTMatVec(w,p0,pFors[k]))
	@test norm(t1-t2)/norm(t1)<1e-10
	print("passed!\n")
end
println("\t---least squares IP: sensitivities OK ---")
