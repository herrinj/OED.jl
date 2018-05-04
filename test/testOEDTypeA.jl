using jInv.Utils
using Base.Test
using OED

n = 10
N = 10
nrhs = 1

A = sdiag(1./linspace(1,N,N).^4)
x = rand(n,nrhs); x[4:7,:] = 1;
D = A*x  + 1e-8* randn(N,nrhs)/1e3

OEDparam = OEDTypeAParam(A,D,zeros(0,0))

display(OEDparam)

(Nt,nt,nrhst) = getProblemSize(OEDparam)
@test Nt==N
@test nt==n
@test nrhst ==nrhs
@test getNumberOfDesignParameters(OEDparam)==N

println("--- test getM and derivatives ---")
v = randn(n)
p0  = rand(N)
pass, = checkDerivative(p-> getM(p,OEDparam)*v, p-> getdpMx(p,v,OEDparam), p0)
@test pass

println("--- test getM^T and derivatives ---")
v = randn(N)
p0  = rand(N)
pass, = checkDerivative(p-> getM(p,OEDparam)'*v, p->getdpMTx(p,v,OEDparam), p0)
@test pass

println("--- test getD and derivatives ---")
pass, = checkDerivative(p->getD(p,OEDparam)[:,1], p->getdpD(p,1,OEDparam),p0)
@test pass


