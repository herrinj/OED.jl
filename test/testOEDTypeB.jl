using jInv.Mesh
using jInv.Utils
using Base.Test
using OED

r = 4;
l = 3;
Mimg  = getRegularMesh([-1 1 -1 1.],[4,4])
Mgrid = getRegularMesh([-1.2 2.4 -1.6 1.4],[8,8])

Ftrue = randn(Mimg.nc,10)
OEDparam = getOEDTypeBParam(r,l,Ftrue,Mimg,Mgrid,noiseLevel=0.)
display(OEDparam)

T = getTransmissionMatrix(OEDparam);

(Nt,nt,nrhst) = getProblemSize(OEDparam)
@test getNumberOfDesignParameters(OEDparam)==l

println("--- test getM and derivatives ---")
v   = randn(Mimg.nc)
p0  = rand(getNumberOfDesignParameters(OEDparam))
pass, = checkDerivative(p-> getM(p,OEDparam)*v, p-> getdpMx(p,v,OEDparam), p0,out=true)
@test pass


println("--- test getM^T and derivatives ---")
v   = randn(Nt)
p0  = rand(getNumberOfDesignParameters(OEDparam))
pass, = checkDerivative(p-> getM(p,OEDparam)'*v, p->getdpMTx(p,v,OEDparam), p0,out=false)
@test pass

println("--- test getD and derivatives ---")
pass, = checkDerivative(p->getD(p,OEDparam)[:,1], p->getdpD(p,1,OEDparam),p0,out=false)
@test pass




