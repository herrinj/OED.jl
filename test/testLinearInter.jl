using OED
using Base.Test
using jInv.Mesh
using jInv.Utils

domain = [-1 1 -1 1.]; n = 4*[8 8]
Mr     = getRegularMesh(domain,n)
x,y    = getCellCenteredAxes(Mr)
xc     = getCellCenteredGrid(Mr)

dataT     = exp.(-xc[:,1].^2- xc[:,2].^2)
dataT     = reshape(xc*rand(2),tuple(n...))

Tc,dT = linearInter(dataT,domain,xc,doDerivative=true)
A0 = getLinearInterMatrix(domain,n,xc)
Tt = A0*vec(dataT)
@test norm(vec(Tc)-vec(dataT))/norm(vec(dataT)) < 1e-8
@test norm(vec(Tt)-vec(dataT))/norm(vec(dataT)) < 1e-8

function linearInterTestFun(x,v=[])
	Tc,dT = linearInter(dataT,domain,x,doDerivative=true)

	if isempty(v)
		return Tc
	else
		return Tc,dT*v
	end
end
x0 = vec(xc+1e-3*randn(size(xc)))
chkDer, = checkDerivative(linearInterTestFun,x0)

@test chkDer
println("check getLinearInterMatrix with points in domain")
xt = reshape(copy(x0),Mr.nc,2);
 xt[:,1]=xc[:,1];  xt[:,2]*=2
A = getLinearInterMatrix(copy(domain),copy(n),copy(xt))
T = linearInter(copy(dataT),copy(domain),copy(xt))
Tt = A*vec(dataT)
@test norm(Tt-T)/norm(T) < 1e-10
println("check getLinearInterMatrix with points out of domain")
A = getLinearInterMatrix(.5*domain,n,x0)
T = linearInter(dataT*0+1,.5*domain,x0)
Tt = A*vec(0*dataT+1)
@test norm(Tt-T)/norm(T) < 1e-10
