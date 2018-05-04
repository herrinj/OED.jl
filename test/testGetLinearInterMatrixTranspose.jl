using OED
using Base.Test
using jInv.Mesh
using jInv.Utils


domain = [-1 1 -1 1.]; n = [8 9]
Mr     = getRegularMesh(domain,n)
xc     = .8*getCellCenteredGrid(Mr)

dataT     = exp.(-xc[:,1].^2- xc[:,2].^2)
dataT     = reshape(dataT,tuple(n...))

Img = randn(prod(n))
At,dAt = getLinearInterMatrixTranspose(domain,n,xc,doDerivative=true,Img=Img)
@test norm(sum(At,1)-1,Inf)<1e-12

function linearInterMatTranspTestFun(x,v=[])
	At,dAT = getLinearInterMatrixTranspose(domain,n,x,doDerivative=true,Img=Img)

	if isempty(v)
		return At*Img
	else
		return At*Img,dAT*v
	end
end
x0 = vec(xc+1e-2*randn(size(xc)))
chkDer, = checkDerivative(linearInterMatTranspTestFun,x0)
@test chkDer

x,y = getCellCenteredAxes(Mr)
A0 = getLinearInterMatrix(domain,n,copy(x0))
A0t= getLinearInterMatrixTranspose(domain,n,x0)
@test norm(A0-A0t',Inf) <= 1e-10
