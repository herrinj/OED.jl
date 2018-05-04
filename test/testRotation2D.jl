using OED
using Base.Test
using jInv.Utils
using jInv.Mesh


domain = [0 10 0 2.]; n = [8 9]
Mr     = getRegularMesh(domain,n)
xc     = getCellCenteredGrid(Mr)
w = 22/pi;
c = vec((domain[2:2:end]+domain[1:2:end])./2)

y = rotation2D(w,xc,c=c)

function testFun(w,v=[])
	y,dy = rotation2D(w,xc,c=c,doDerivative=true)
	if isempty(v)
		return y
	else
		return y,dy*v
	end
end
chkDer, = checkDerivative(testFun,0.1)
@test chkDer
