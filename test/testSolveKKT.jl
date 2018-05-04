using Base.Test
using OED

println("\t\t--- test different KKT solvers ---")
n = 2

G = eye(n)
c =randn(n)

Ae = [1 1.]
be = 1.2
me = size(Ae,1)

Ai = eye(n)
bi = rand(n)
mi = size(Ai,1)

x0  = rand(n)+3.4
y0  = rand(mi)+1.2
le0 = rand(me)
li0 = rand(mi)

# solveFullKKT(G,Ae,Ai,y,li,rd,rpe,rpi,rli;KKTinv=[],doClear=true)

tau0 = randn()
rd,rpe,rpi,rli, = OED.getF(G,c,Ae,[be],Ai,bi,x0,y0,le0,li0,tau=tau0)
d1, = OED.solveFullKKT(G,Ae,Ai,y0,li0,rd,rpe,rpi,rli)
(dx,dy,de,di)   = OED.splitVecIPQP(d1,n,1,n)
@test norm(Ae*(x0+dx) - be) < 1e-14
d2, = OED.solveAugKKT(G,Ae,Ai,y0,li0,rd,rpe,rpi,rli)
(dx,dy,de,di)   = OED.splitVecIPQP(d2,n,1,n)
@test norm(Ae*(x0+dx) - be) < 1e-14
d3, = OED.solveNormalEqKKT(G,Ae,Ai,y0,li0,rd,rpe,rpi,rli)
(dx,dy,de,di)   = OED.splitVecIPQP(d3,n,1,n)
@test norm(Ae*(x0+dx) - be) < 1e-14

@test norm(d1-d2)/norm(d1)<1e-14
@test norm(d1-d3)/norm(d1)<1e-14

println("\t--- KKT solvers OK ---")
