# using OED
using jInv.Utils
using Base.Test
using OED

println("=== solve 2D QP ===")
n = 2

G = eye(n)
c = 0*[1;1.]

Ae = [1 1.]
be = [1.2]
me = size(Ae,1)

Ai = eye(n)
bi = 0.0*ones(n)
mi = size(Ai,1)

x0  = rand(n)+3.4
y0  = rand(mi)+0*1.2
le0 = rand(me)
li0 = rand(mi)

println("\n\t---test Interior Point Method for QP ---")

print("\t\tsolve 2D QP with inequality and equality constraints...\n\t\t")
xas,yas,leas,lias,hisas = ipqp(G,c,Ae,be,Ai,bi,x0,y0,le0,li0,
			maxIter=20,ftol=1e-30,solveKKT=OED.solveFullKKT)
@test all(Ai*xas-bi .>0)
# @test norm(xas-0) < 1e-3
@test all(yas .> 0)
@test all(lias .> 0)
@test norm(Ai*xas-yas-bi) < 1e-14
@test maximum(yas.*lias) < 1e-2
@test maximum(G*xas - Ae'*leas - Ai'*lias + c)< 1e-15

xas,yas,leas,lias,hisas = ipqp(G,c,Ae,be,Ai,bi,x0,y0,le0,li0,
			maxIter=20,ftol=1e-30)
@test all(Ai*xas-bi .>0)
# @test norm(xas-0) < 1e-3
@test all(yas .> 0)
@test all(lias .> 0)
@test norm(Ai*xas-yas-bi) < 1e-14
@test maximum(yas.*lias) < 1e-2
@test maximum(G*xas - Ae'*leas - Ai'*lias + c)< 1e-15

xas,yas,leas,lias,hisas = ipqp(G,c,Ae,be,Ai,bi,x0,y0,le0,li0,
			maxIter=20,ftol=1e-30,solveKKT=OED.solveAugKKT)
@test all(Ai*xas-bi .>0)
# @test norm(xas-0) < 1e-3
@test all(yas .> 0)
@test all(lias .> 0)
@test norm(Ai*xas-yas-bi) < 1e-14
@test maximum(yas.*lias) < 1e-2
@test maximum(G*xas - Ae'*leas - Ai'*lias + c)< 1e-15

print("test passed!\n")

print("\t\tsolve 2D QP without equality constraints...\n\t\t")
n = 2

G = eye(n)
c = 0*[1;1.]

Ae = zeros(0,n)
be = zeros(0)
me = size(Ae,1)

Ai = eye(n)
bi = 0.0*ones(n)
mi = size(Ai,1)

x0  = rand(n)+3.4
y0  = rand(mi)+0*1.2
le0 = rand(me)
li0 = rand(mi)

xas,yas,leas,lias,hisas = ipqp(G,c,Ae,be,Ai,bi,x0,y0,le0,li0,out=1,
			maxIter=20,ftol=1e-30,storeInterm=true,solveKKT=OED.solveFullKKT)
@test all(Ai*xas[:,end]-bi .>0)
# @test norm(xas-0) < 1e-3
@test norm(Ai*xas[:,end]-yas[:,end]-bi) < 1e-14
@test maximum(yas[:,end].*lias[:,end]) < 1e-2
@test all(yas .> 0)
@test all(lias .> 0)
@test maximum(G*xas[:,end] -  Ai'*lias[:,end] + c)< 1e-15

xas,yas,leas,lias,hisas = ipqp(G,c,Ae,be,Ai,bi,x0,y0,le0,li0,out=1,
			maxIter=20,ftol=1e-30,storeInterm=true)
@test all(Ai*xas[:,end]-bi .>0)
# @test norm(xas-0) < 1e-3
@test norm(Ai*xas[:,end]-yas[:,end]-bi) < 1e-14
@test maximum(yas[:,end].*lias[:,end]) < 1e-2
@test all(yas .> 0)
@test all(lias .> 0)
@test maximum(G*xas[:,end] -  Ai'*lias[:,end] + c)< 1e-15

xas,yas,leas,lias,hisas = ipqp(G,c,Ae,be,Ai,bi,x0,y0,le0,li0,out=1,
			maxIter=20,ftol=1e-30,storeInterm=true,solveKKT=OED.solveAugKKT)
@test all(Ai*xas[:,end]-bi .>0)
# @test norm(xas-0) < 1e-3
@test norm(Ai*xas[:,end]-yas[:,end]-bi) < 1e-14
@test maximum(yas[:,end].*lias[:,end]) < 1e-2
@test all(yas .> 0)
@test all(lias .> 0)
@test maximum(G*xas[:,end] -  Ai'*lias[:,end] + c)< 1e-15
print("test passed!\n")
println("\n\t---[IPQP] all tests passed---")


