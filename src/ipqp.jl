export ipqp

"""
	function stepLength(x,y,dx,dy;eta)

	compute step length such that x+alpha*dx and y+alpha*dy are feasible

	then (optional) backtrack a little by multiplying with eta
"""
function stepLength{T<:AbstractFloat}(x::Vector{T},y::Vector{T},dx::Vector{T},dy::Vector{T};eta::T=0.999)
    alphax = -1./min(minimum(dx./x),-1)
    alphax = min(1,alphax)

    alphay = -1./min(minimum(dy./y),-1)
    alphay = min(1,alphay)

    return eta*min(alphax,alphay)
end

"""
	function splitVecIPQP(v,n,me,mi)

	get parts of v=[x;y;le;li]. used in ipqp.
"""
function splitVecIPQP(v::Vector{Float64},n::Int64,me::Int64,mi::Int64)
    x   = v[1:n]
    y   = v[n+(1:mi)]
    le  = v[n+mi+(1:me)]
	li  = v[n+mi+me+1:end]
    return x,y,le,li
end

"""
	function getF(G,c,Ae,be,Ai,bi,x,y,le,li;tau::Real=0.0)

	gets right hand side of KKT system in interior point method

"""
function getF{T}(G::AbstractArray{T},c::Vector{T},Ae::AbstractArray{T},be::Vector{T},Ai::AbstractArray{T},bi::Vector{T},
	             x::Vector{T},y::Vector{T},le::Vector{T},li::Vector{T};
	             dy::Vector{T}=zeros(T,length(y)),dli::Vector{T}=zeros(T,length(li)),tau::T=0.0)
    rd    = G*x - Ae'*le - Ai'*li + c
    rpe   = Ae*x - be
    rpi   = Ai*x - y - bi
    rli   = y.*li + dy.*dli - tau
    return rd,rpe,rpi,rli,y.*li
end

"""
	function getJ(G,Ae,Ai,y,li)

	gets KKT matrix for interior point method for solving

	      | G   0   -Ae'  -Ai'  |
	      | Ae  0    0     0    |
	KKT = | Ai -I    0     0    |
	      | 0   Li   0     Y    |

	where Li=sdiag(li) and Y=sdiagg(y) are diagonal matrices

		min_x .5*x'*G*c + c'*x  subject to Ae*x = be and Ai*x >= bi

"""
function getJ{T}(G::AbstractArray{T},Ae::AbstractArray{T},Ai::AbstractArray{T},y::Vector{T},li::Vector{T})
    n  = size(G,1)
    me = size(Ae,1)
    mi = size(Ai,1)

    In    = speye(n)
    Imi   = speye(mi)
    Znmi  = spzeros(n,mi)
	Znme  = spzeros(n,me)
	Znn   = spzeros(n,n)
    Zmimi = spzeros(mi,mi)
	Zmeme = spzeros(me,me)
	Zmime = spzeros(mi,me)

    Y   = sdiag(y)
    Li  = sdiag(li)

    return [G Znmi -Ae' -Ai'; Ae Zmime' Zmeme Zmime'; Ai -Imi Zmime Zmimi; Znmi' Li Zmime Y]
end

"""
	function solveNormalEqKKT(G,c,Ae,be,Ai,bi,x,y,le,li,rd,rpe,rpi,rli;KKTinv=[],doClear=true)

	brings KKT system to \"normal equation\" form and solves it.

"""
function solveNormalEqKKT{T}(G::AbstractArray{T},Ae::AbstractArray{T},Ai::AbstractArray{T},                  y::Vector{T},li::Vector{T},rd::Vector{T},rpe::Vector{T},rpi::Vector{T},rli::Vector{T};KKTinv=cholfact(speye(1)),doClear=true)

	# build normal eq matrix
    n  = size(G,1)
    me = size(Ae,1)
    mi = size(Ai,1)


 	liy = sdiag(li./y)
	yLi = sdiag(y./li)

    rh1  = -(rd+Ai'*liy*(rpi+(rli./li)));
	rh2  = -rpe

	if size(KKTinv,1)!=size(rh1,1) || doClear
		try
			KKTinv = cholfact(G+Ai'*liy*Ai)
		catch E
			if isa(E,Base.LinAlg.PosDefException)
				KKTinv = pinv(G+Ai'*liy*Ai)
			else
				throw(E)
			end
		end
	end

	if isa(KKTinv,Array)
		if me==0
			dle = zeros(0)
			dx  = KKTinv*rh1
		else
			dle = (Ae*(KKTinv*Ae'))\(rh2-Ae*(KKTinv*rh1))
			dx  = KKTinv*(rh1+Ae'*dle)::Vector{T}
		end
	else
		if me==0
			dle =zeros(0)
			dx  = KKTinv\rh1
		else
			dle =  (Ae*(KKTinv\Ae'))\(rh2-Ae*(KKTinv\rh1))
			dx  = KKTinv\(rh1+Ae'*dle)::Vector{T}
		end
	end
	dli  = yLi\(-rpi-(rli./li)-Ai*dx)
	dy   = (-rli-y.*dli)./li
	return vec([dx;dy;dle;dli])::Vector{T}, KKTinv
end



"""
	function solveAugKKT(G,c,Ae,be,Ai,bi,x,y,le,li,rd,rpe,rpi,rli;KKTinv=[],doClear=true)

	builds and solves the augmented system. works so-so in practice
"""
function solveAugKKT{T}(G::AbstractArray{T},Ae::AbstractArray{T},Ai::AbstractArray{T},
	                  y::Vector{T},li::Vector{T},rd::Vector{T},rpe::Vector{T},rpi::Vector{T},rli::Vector{T};KKTinv=lufact(speye(1)),doClear=true)
	# build augmented matrix
    n  = size(G,1)
    me = size(Ae,1)
    mi = size(Ai,1)

    In  = speye(n)
    Imi = speye(mi)
    Znmi = spzeros(n,mi)
	Znme = spzeros(n,me)
	Znn = spzeros(n,n)
    Zmimi = spzeros(mi,mi)
	Zmeme = spzeros(me,me)
	Zmime = spzeros(mi,me)

    Y   = sdiag(y)
    Li  = sdiag(li)
	yLi = sdiag(y./li)

	rhs  = -[rd;rpe;rpi+(rli./li)]
	if size(KKTinv,1)!=size(rhs,1) || doClear
		# println("solveAugKKT: rebuild!")
    	KKT  =  [G -Ae' -Ai'; Ae  Zmeme Zmime'; Ai  Zmime yLi]
		KKTinv = lufact(KKT)
	end
	dt   = KKTinv\rhs::Vector{T}

	# solve for dy
	dx   = dt[1:n]
	dle  = dt[n+(1:me)]
	dli  = dt[n+me+1:end]
	dy   = (-rli-y.*dli)./li

	return vec([dx;dy;dle;dli])::Vector{T},KKTinv
end

"""
	function solveFullKKT(G,c,Ae,be,Ai,bi,x,y,le,li,rd,rpe,rpi,rli;KKTinv=[],doClear=true)

	builds and solves the full KKT system. Not recommended in practice
"""
function solveFullKKT{T}(G::AbstractArray{T},Ae::AbstractArray{T},Ai::AbstractArray{T},
	                  y::Vector{T},li::Vector{T},rd::Vector{T},rpe::Vector{T},rpi::Vector{T},rli::Vector{T};
	                  KKTinv=lufact(speye(1)),doClear=true)
	rhs = -[rd;rpe;rpi;rli]
	if size(KKTinv,1)!=size(rhs,1) || doClear
		# println("solveFull: rebuild!")
		KKT    = getJ(G,Ae,Ai,y,li)
		KKTinv = lufact(KKT)
	end

	return vec(KKTinv\rhs)::Vector{T}, KKTinv
end

"""
	function ipqp(G,c,Ae,be,Ai,bi,x,y,le,li;kwargs...)

	interior point solver for quadratic program based on Mehrotra's predictor-corrector approach

		min_x .5*x'*G*c + c'*x  subject to Ae*x = be and Ai*x >= bi

	Implementation is based on the presentation in Chapter 16 of

	Nocedal, J., & Wright, S. (2006). Numerical Optimization.
		Springer Science & Business Media. http://doi.org/10.1007/978-0-387-40065-5

	Optional keyword arguments:

		solveKKT          - function handle for linear solver. Options are currently
		                        solveNormalEqKKT (default)  (solves spsd linear system)
		                        solveAugKKT (symmetric indefinite solver)
		                        solveFullKKT
		maxIter::Int      - maximum number of iterations (default=10)
		ftol::Real        - stopping tolerance for norm of gradient of Lagrangian (default:1e-4)
		mutol::Real       - stopping tolerance on complementarity measure (default=1e-4)
		storeInterm::Bool - flag for storing intermediate iterates (default=false)
		out::Int          - controls verbosity
		                    (<1 --> no output, ==1 --> print final status, ==2 --> print status in each iteration)
		eta::Real         - reduction factor for steps (to remain interior) default=0.999

"""
function ipqp{T}(G::AbstractArray{T},c::Vector{T},Ae::AbstractArray{T},be::Vector{T},Ai::AbstractArray{T},bi::Vector{T},
	             x::Vector{T},y::Vector{T},le::Vector{T},li::Vector{T};maxIter::Int=10,ftol::Float64=1e-4,mutol::Float64=1e-4,
                 storeInterm::Bool=false,out::Int=1,eta::Float64=.999,solveKKT::Function=solveNormalEqKKT)

    n  = size(G,1)
    me = size(Ae,1)
    mi = size(Ai,1)

    his = zeros(maxIter+2,7)
    if storeInterm
        X  = zeros(n,maxIter)
        Y  = zeros(mi,maxIter)
        Le = zeros(me,maxIter)
        Li = zeros(mi,maxIter)
    end

	# 0) compute good starting guess (make sure y and li are sufficiently far in interior)
	mu  = dot(y,li)/mi
	rd,rpe,rpi,rli, = getF(G,c,Ae,be,Ai,bi,x,y,le,li)
	d,KKTinv  = solveKKT(G,Ae,Ai,y,li,rd,rpe,rpi,rli)
    (dx0,dy0,dle0,dli0) = splitVecIPQP(d,n,me,mi)
    y  = max.(1,abs.(y+dy0))
    li = max.(1,abs.(li+dli0))
	nrmrpe = (isempty(rpe)) ? 0.0 : norm(rpe)
   his[1,:] = [norm([rd;rpe;rpi;rli]) mu nrmrpe norm(rpi) norm(rd) norm(rli)  0]

    if out==2
        @printf "=== Interior Point for QP (n=%d,me=%d,mi=%d) ===\n" n me mi
        @printf "%4s\t%3s\t\t%2s\t\t%5s\t\t%5s\t\t%4s\t\t%4s\t\t%5s\n" "iter" "|F|" "mu" "|rpe|" "|rpi|" "|rd|" "|rl|" "alpha"
        @printf "%4d\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\n" -1 his[1,1] his[1,2] his[1,3] his[1,4] his[1,5]  his[1,6]
   end

    flag = -1; iter = 1;
    for iter=1:maxIter
		# get centrality and optimality conditions
		mu  = dot(y,li)/mi
		rd,rpe,rpi,rli, = getF(G,c,Ae,be,Ai,bi,x,y,le,li)

		# save and print history
        nrmrpe = (isempty(rpe))? 0.0 : norm(rpe)
		his[iter+1,1:6] = [norm([rd;rpe;rpi;rli]) mu nrmrpe norm(rpi) norm(rd) norm(rli)]
        if out==2
            @printf "%4d\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.2e" iter-1 his[iter,1] his[iter,2] his[iter,3] his[iter,4] his[iter,5]  his[iter,6]
        end

		# check stopping criteria
		if his[iter,1]<ftol && (his[iter,2]<mutol)
			(out==2) && @printf "\n"
            flag = 0
			iter -=1
            break
        end

		# 1) compute an affine scaling step (sigma=0)
		d,KKTinv           = solveKKT(G,Ae,Ai,y,li,rd,rpe,rpi,rli,KKTinv=KKTinv,doClear=true)
        (dxaff,dyaff,dleaff,dliaff) = splitVecIPQP(d,n,me,mi)
        alphaaff = stepLength(y,li,dyaff,dliaff;eta=eta)

 		# 2) compute centrality parameter
        muaff = dot(y+alphaaff*dyaff, li+alphaaff*dliaff)/mi
		sigma = (muaff/mu)^3

		# 3) solve for new step
		rd,rpe,rpi,rli, = getF(G,c,Ae,be,Ai,bi,x,y,le,li,dy=dyaff,dli=dliaff,tau=sigma*mu)
		d,KKTinv        = solveKKT(G,Ae,Ai,y,li,rd,rpe,rpi,rli,KKTinv=KKTinv,doClear=false)
        (dx,dy,dle,dli) = splitVecIPQP(d,n,me,mi)
        alpha = stepLength(y,li,dy,dli;eta=0.999)


        # 4) update
        x   += alpha*dx
        y   += alpha*dy
        le  += alpha*dle
        li  += alpha*dli

		his[iter,7] = alpha
		(out==2) && @printf "\t%1.2e\n" alpha

        if storeInterm
            X[:,iter]  = x
            Y[:,iter]  = y
            Le[:,iter] = le
            Li[:,iter] = li
        end
    end

    if out>=0
        if flag==-1
            @printf "ipqp iterated maxIter(=%d) times but reached only |F(x,s,u,v)|=%1.2e (ftol=%1.2e) and mu=%1.2e (mutol=%1.2e)\n" maxIter his[iter+1,1] ftol his[iter+1,2] mutol
        elseif flag==0 && out>=1
            @printf "ipqp achieved ftol=%1.2e and mutol=%1.2e at iteration %d. Returned result satisfies |F(x,y,s)|=%1.2e and mu=%1.2e.\n" ftol  mutol iter his[iter,1] his[iter,2]
        end
    end

    (x,y,le,li) = (storeInterm) ? (X[:,1:iter],Y[:,1:iter],Le[:,1:iter],Li[:,1:iter]) : (x,y,le,li)
    return x,y,le,li,his[1:iter+1,:],KKTinv
end
