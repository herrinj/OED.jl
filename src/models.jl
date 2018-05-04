export parallelTomoMod

"""
function parallelTomoMod(w,n)

model for tomography problems where for each angle n parallel 
projections are used.

Input:
    
    w - weights for each angle
    n - number of projections per angle

Output:
    sig,dsig - model and derivative 

"""
function parallelTomoMod(w::Vector{Float64},n::Int)
    sig  = kron(w,ones(n))
    dsig = kron(speye(length(w)),ones(n))
    return sig,dsig
end