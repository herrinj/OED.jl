export linearInter

function linearInter(T::Array{Float64},domain::Array{Float64},xin::Array{Float64};doDerivative::Bool=false)

	dim = round(Int,length(domain)/2)
	n   = (dim==1) ? length(T) : size(T)
	h   = sdiag(1./collect(n)) * vec(domain[2:2:end]-domain[1:2:end])
	np  = round(Int,length(xin)/dim)
	x   = reshape(copy(xin),np,dim)


	# map x from [h/2,domain-h/2] -> [1,m],
	for i=1:dim
		x[:,i] = (x[:,i]-domain[2*i-1])/h[i] + 0.5;
	end
	Tc = zeros(np)                                # initialize output
	Valid(j) = (0.<x[:,j]) .& (x[:,j].<(n[j]+1))      # determine indices of valid points

	if dim==1
		valid = find(Valid(1))
	elseif dim==2
		valid = find( Valid(1) .& Valid(2) );
	elseif dim==3
		valid = find( Valid(1) .& Valid(2) .& Valid(3) );
	end

	if isempty(valid)
		if !doDerivative
			return Tc
		else
			return Tc, spzeros(np,dim*np)
		end
	end

	pad = 1; TP = zeros(tuple(collect(n)+2*pad...));                 # pad data to reduce cases
	P     = round.(Int,floor.(x)); x = x-P;                        # split x into integer/remainder
	p(j)  = P[valid,j]
	xi(j) = x[valid,j]

	# increments for linearized ordering
	i1 = 1; i2 = size(T,1)+2*pad; i3 = (size(T,1)+2*pad)*(size(T,2)+2*pad);

	if dim==1
		TP[pad+(1:n)] = vec(T)
		p = pad + p(1);
	    Tc[valid] = TP[p].* (1-xi(1)) + TP[p+1].*xi(1);   # compute weighted sum

	    if doDerivative
			dT = zeros(np)
			dT[valid] = TP[p+1]-TP[p]
			dT = sdiag(dT[:,1]/h[1])
		end;
	elseif dim==2
	    TP[pad+(1:n[1]),pad+(1:n[2])] = T
	    TP = vec(TP)
	    p  = (pad + p(1)) + i2*(pad + p(2) - 1)


	    # compute Tc as weighted sum
	    Tc[valid] = (TP[p]   .* (1-xi(1))  + TP[p+i1]    .*xi(1)) .* (1-xi(2))  + (TP[p+i2] .* (1-xi(1)) + TP[p+i1+i2] .*xi(1)) .* xi(2)

		if doDerivative
			dT = zeros(np,dim)
			dT[valid,1] = (TP[p+i1]-TP[p]).*(1-xi(2)) + (TP[p+i1+i2]-TP[p+i2]).*xi(2)
			dT[valid,2] = (TP[p+i2]-TP[p]).*(1-xi(1)) + (TP[p+i1+i2]-TP[p+i1]).*xi(1)
			dT = spdiagm((dT[:,1]/h[1], dT[:,2]/h[2]),(0,np), np, 2*np)
			# dT = [sdiag(dT[:,1]/h[1])  sdiag(dT[:,2]/h[2])]
		end
	  elseif dim==3
	    TP[pad+(1:n[1]),pad+(1:n[2]),pad+(1:n[3])] = T
	    p  = (pad + p(1)) + i2*(pad + p(2) - 1) + i3*(pad + p(3) -1);
	    # compute Tc as weighted sum
	    Tc[valid] = ((TP[p].*(1-xi(1))+TP[p+i1].*xi(1)).*(1-xi(2)) +
		  (TP[p+i2].*(1-xi(1))+TP[p+i1+i2].*xi(1)).*(xi(2))).*(1-xi(3)) +
	      +((TP[p+i3].*(1-xi(1))+TP[p+i1+i3].*xi(1)).*(1-xi(2)) +
	      +(TP[p+i2+i3].*(1-xi(1))+TP[p+i1+i2+i3].*xi(1)).*(xi(2))).*(xi(3))

	    if doDerivative
			dT = zeros(np,dim)
			dT[valid,1] = ((TP[p+i1]-TP[p]).*(1-xi(2))+(TP[p+i1+i2]-TP[p+i2]).*xi(2)).*(1-xi(3)) +
	                      ((TP[p+i1+i3]-TP[p+i3]).*(1-xi(2))+(TP[p+i1+i2+i3]-TP[p+i2+i3]).*xi(2)).*(xi(3))
		    dT[valid,2] = ((TP[p+i2]-TP[p]).*(1-xi(1))+(TP[p+i1+i2]-TP[p+i1]).*xi(1)).*(1-xi(3)) +
	                      ((TP[p+i2+i3]-TP[p+i3]).*(1-xi(1))+(TP[p+i1+i2+i3]-TP[p+i1+i3]).*xi(1)).*(xi(3))
		    dT[valid,3] = ((TP[p+i3].*(1-xi(1))+TP[p+i1+i3].*xi(1)).*(1-xi(2)) +
	                      (TP[p+i2+i3].*(1-xi(1))+TP[p+i1+i2+i3].*xi(1)).*(xi(2))) -
	                      ((TP[p].*(1-xi(1))+TP[p+i1].*xi(1)).*(1-xi(2)) +
	                      (TP[p+i2].*(1-xi(1))+TP[p+i1+i2]*xi(1)).*(xi(2)));
			dT = [sdiag(dT[:,1]/h[1])  sdiag(dT[:,2]/h[2])  sdiag(dT[:,3]/h[3])]
		end
	end
	if doDerivative
		return Tc,dT
	else
		return Tc
	end
end
