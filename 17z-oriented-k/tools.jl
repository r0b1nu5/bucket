using LinearAlgebra

function loc_r(th::Array{Float64,1}, A::Array{Float64,2})
	n = length(th)

	d = vec(sum(A,dims=2))
	z = exp.(im*A*x)./d
	
	r = norm.(z)
	phi = angle.(z)

	return z,r,phi
end


function L2B_dir(L::Array{Float64,2})
	n = size(L)[1]
	
	S = Array{Float64,2}(undef,n,0) # "source matrix", S_{ie}=1 iff i is the source of e.
	T = Array{Float64,2}(undef,n,0) # "target matrix", T_{ie}=1 iff i is the target of e.
	w = Array{Float64,1}(undef,0)
	for i in 1:n
		for j in 1:n
			if (L[i,j] != 0.0) && (i != j)
				s = zeros(n)
				t = zeros(n)
				s[i] = 1.0
				t[j] = 1.0
				S = [S s]
				T = [T t]
				push!(w,-L[i,j])
			end
		end
	end

	B = S - T
	
	return B,S,T,w
end
			


