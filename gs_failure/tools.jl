using LinearAlgebra

function jacobian(th::Array{Float64,1}, LL::Array{Float64,2}=zeros(1,1))
	n = length(th)

	if size(LL) == (1,1)
		L = n*diagm(0 => ones(n)) - ones(n,n)
		A = 1 .- diagm(0 => ones(n))
	else
		L = LL
		A = -L.*(1 .- diagm(0 => ones(n)))
	end

	J0 = A.*cos.(th*ones(1,n) - ones(n)*th')
	J = J0 - diagm(0 => vec(sum(J0,dims=2)))

	return J
end



function span_tree(B::Matrix{Float64})
	n,m = size(B)

	ids = [1,]
	T = B[:,ids]
	i = 1

	while size(T)[2] < n-1
		i += 1
		if rank([T B[:,i]]) == size(T)[2]+1
			T = [T B[:,i]]
			push!(ids,i)
		end
	end

	return T,ids
end


function dcc(x::Float64)
	return mod(x + π,2π) - π
end

function dcc(x::Vector{Float64})
	return [dcc(x[i]) for i in 1:length(x)]
end

function dcc(x::Matrix{Float64})
	n,m = size(x)
	d = zeros(n,m)
	for i in 1:m
		d[:,i] = dcc(x[:,i])
	end
	return d
end


