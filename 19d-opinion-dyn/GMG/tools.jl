using LinearAlgebra

# n agents, p parties, uniform initial distribution, no re-ordering.
function gen_opinion_1(n::Int64, p::Int64)
	X = rand(n,p)
	S = sum(X,dims=2)
	X ./= repeat(S,1,p)

	return X
end

# n agents, p parties, each opinion is re-drawn if it does not satisfy the ordering
function gen_opinion_2(n::Int64, p::Int64)
	X = zeros(0,p)
	for k in 1:n
		x = rand(p)
		j = findmax(x)[2]
		if j == 1
			test = prod([x[i-1] > x[i] for i in j+1:p])
		elseif j == p
			test = prod([x[i] < x[i+1] for i in 1:j-1])
		else
			test = prod([x[i] < x[i+1] for i in 1:j-1])*prod([x[i-1] > x[i] for i in j+1:p])
		end
		
		while !test
			x = rand(p)
			j = findmax(x)[2]
			if j == 1
				test = prod([x[i-1] > x[i] for i in j+1:p])
			elseif j == p
				test = prod([x[i] < x[i+1] for i in 1:j-1])
			else
				test = prod([x[i] < x[i+1] for i in 1:j-1])*prod([x[i-1] > x[i] for i in j+1:p])
			end
		end

		x ./= sum(x)
		X = [X;x']
	end

	return X
end

# n agents, each party has a number of partisan equal to props[i]*n, ordering is satisfied.
function gen_opinion_3(n::Int64, props::Vector{Float64})
	p = length(props)
	prop = props/sum(props)
	ns = round.(Int64,prop*n)
	
	X = zeros(0,p)
	j = 1
	for k in 1:ns[j]
		x = rand(p)
		test = prod([x[i-1] > x[i] for i in j+1:p])
		while !test
			x = rand(p)
			test = prod([x[i-1] > x[i] for i in j+1:p])
		end
		x ./= sum(x)
		X = [X;x']
	end

	for j in 2:p-1
		for k in 1:ns[j]
			x = rand(p)
			test = prod([x[i] < x[i+1] for i in 1:j-1])*prod([x[i-1] > x[i] for i in j+1:p])

			while !test
				x = rand(p)
				test = prod([x[i] < x[i+1] for i in 1:j-1])*prod([x[i-1] > x[i] for i in j+1:p])
			end

			x ./= sum(x)
			X = [X;x']
		end
	end

	j = p
	for k in 1:n-sum(ns[1:p-1])
		x = rand(p)
		test = prod([x[i] < x[i+1] for i in 1:j-1])
		while !test
			x = rand(p)
			test = prod([x[i] < x[i+1] for i in 1:j-1])
		end
		x ./= sum(x)
		X = [X;x']
	end

	return X
end


function order_opinion(X::Matrix{Float64})
	n,p = size(X)

	for i in 1:n
		m,j = findmax(X[i,:])
		x = [sort(X[i,1:j-1]);X[i,j];sort(X[i,j+1:p],rev=true)]
		X[i,:] = x
	end

	return X
end

function Lϵ(x0::Array{Float64,1}, ϵ::Float64)
	n = length(x0)

	A = Float64.((0 .< abs.(repeat(x0,1,n) - repeat(x0',n,1)) .< ϵ))
	D = diagm(0 => vec(sum(A,dims=1)))
	L = D - A

	return L, A, D
end

function Lϵ(X::Matrix{Float64}, ϵ::Float64, q::Int64=1)
	n,p = size(X)

	x = [X[i,:] for i in 1:n]

	A = zeros(n,n)
	for i in 1:n-1
		for j in i+1:n
			A[i,j] = Float64(norm(x[i]-x[j],q) < ϵ)
			A[j,i] = A[i,j]
		end
	end

	d = vec(sum(A,dims=1))
	D = diagm(0 => d)
	L = D - A
	
	Di = pinv(D)
	M = inv(diagm(0 => ones(n)) + Di*L)

	return L, A, D, M
end

function vote(X::Matrix{Float64})
	return Float64.(X .>= repeat(maximum(X,dims=2),1,p))
end

function party_vote(X::Matrix{Float64})
	n,p = size(X)

	return [findmax(X[i,:])[2] for i in 1:n]
end

function election_result(X::Matrix{Float64})
	n,p = size(X)

	V = vote(X)

	return vec(sum(V,dims=1))./n
end

