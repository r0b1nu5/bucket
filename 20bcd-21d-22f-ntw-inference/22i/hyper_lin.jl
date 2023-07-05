using LinearAlgebra, DelimitedFiles, Distributions

include("cnoise.jl")

function hyper_lin(M2::Array{Float64,2}, M3::Array{Float64,3}, x0::Vector{Float64}, τ0::Float64, ξ0::Float64, h::Float64=.01, max_iter::Int64=10000, tol::Float64=1e-6)
	n = length(x0)

	x = x0
	xs = x0
	dxs = zeros(n,0)

	iter = 0
	c = 0
	err = 1000.

	ξ = rand(Normal(0.,1.),n)

	while iter < max_iter && err > tol
		iter += 1
		if iter%1000 == 0
			@info "iter: $iter"
			c += 1
			writedlm("temp/xs_$c.csv",xs[:,1:end-1],',')
			xs = xs[:,end]
			writedlm("temp/dxs_$c.csv",dxs[:,1:end],',')
			dxs = zeros(n,0)
		end

		ξ = [cnoise(ξ[i],τ0/h) for i in 1:n]

		k1 = M2*x + prod3(M3,x) + ξ0*ξ
		k2 = M2*(x+h/2*k1) + prod3(M3,x+h/2*k1) + ξ0*ξ
		k3 = M2*(x+h/2*k2) + prod3(M3,x+h/2*k2) + ξ0*ξ
		k4 = M2*(x+h*k3) + prod3(M3,x+h*k3) + ξ0*ξ

		dx = (k1 + 2*k2 + 2*k3 + k4)/6
		x += h*dx

		xs = [xs x]
		dxs = [dxs dx]

		err = maximum(abs.(dx))
	end

	Xs = zeros(n,0)
	dXs = zeros(n,0)

	for i in 1:c
		Xs = [Xs readdlm("temp/xs_$i.csv",',')]
		rm("temp/xs_$i.csv")
		dXs = [dXs readdlm("temp/dxs_$i.csv",',')]
		rm("temp/dxs_$i.csv")
	end
	Xs = [Xs xs]
	dXs = [dXs dxs]

	return Xs,dXs
end


function prod3(M::Array{Float64,3}, x::Vector{Float64})
	n = length(x)
	p = Float64[]
	for i in 1:n
		t = 0.
		for j in 1:n
			for k in 1:n
				t += M[i,j,k]*x[j]*x[k]
			end
		end
		push!(p,t)
	end

	return p
end

# generate a symetric tensor of order 3.
function gen_sym_tensor(n::Int64, tmin::Float64=0., tmax::Float64=1.)
	T = zeros(n,n,n)
	for i in 1:n
		for j in i:n
			for k in j:n
				x = (tmax-tmin)*rand() + tmin
				T[i,j,k] = x
				T[i,k,j] = x
				T[j,i,k] = x
				T[j,k,i] = x
				T[k,i,j] = x
				T[k,j,i] = x
			end
		end
	end
	#=
	for i in 1:n-1
		for j in i+1:n
			x = (tmax-tmin)*rand() + tmin
			T[i,i,j] = x
			T[i,j,i] = x
			T[j,i,i] = x
		end
	end
	for i in 1:n
		x = (tmax-tmin)*rand() + tmin
		T[i,i,i] = x
	end
=#
	return T
end


# generate a symetric tensor of order 3.
# Diagonal has opposite sign as the off-diagonals.
function gen_sym_tensor2(n::Int64, tmin::Float64=0., tmax::Float64=1.)
	T = zeros(n,n,n)
	for i in 1:n
		for j in i:n
			for k in j:n
				x = (tmax-tmin)*rand() + tmin
				T[i,j,k] = x
				T[i,k,j] = x
				T[j,i,k] = x
				T[j,k,i] = x
				T[k,i,j] = x
				T[k,j,i] = x
			end
		end
	end

	for i in 1:n
		x = -((tmax-tmin)*rand() + tmin)
		T[i,i,i] = x
	end

	#=
	for i in 1:n-1
		for j in i+1:n
			x = (tmax-tmin)*rand() + tmin
			T[i,i,j] = x
			T[i,j,i] = x
			T[j,i,i] = x
		end
	end
	for i in 1:n
		x = (tmax-tmin)*rand() + tmin
		T[i,i,i] = x
	end
=#
	return T
end


# generate a symetric tensor of order 3.
# Diagonal has opposite sign compared to and larger magnitude than off-diagonals.
function gen_sym_tensor3(n::Int64, tmin::Float64=0., tmax::Float64=1.)
	T = zeros(n,n,n)
	for i in 1:n
		for j in i:n
			for k in j:n
				x = (tmax-tmin)*rand() + tmin
				T[i,j,k] = x
				T[i,k,j] = x
				T[j,i,k] = x
				T[j,k,i] = x
				T[k,i,j] = x
				T[k,j,i] = x
			end
		end
	end

	for i in 1:n
		x = -100*((tmax-tmin)*rand() + tmin)
		T[i,i,i] = x
	end

	#=
	for i in 1:n-1
		for j in i+1:n
			x = (tmax-tmin)*rand() + tmin
			T[i,i,j] = x
			T[i,j,i] = x
			T[j,i,i] = x
		end
	end
	for i in 1:n
		x = (tmax-tmin)*rand() + tmin
		T[i,i,i] = x
	end
=#
	return T
end



