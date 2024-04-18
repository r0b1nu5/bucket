using DelimitedFiles

function lv(x0::Vector{Float64}, r::Vector{Float64}, θ::Vector{Float64}, B::Matrix{Float64}, h::Float64=.01, save::Bool=false, tol::Float64=1e-6, maxiter::Int64=10000)
	x = copy(x0)
	xs = copy(x0)

	corr = 1e10
	iter = 0

	while corr > tol && iter < maxiter
		iter += 1
		if iter%100 == 0
			@info "iter: $iter, log-err: $(round(log(corr))))"
		end

		k1 = f_lv(x,r,θ,B)
		k2 = f_lv(x + h*k1/2,r,θ,B)
		k3 = f_lv(x + h*k2/2,r,θ,B)
		k4 = f_lv(x + h*k3,r,θ,B)

		dx = (k1+2*k2+2*k3+k4)/6
		corr = maximum(abs.(dx))

		x += h*dx
		if save 
			xs = [xs x]
		end
	end

	if save 
		return xs
	else
		return x
	end
end


function f_lv(x::Vector{Float64}, r::Vector{Float64}, θ::Vector{Float64}, B::Matrix{Float64})
	return x.*(r - θ.*x + B*x)
end

# From "notes-pj-240417.pdf"
function lv_bunin(N0::Vector{Float64}, A::Matrix{Float64}, κ::Vector{Float64}, μ::Float64, σ::Float64, tol::Float64=1e-6, min_iter::Int64=1000, max_iter::Int64=10000, h::Float64=.01)
	S = length(N0)

	N = N0
	Ns = N0
	iter = 1
	c = 0
	err = 1000.

	while iter < min_iter
		iter += 1

		k1 = f_lv_bunin(N,A,κ,μ,σ)
		k2 = f_lv_bunin(N + h/2*k1,A,κ,μ,σ)
		k3 = f_lv_bunin(N + h/2*k2,A,κ,μ,σ)
		k4 = f_lv_bunin(N + h*k3,A,κ,μ,σ)

		dN = (k1+2*k2+2*k3+k4)/6
		N += h*dN
		N .*= (N .> 1e-15)
		Ns = [Ns N]

		if (iter-1)%1000 == 0
			@info "iter: $iter"
			c += 1
			writedlm("data/Ns-$c.csv",Ns[:,1:end-1],',')
			Ns = N
		end
	end

	Ss = sum(N .> 1e-15)
	S0 = length(N)+1

	while Ss != S0 && iter < max_iter
		S0 = sum(N .> 1e-15)

		for i in 1:min_iter
			iter += 1

			k1 = f_lv_bunin(N,A,κ,μ,σ)
			k2 = f_lv_bunin(N + h/2*k1,A,κ,μ,σ)
			k3 = f_lv_bunin(N + h/2*k2,A,κ,μ,σ)
			k4 = f_lv_bunin(N + h*k3,A,κ,μ,σ)
	
			dN = (k1+2*k2+2*k3+k4)/6
			N += h*dN
			N .*= (N .> 1e-15)
			Ns = [Ns N]
	
			if (iter-1)%1000 == 0
				@info "iter: $iter"
				c += 1
				writedlm("data/Ns-$c.csv",Ns[:,1:end-1],',')
				Ns = N
			end
		end

		Ss = sum(N .> 1e-15)
	end

	c += 1
	writedlm("data/Ns-$c.csv",Ns,',')

	Ns = zeros(length(N),0)
	for i in 1:c
		Ns = [Ns readdlm("data/Ns-$i.csv",',')]
		rm("data/Ns-$i.csv")
	end

	return Ns
end

function lv_bunin(N0::Vector{Float64}, A::Matrix{Float64}, κ::Float64, μ::Float64, σ::Float64, tol::Float64=1e-6, min_iter::Int64=1000, max_iter::Int64=10000, h::Float64=.01)
	return lv_bunin(N0,A,κ*ones(length(N0)),μ,σ,tol,min_iter,max_iter,h)
end

function f_lv_bunin(N::Vector{Float64}, A::Matrix{Float64}, κ::Vector{Float64}, μ::Float64, σ::Float64)
	S = length(N)

	return N.*(κ - N .- μ/S*sum(N) - σ/sqrt(S)*A*N)
end




