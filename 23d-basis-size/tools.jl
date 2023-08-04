using LinearAlgebra, Statistics, SpecialFunctions, Distributions

function L2B(L::Matrix{Float64})
	n = size(L)[1]

	Bout = zeros(n,0)
	Bin = zeros(n,0)
	w = Float64[]

	for i in 1:n
		for j in 1:n
			if i !== j && abs(L[i,j]) > 1e-6
				Bout = [Bout I[1:n,i]]
				Bin = [Bin I[1:n,j]]
				push!(w,-L[i,j])
			end
		end
	end

	return Bout,Bin,w
end

function winding(θ::Vector{Float64}, σ::Vector{Int64})
	Δ = θ[[σ[2:end];σ[1]]] - θ[σ]
	return round(Int64,sum(mod.(Δ .+ π,2π) .- π)/(2π))
end

function gen_cycle_undir(n::Int64)
	return diagm(0 => 2*ones(n)) - diagm(1 => ones(n-1)) - diagm(-1 => ones(n-1)) - diagm(n-1 => ones(1)) - diagm(1-n => ones(1))
end

function gen_cycle_dir(n::Int64)
	return diagm(0 => ones(n)) - diagm(1 => ones(n-1)) - diagm(1-n => ones(1))
end

function gof_normal(x::Vector{<:Any}, qmax::Int64, iter::Int64=5000)
	μ0 = mean(x)
	σ0 = std(x)
	n = length(x)

	d0 = ks_discr_normal(x,qmax,μ0,σ0)

	c = 0 

	for i in 1:iter
		y = round.(rand(Normal(μ0,σ0),n))
		μ1 = mean(y)
		σ1 = std(y)

		d1 = ks_discr_normal(y,qmax,μ1,σ1)
		
		if d1 > d0
			c +=1
#			@info "($d0,$d1) <-------------"
		else
#			@info "($d0,$d1)"
		end
	end

	return c/iter
end

function ks_discr_normal(x::Vector{<:Any}, qmax::Int64, μ::Float64=0., σ::Float64=1.)
	cdft = [(.5+.5*erf((q+.5-μ)/(sqrt(2)*σ))) for q in -qmax:qmax]
	cdfe = [sum(x .<= q)/length(x) for q in -qmax:qmax]

	return norm(cdft-cdfe,Inf)
end

function qqplot(x::Vector{<:Any}, qmax::Int64, μ::Float64=0., σ::Float64=1.; color::String)
	cdft = [(.5+.5*erf((q+.5-μ)/(sqrt(2)*σ))) for q in -qmax:qmax]
	cdfe = [sum(x .<= q)/length(x) for q in -qmax:qmax]
	
	figure("Q-Q plot")
	PyPlot.plot(sort(cdft),sort(cdfe),".",color=color)
	xlabel("Theoretical")	
	ylabel("Empirical")
	
end




