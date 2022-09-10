using LinearAlgebra, PyPlot 

include("L2B.jl")

function analyze_sol(θ::Array{Float64,1}, L::Array{Float64,2}, ω::Array{Float64,1}, α::Float64, plt::Bool=false, r::Float64=1.)
	B,w = L2B(L)
	W = diagm(0 => w)
	n,m = size(B)
	B1 = B.*(B .> 0)
	B2 = -B.*(B .< 0)
	B12 = [B1 B2]
	BB = [B -B]
	WW = [W zeros(m,m);zeros(m,m) W]


	dθ = ω - B12*WW*sin.(BB'*θ .- α)

	if plt
		figure("analysis")
		subplot(1,2,1)
		PyPlot.plot(r*cos.(θ),r*sin.(θ),"o")
		subplot(1,2,2)
		PyPlot.plot(r*ones(n),dθ,"o")
	end

	return dθ, maximum(dθ), minimum(dθ)
end

function analyze_sol(θ::Array{Float64,2}, L::Array{Float64,2}, ω::Array{Float64,1}, α::Float64, plt::Bool=false)
	n,s = size(θ)
	
	dθs = Array{Float64,2}(undef,n,0)
	mas = Array{Float64,1}()
	mis = Array{Float64,1}()

	if s > 1
		rs = LinRange(1.,2.,s)
	else
		rs = [1.,]
	end

	for i in 1:s
		xxx = analyze_sol(θ[:,i],L,ω,α,plt,rs[i])

		dθs = [dθs xxx[1]]
		push!(mas,xxx[2])
		push!(mis,xxx[3])
	end

	return dθs, mas, mis
end



