using PyPlot

include("kuramoto.jl")
include("splay_states.jl")

# Runs the reversed KM from the initial state θ0, until it reaches the neighborhood of a splay state. Then jumps to the symmetric point wrt the splay state and runs the standard KM from there. 

function switch_dyn(θ0::Vector{Float64}, n_switches::Int64, L::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, ω::Vector{Float64}, plot::Bool=false, h::Float64=.01, tol::Float64=1e-4)
	n = length(θ0)
	δ = 1e-4

	θ = copy(θ0)

	c = 0

	while c < n_switches
		c += 1
		θs = kuramoto_series((-1)^c*L,(-1)^c*ω,θ,h,10000,tol)
		v = θs[:,end] - θs[:,end-1]
		θ = θs[:,end] + δ*v
		writedlm("temp_data/ths_$c.csv",θs,',')
	end

	Θs = Matrix{Float64}(undef,n,0)
	for d in 1:c
		Θs = [Θs readdlm("temp_data/ths_$d.csv",',')]
		rm("temp_data/ths_$d.csv")
	end

	dΘs = [Θs[2:n,:];Θs[[1,],:]] - Θs
	dmax = [maximum(abs.(dΘs[:,i])) for i in 1:size(Θs)[2]]

	qs = [winding(Θs[:,i],Vector(1:n)) for i in 1:size(Θs)[2]]

	if plot
		subplot(1,3,1)
		for i in 1:n
			PyPlot.plot(Θs[i,:])
		end
		
		subplot(1,3,2)
		PyPlot.plot(dmax)

		subplot(1,3,3)
		PyPlot.plot(qs)
	end

	return Θs
end






