using PyPlot, Statistics

include("tools.jl")

function plot_traj_lin(ths::Array{Float64,2}, dhs::Array{Float64,2}, Cs::Array{Array{Int64,1},1}, h::Float64=.01, rot_frame::Bool=true)
	n,T = size(dhs)
	w = mean(dhs[:,end])
	
	figure()

	if rot_frame
		for i in 1:n
			subplot(3,1,1)
			PyPlot.plot(h*(0:T),ths[i,:] - w*h*(0:T))
			subplot(3,1,2)
			PyPlot.plot(h*(1:T),dhs[i,:] .- w)
		end
	else
		for i in 1:n
			subplot(3,1,1)
			PyPlot.plot(h*(0:T),ths[i,:])
			subplot(3,1,2)
			PyPlot.plot(h*(1:T),dhs[i,:])
		end
	end
	
	qs = Array{Int64,2}(undef,0,T+1)
	for C in Cs
		qs = [qs;[winding(ths[:,t],C) for t in 1:T+1]']
		subplot(3,1,3)
		PyPlot.plot(h*(0:T),qs[end,:])
	end
	xlabel("t")
	ylabel("q")

	subplot(3,1,1)
	ylabel("θ")

	subplot(3,1,2)
	ylabel("ω")
end


