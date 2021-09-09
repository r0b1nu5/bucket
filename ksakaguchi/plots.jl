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

function plot_cyc_vid(αs::Vector{Float64}, βmax::Vector{Float64}, βmin::Vector{Float64}, Θs0::Matrix{Float64}, α::Float64, β::Float64, θ0::Vector{Float64}, figid::String)
	n = length(θ0)
	θ = θ0 .- θ0[1]

	Θs = Θs0 - repeat(Θs0[[1,],:],n,1)
	Δs = mod.(Θs - [Θs[2:end,:];Θs[[1,],:]] .+ π,2π) .- π

	mz = 10.
	lw = 2.
	cmap = get_cmap("plasma")

	figure(figid,(16.,5.))
	
	subplot(1,3,1)
	PyPlot.fill([αs;αs[end:-1:1]],[βmax;βmin[end:-1:1]])
	PyPlot.plot(α,β,"o",color="C3")

	xlabel("α")
	ylabel("β")

	subplot(1,3,2)
	PyPlot.plot(sin.(2π*(0:n)/n),cos.(2π*(0:n)/n),"-k",linewidth=lw)
	for i in 1:n
		PyPlot.arrow(sin(2π*i/n),cos(2π*i/n),.25*cos(θ[i]),.25*sin(θ[i]),color="C0",width=.01)
		if i == 1
			PyPlot.plot(sin.(2π*i/n),cos.(2π*i/n),"s",color=cmap((i-1)/(n-1)),markersize=mz)
		elseif abs(i - 1 - n/3) < 1e-2
			PyPlot.plot(sin.(2π*i/n),cos.(2π*i/n),"v",color=cmap((i-1)/(n-1)),markersize=mz)
		else
			PyPlot.plot(sin.(2π*i/n),cos.(2π*i/n),"o",color=cmap((i-1)/(n-1)),markersize=mz)
		end
	end

	axis([-1.3,1.3,-1.3,1.3])
	xticks([])
	yticks([])

	subplot(1,3,3)
	PyPlot.plot([αs[1],αs[end]],[π/2,π/2],"--k")
	PyPlot.plot([αs[1],αs[end]],[-π/2,-π/2],"--k")
	for i in 1:n
		PyPlot.plot(αs,Δs[i,:],color=cmap((i-1)/(n-1)))
	end
	PyPlot.plot([α,α],[-π,π],"k")

	axis([αs[1],αs[end],-π,π])
	xlabel("α")
	ylabel("Δ_i")
end



