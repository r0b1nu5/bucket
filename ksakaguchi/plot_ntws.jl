using PyPlot

function plot_ntws(ntw::String, ω::Vector{Float64}, ms::Float64=10., lw::Float64=2., name_cmap::String="plasma")
	cmap = get_cmap(name_cmap)

	ωmax = maximum(ω)
	ωmin = minimum(ω)
	ωcol = (ω .- ωmin)./(ωmax - ωmin)

	if ntw == "cyc1_18"
		n = 18
		t = 2π*(0:n)/n .- π/3

		PyPlot.plot(sin.(t),cos.(t),"-k",linewidth=lw)
		for i in 1:n
			PyPlot.plot(sin(t[i]),cos(t[i]),"o",color=cmap(ωcol[i]),markersize=ms)
		end
	elseif ntw == "cyc2_18"
		n = 18
		t = 2π*(0:n)/n .- π/4

		PyPlot.plot(sin.(t),cos.(t),"-k",linewidth=lw)
		PyPlot.plot(sin.(t[[1,10]]),cos.(t[[1,10]]),"-k",linewidth=lw)
		for i in 1:n
			PyPlot.plot(sin(t[i]),cos(t[i]),"o",color=cmap(ωcol[i]),markersize=ms)
		end
	elseif ntw == "cyc2_12"
		n = 12
		t = 2π*(0:n)/n

		PyPlot.plot(sin.(t),cos.(t),"-k",linewidth=lw)
		PyPlot.plot(sin.(t[[1,7]]),cos.(t[[1,7]]),"-k",linewidth=lw)	
		for i in 1:n
			PyPlot.plot(sin(t[i]),cos(t[i]),"o",color=cmap(ωcol[i]),markersize=ms)
		end
	elseif ntw == "cyc1_12"
		n = 12
		t = 2π*(0:n)/n

		PyPlot.plot(sin.(t),cos.(t),"-k",linewidth=lw)
		for i in 1:n
			PyPlot.plot(sin(t[i]),cos(t[i]),"o",color=cmap(ωcol[i]),markersize=ms)
		end
	else
		@info "No network found..."
	end

	xticks([])
	yticks([])
	
	return nothing
end

function plot_ntws(ntw::String, ω::Vector{Float64}, θ::Vector{Float64}, ms::Float64=10., lw::Float64=2., name_cmap::String="plasma")
	if ntw == "cyc1_18"
		n = 18
		t = 2π*(0:n-1)/n .- π/3

		for i in 1:n
			PyPlot.plot(sin(t[i]) .+ [0.,.3*cos(θ[i])],cos(t[i]) .+ [0.,.3*sin(θ[i])],"k")
		end
	elseif ntw == "cyc2_18"
		n = 18
		t = 2π*(0:n-1)/n .- π/4
		
		for i in 1:n
			PyPlot.plot(sin(t[i]) .+ [0.,.3*cos(θ[i])],cos(t[i]) .+ [0.,.3*sin(θ[i])],"k")
		end
	end

	plot_ntws(ntw,ω,ms,lw,name_cmap)
end



