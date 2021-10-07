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



