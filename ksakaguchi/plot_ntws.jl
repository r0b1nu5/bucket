using PyPlot

function plot_ntws(ntw::String, ms::Float64=10., lw::Float64=2.)
	if ntw == "cyc1_18"
		n = 18
		t = 2π*(0:n)/n

		PyPlot.plot(sin.(t),cos.(t),"ok",markersize=ms)
		PyPlot.plot(sin.(t),cos.(t),"-k",linewidth=lw)
		PyPlot.plot(sin(t[1]),cos(t[1]),"o",color="C2",markersize=1.1*ms)
		PyPlot.plot(sin(t[7]),cos(t[7]),"o",color="C3",markersize=1.1*ms)
	elseif ntw == "cyc2_18"
		n = 18
		t = 2π*(0:n)/n

		PyPlot.plot(sin.(t),cos.(t),"ok",markersize=ms)
		PyPlot.plot(sin.(t),cos.(t),"-k",linewidth=lw)
		PyPlot.plot(sin.(t[[1,10]]),cos.(t[[1,10]]),"-k",linewidth=lw)
		PyPlot.plot(sin(t[1]),cos(t[1]),"o",color="C2",markersize=1.1*ms)
		PyPlot.plot(sin(t[10]),cos(t[10]),"o",color="C3",markersize=1.1*ms)
	elseif ntw == "cyc2_12"
		n = 12
		t = 2π*(0:n)/n

		PyPlot.plot(sin.(t),cos.(t),"ok",markersize=ms)
		PyPlot.plot(sin.(t),cos.(t),"-k",linewidth=lw)
		PyPlot.plot(sin.(t[[1,7]]),cos.(t[[1,7]]),"-k",linewidth=lw)
		PyPlot.plot(sin(t[1]),cos(t[1]),"o",color="C2",markersize=1.1*ms)
		PyPlot.plot(sin(t[7]),cos(t[7]),"o",color="C3",markersize=1.1*ms)

	else
		@info "No network found..."
	end

	xticks([])
	yticks([])
	
	return nothing
end



