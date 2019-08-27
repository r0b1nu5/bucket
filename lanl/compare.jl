using FFTW, PyPlot, Statistics

include("find_trig.jl")


function compare(Xs::Array{Float64,2}, dt::Float64)
	nn,T = size(Xs)
	n = Int(nn/2)

	fs = (-T/2:T/2-1)./(dt*T)
	
	ths = vec(sum(Xs[1:n,:],dims=1))
	oms = vec(sum(Xs[n+1:2*n,:],dims=1))

	a1,f1,p1,sh1 = find_sin(ths,dt)
	a2,f2,p2,sh2 = find_cos(oms,dt)
	
	Fths = norm.(fft(ths))
	Foms = norm.(fft(oms))
	
	figure()
	subplot(121)
	PyPlot.plot(fs,fftshift(Fths),"-",color="C0")
	PyPlot.plot([f1,f1],[1e-8,1.1*maximum(Fths)],"o",color="C0")
	PyPlot.plot([f2,f2],[1e-8,1.1*maximum(Fths)],"--",color="C1")
	axis([-.02,.02,0,1.1*maximum(Fths)])
	subplot(122)
	PyPlot.plot(fs,fftshift(Fths),"-",color="C0")
	PyPlot.plot([f1,f1],[1e-8,1.1*maximum(Fths)],"o",color="C0")
	PyPlot.plot([f2,f2],[1e-8,1.1*maximum(Fths)],"--",color="C1")

	figure()
	subplot(121)
	PyPlot.plot(fs,fftshift(Foms),"-",color="C1")
	PyPlot.plot([f2,f2],[1e-8,1.1*maximum(Foms)],"o",color="C1")
	PyPlot.plot([f1,f1],[1e-8,1.1*maximum(Foms)],"--",color="C0")
	axis([-.02,.02,0,1.1*maximum(Foms)])
	subplot(122)
	PyPlot.plot(fs,fftshift(Foms),"-",color="C1")
	PyPlot.plot([f2,f2],[1e-8,1.1*maximum(Foms)],"o",color="C1")
	PyPlot.plot([f1,f1],[1e-8,1.1*maximum(Foms)],"--",color="C0")

end

