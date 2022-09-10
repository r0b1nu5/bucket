using PyPlot

function plot_comparison(ths::Array{Float64,2}, L::Array{Float64,2}, h::Float64=.01)
	n,T = size(ths)

	figure()
	for i in 1:n
		subplot(2,1,1)
		PyPlot.plot(h*(1:T),ths[i,:] .- ths[i,1])
		subplot(2,1,2)
		PyPlot.plot(h*(1:T),(L*ths)[i,:] .- (L*ths[:,1])[i])
	end

	subplot(2,1,1)
	ylabel("x(t) - x*")

	subplot(2,1,2)
	xlabel("t")
	ylabel("ψ(t) - ψ*")
end


