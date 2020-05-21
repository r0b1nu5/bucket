using PyPlot

function plot_circ_lin(ths::Array{Float64,2})
	n,T = size(ths)

	figure("lin",(8,8))
	ts = LinRange(0,2pi,200)
	tt = LinRange(0,T,5)
	for t in tt[2:end]
		PyPlot.plot(t*cos.(ts), t*sin.(ts), ":k")
	end
	PyPlot.plot([0,0], [-T,T], ":k")
	PyPlot.plot([-T,T], [0,0], ":k")
	PyPlot.plot([-T,T]/sqrt(2), [-T,T]/sqrt(2), ":k")
	PyPlot.plot([-T,T]/sqrt(2), [T,-T]/sqrt(2), ":k")

	for i in 1:n
		PyPlot.plot((1:T).*cos.(ths[i,:]), (1:T).*sin.(ths[i,:]))
	end

	r = vec(sum(exp.(im*ths),dims=1)/n)
	PyPlot.plot((1:T).*real.(r), (1:T).*imag.(r), "--k", linewidth=2)
end


function plot_circ_log(ths::Array{Float64,2})
	n,T = size(ths)

	figure("log",(8,8))
	ts = LinRange(0,2pi,200)
	t = 1
	while t < T
		PyPlot.plot(log(t+1)*cos.(ts), log(t+1)*sin.(ts), ":k")
		t *= 10
	end
	PyPlot.plot([0,0], [-log(T+1),log(T+1)], ":k")
	PyPlot.plot([-log(T+1),log(T+1)], [0,0], ":k")
	PyPlot.plot([-log(T+1),log(T+1)]/sqrt(2), [-log(T+1),log(T+1)]/sqrt(2), ":k")
	PyPlot.plot([-log(T+1),log(T+1)]/sqrt(2), [log(T+1),-log(T+1)]/sqrt(2), ":k")

	for i in 1:n
		PyPlot.plot(log.((1:T).+1).*cos.(ths[i,:]), log.((1:T).+1).*sin.(ths[i,:]))
	end

	r = vec(sum(exp.(im*ths),dims=1)/n)
	PyPlot.plot(log.((1:T).+1).*real.(r), log.((1:T).+1).*imag.(r), "--k", linewidth=2)
end

       
