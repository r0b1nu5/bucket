using PyPlot, DelimitedFiles

include("scripts.jl")

nx = 4
n0 = 0
#n0 = 4
emi = 0.
ema = 1.
#ema = .4
ne = 40
epss = Array(LinRange(emi,ema,ne))
n_run = 100
d = 1.
#d = 0.
sig = .2
n1 = 1000
n2 = 1001
n = n1 + n2

figure(111)

for i in n0+1:n0+nx
	x0 = vec(readdlm("data/x$i.csv",','))
	er = Array{Float64,2}(undef,n_run,0)
	ef = Array{Float64,2}(undef,1,0)
	em = Array{Float64,1}()
	os = Array{Float64,1}()
	Cns = Array{Float64,1}()
	Cps = Array{Float64,1}()
	C0s = Array{Float64,1}()
	for eps in epss
		e = Array{Float64,1}()
		for j in 1:n_run
			push!(e,readdlm("data/eff_x$(i)_rand$(j)_$eps.csv",',')[1])
		end
		er = [er e]
		ef = [ef readdlm("data/eff_x$(i)_fiedler_$eps.csv",',')[1]]
		push!(em,readdlm("data/eff_x$(i)_mini_$eps.csv",',')[1])
		push!(os,readdlm("data/o_x$(i)_$eps.csv",',')[1])

		A = Float64.((0 .< abs.(repeat(x0,1,n) - repeat(x0',n,1)) .< eps))
		L = diagm(0 => vec(sum(A,dims=1))) - A
		cs = clusterings(L,x0)

		push!(Cns,cs[1])
		push!(Cps,cs[2])
		push!(C0s,cs[3])
	end
	
	ee = eps_connect(x0,epss)
	xima = max(maximum(abs.(er)),maximum(abs.(ef)),maximum(abs.(em)))
	
	subplot(3,nx,i-n0)
	PyPlot.plot([ee,ee],[0,1.1*xima],"--k")
	plot_mean(abs.(er),epss)
	plot_fiedler(abs.(ef),epss)
	plot_mini(abs.(em),epss)
	xlabel("ε")
	ylabel("ξ")

	subplot(3,nx,nx+i-n0)
	PyPlot.plot(epss,os)
	xlabel("ε")
	ylabel("initial outcome")

	subplot(3,nx,2*nx+i-n0)
	PyPlot.semilogy(epss,Cns,label="C-")
	PyPlot.semilogy(epss,Cps,label="C+")
	PyPlot.semilogy(epss,C0s,label="C0")
	xlabel("ε")
	ylabel("clustering")
	legend()
end





			



