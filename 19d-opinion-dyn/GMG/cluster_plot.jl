using PyPlot, DelimitedFiles

include("scripts.jl")

# Generates plots, based on the data obtained with "cluster_run.jl". 
# (1) effort vs. eps
# (2) initial outcome vs. eps
# (3) clusterings vs. eps

nx = 3
n0 = 20
emi = .05
ema = .6
ne = 30
epss = Array(LinRange(emi,ema,ne))
n_run = 20
d = 1.
sig = .2
n1 = 1000
n2 = 1001
n = n1 + n2
n_modes = [2,3,4]

figure(111)

for i in n0+1:n0+nx
	x0 = sort(vec(readdlm("data/x$i.csv",',')))
	er = Array{Float64,2}(undef,n_run,0)
	ef = Array{Float64,2}(undef,length(n_modes),0)
	em = Array{Float64,1}()
	ec = Array{Float64,1}()
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
		ef = [ef vec(readdlm("data/eff_x$(i)_fiedler_$eps.csv",','))]
		push!(em,readdlm("data/eff_x$(i)_mini_$eps.csv",',')[1])
		push!(ec,readdlm("data/eff_x$(i)_cent_$eps.csv",',')[1])
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
	m = size(ef)[1]

	subplot(3,nx,i-n0)
	PyPlot.plot([ee,ee],[0,1.1*xima],"--k")
	plot_mean(abs.(er),epss)
	plot_fiedler(abs.(ef),epss,["C3","C4","C5","C6","C7","C8","C9"],Array(1:m))
	plot_mini(abs.(em),epss)
	plot_cent(abs.(ec),epss)
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





			



