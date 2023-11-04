using Random, Dates, ROC

include("hyper_inf.jl")
include("hyper_kuramoto.jl")
include("gen_rand_hyperg.jl")
include("tools_hyper.jl")

ns = [5,10,15,20]
ns = [5,6]
iter = 10

ooi = [2,3]
T = 100
amplitude = .1
ξ0 = .0005

p1 = .1
p2 = .4

ξ = 1e-10

tt = zeros(Int64,iter,0)

figure("ntw")

for n in ns
	@info "n = $n"
	t = Int64[]
	for i in 1:iter	
		figure("ntw")
		clf()
		A2,A3 = gen_hyper_er(n,p1,p2,true)
		adj = cat_As(A2,A3)
	
		X = amplitude*(rand(n,T) .- .5)
		Y = f_kuramoto_3rd(X,A2,A3,zeros(n),π/4,π/4) + ξ0*randn(size(X))
	
		t0 = CPUtime_us()
		xxx = hyper_inf(X,Y,ooi,4,1e-1)
		t1 = CPUtime_us()
		push!(t,Int64(t1-t0)) # nanoseconds
		push!(c,Int64(c1-c0)) # microseconds
	
		A2h = inferred_adj_2nd(xxx[1][2],n)[2]
		A3h = inferred_adj_3rd(xxx[1][3],n)[2]
	#	adh = get_adj_3rd(A2h,A3h)[1]
		adh = cat_As(A2h,A3h)
	
		r1 = roc(abs.(adh) + ξ*rand(Float64,size(adh)), adj)
		r2 = roc(abs.(A2h) + ξ*rand(n,n),A2)
		r3 = roc(abs.(A3h) + ξ*rand(n,n,n),A3)
	
		figure(111)
		subplot(1,3,1)
		PyPlot.plot(r1.FPR,r1.TPR)
		subplot(1,3,2)
		PyPlot.plot(r2.FPR,r2.TPR)
		subplot(1,3,3)
		PyPlot.plot(r3.FPR,r3.TPR)
	end
	global tt = [tt t]
end

figure()
PyPlot.plot(ns,vec(mean(cc,dims=1)).*1e-6,"o")
xlabel("n")
ylabel("computation time [s]")



