using Random, Dates, ROC

include("hyper_inf.jl")
include("hyper_kuramoto.jl")
include("gen_rand_hyperg.jl")
include("tools_hyper.jl")

n = 7
iter = 10

ooi = [2,3]
T = 100
amp = .1
ξ0s = [.0005,.001,.002,.005,.01,.02,.05]

p1 = .05
p2 = .4

ξ = 1e-10

a1s = zeros(Float64,iter,0)
a2s = zeros(Float64,iter,0)
a3s = zeros(Float64,iter,0)

figure("ntw")

for ξ0 in ξ0s
@info "ξ0 = $ξ0"
a1 = Float64[]
a2 = Float64[]
a3 = Float64[]
for i in 1:iter
	figure("ntw")
	clf()
	A2,A3 = gen_hyper_er(n,p1,p2,true)
	adj = cat_As(A2,A3)

	X = amp*(rand(n,T) .- .5)
	Y = f_kuramoto_3rd(X,A2,A3,zeros(n),π/4,π/4) + ξ0*randn(size(X))

	xxx = hyper_inf(X,Y,ooi,4,-1e-4)

	A2h = inferred_adj_2nd(xxx[1][2],n)[2]
	A3h = inferred_adj_3rd(xxx[1][3],n)[2]
#	adh = get_adj_3rd(A2h,A3h)[1]
	adh = cat_As(A2h,A3h)

	r1 = roc(abs.(adh) + ξ*rand(Float64,size(adh)), adj)
	r2 = roc(abs.(A2h) + ξ*rand(n,n),A2)
	r3 = roc(abs.(A3h) + ξ*rand(n,n,n),A3)

	push!(a1,AUC(r1))
	push!(a2,AUC(r2))
	push!(a3,AUC(r3))
	
	figure(111)
	subplot(1,3,1)
	PyPlot.plot(r1.FPR,r1.TPR)
	subplot(1,3,2)
	PyPlot.plot(r2.FPR,r2.TPR)
	subplot(1,3,3)
	PyPlot.plot(r3.FPR,r3.TPR)
end
global a1s = [a1s a1]
global a2s = [a2s a2]
global a3s = [a3s a3]
end

figure()
for i in 1:length(ξ0s)
	PyPlot.plot([ξ0s[i],ξ0s[i]],[quantile(a1s[:,i],0.),quantile(a1s[:,i],1.)],"x",color="C0")
	PyPlot.plot([ξ0s[i],ξ0s[i]],[quantile(a1s[:,i],.25),quantile(a1s[:,i],.75)],color="C0")
	PyPlot.plot(ξ0s[i],median(a1s[:,i]),"o",color="C0")
	PyPlot.plot([ξ0s[i],ξ0s[i]],[quantile(a2s[:,i],0.),quantile(a2s[:,i],1.)],"x",color="C1")
	PyPlot.plot([ξ0s[i],ξ0s[i]],[quantile(a2s[:,i],.25),quantile(a2s[:,i],.75)],color="C1")
	PyPlot.plot(ξ0s[i],median(a2s[:,i]),"o",color="C1")
	PyPlot.plot([ξ0s[i],ξ0s[i]],[quantile(a3s[:,i],0.),quantile(a3s[:,i],1.)],"x",color="C2")
	PyPlot.plot([ξ0s[i],ξ0s[i]],[quantile(a3s[:,i],.25),quantile(a3s[:,i],.75)],color="C2")
	PyPlot.plot(ξ0s[i],median(a3s[:,i]),"o",color="C2")
end
xlabel("noise amplitude")
ylabel("AUCs")

figure()
PyPlot.plot(vec(ones(iter)*ξ0s'),vec(a1s),"xk")
xlabel("noise amplitude")
ylabel("AUCs")


