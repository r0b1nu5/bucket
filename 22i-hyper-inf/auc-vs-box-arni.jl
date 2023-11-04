using Random, Dates, ROC

include("hyper_inf.jl")
include("hyper_kuramoto.jl")
include("gen_rand_hyperg.jl")
include("tools_hyper.jl")

include("arni-reconstruct.jl")
include("arni-reconstruct-3rd.jl")

n = 7
iter = 10

ooi = [2,3]
T = 100
amps = [.01,.02,.05,.1,.2,.5,1.,2.,5.]
ξ0 = .0005

p1 = .05
p2 = .4

ξ = 1e-10

a1s = zeros(Float64,iter,0)
a2s = zeros(Float64,iter,0)
a3s = zeros(Float64,iter,0)

figure("ntw")

for amp in amps
@info "amp = $amp"
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
	
	adjarni = zeros(n,Int64(n*(n-1)/2))
	for i in 1:n
		w = reconstruct_3rd(X,Y,i,adj,1e-6,"power_series")
		adjarni[i,:] = w[1]
	end

	A2h,A3h = adj2As(adjarni)
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
for i in 1:length(amps)
	PyPlot.plot([amps[i],amps[i]],[quantile(a1s[:,i],0.),quantile(a1s[:,i],1.)],"x",color="C0")
	PyPlot.plot([amps[i],amps[i]],[quantile(a1s[:,i],.25),quantile(a1s[:,i],.75)],color="C0")
	PyPlot.plot(amps[i],median(a1s[:,i]),"o",color="C0")
	PyPlot.plot([amps[i],amps[i]],[quantile(a2s[:,i],0.),quantile(a2s[:,i],1.)],"x",color="C1")
	PyPlot.plot([amps[i],amps[i]],[quantile(a2s[:,i],.25),quantile(a2s[:,i],.75)],color="C1")
	PyPlot.plot(amps[i],median(a2s[:,i]),"o",color="C1")
	PyPlot.plot([amps[i],amps[i]],[quantile(a3s[:,i],0.),quantile(a3s[:,i],1.)],"x",color="C2")
	PyPlot.plot([amps[i],amps[i]],[quantile(a3s[:,i],.25),quantile(a3s[:,i],.75)],color="C2")
	PyPlot.plot(amps[i],median(a3s[:,i]),"o",color="C2")
end
xlabel("box size")
ylabel("AUCs")



