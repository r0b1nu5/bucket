using Random

include("hyper_inf.jl")
include("hyper_kuramoto.jl")
include("hyper_ktanh.jl")
include("gen_rand_hyperg.jl")
include("tools_hyper.jl")

include("../../ARNI/reconstruct.jl")

# Generating the hypergraph.
n = 7
 #=
ntw = "wheel"
p1 = 0.
p2 = 0.05
p3 = .3
# =# 
# #=
ntw = "er"
p1 = 0.
p2 = .99
# =#

amplitude = .2


if ntw == "wheel"
	A2,A3 = gen_rand_hyperwheel(n,p1,p2,p3,true)
	A4 = zeros(n,n,n,n)
elseif ntw == "er"
	A2,A3 = gen_hyper_er(n,p1,p2)
	A4 = zeros(n,n,n,n)
end

figure("ROCs")
subplot(2,3,1)
title("us")
xlabel("FPR")
ylabel("TPR")
subplot(2,3,2)
title("ARNI (polynomial)")
xlabel("FPR")
ylabel("TPR")
subplot(2,3,3)
title("ARNI (polynomial diff)")
xlabel("FPR")
ylabel("TPR")
subplot(2,3,4)
title("ARNI (fourier)")
xlabel("FPR")
ylabel("TPR")
subplot(2,3,5)
title("ARNI (fourier diff)")
xlabel("FPR")
ylabel("TPR")
subplot(2,3,6)
title("ARNI (poly) - us")
xlabel("FPR")
ylabel("TPR (ARNI-us)")

cmapme = get_cmap("RdPu")
cmaparni = get_cmap("GnBu")
cmapdiff = get_cmap("Greys")
cmapme = get_cmap("plasma")
cmaparni = get_cmap("viridis")
c0 = 0.
c1 = 1.

adj = get_adj_3rd(A2,A3)

X = amplitude*(rand(n,400) .- .5)
Y = f_kuramoto_3rd(X,A2,A3,zeros(n),π/4,π/4)

writedlm("data/test-arni-Xs.csv",X,',')
writedlm("data/test-arni-Ys.csv",Y,',')
writedlm("data/test-arni-adj.csv",A2,',')

sen2 = Float64[]
spe2 = Float64[]
sen3 = Float64[]
spe3 = Float64[]
sen4 = Float64[]
spe4 = Float64[]

fprarni1 = zeros(n^2+1,0)
tprarni1 = zeros(n^2+1,0)
aucarni1 = Float64[]
fprarni2 = zeros(n^2+1,0)
tprarni2 = zeros(n^2+1,0)
aucarni2 = Float64[]
fprarni3 = zeros(n^2+1,0)
tprarni3 = zeros(n^2+1,0)
aucarni3 = Float64[]
fprarni4 = zeros(n^2+1,0)
tprarni4 = zeros(n^2+1,0)
aucarni4 = Float64[]
fprme = zeros(n^2+1,0)
tprme = zeros(n^2+1,0)
aucme = Float64[]

# Compute the sensitivity and specificity of the inference for various lengths of time series.
iters = 10:5:200
iters = 10:5:80
ooi = [2,]
c = 0
for iter in iters
	global c += 1
	@info "Run $c/$(length(iters))"

	xxx = hyper_inf(X[:,1:iter],Y[:,1:iter],ooi,4,-1e-4)
	adjme = inferred_adj_2nd(xxx[1][2],n)[1]
	
	adjarni1 = zeros(n,n)
	adjarni2 = zeros(n,n)
	adjarni3 = zeros(n,n)
	adjarni4 = zeros(n,n)
	for i in 1:n
		w = reconstruct(X[:,1:iter],Y[:,1:iter],i,A2,"polynomial")
		adjarni1[i,:] = w[1]
		w = reconstruct(X[:,1:iter],Y[:,1:iter],i,A2,"polynomial_diff")
		adjarni2[i,:] = w[1]
		w = reconstruct(X[:,1:iter],Y[:,1:iter],i,A2,"fourier")
		adjarni3[i,:] = w[1]
		w = reconstruct(X[:,1:iter],Y[:,1:iter],i,A2,"fourier_diff")
		adjarni4[i,:] = w[1]
	end
		
	ξ = 1e-10
	rocme = roc(adjme .+ ξ*rand(n,n),A2) # Adding rand(n,n) makes sure that FPR and TPR have the same dim as fprme and tprme.
	x0 = rocme.FPR
	y0 = rocme.TPR
	global fprme = [fprme rocme.FPR]
	global tprme = [tprme rocme.TPR]
	push!(aucme,AUC(rocme))
	
	rocarni = roc(adjarni1 .+ ξ*rand(n,n),A2)
	x1 = rocarni.FPR
	y1 = rocarni.TPR
	global fprarni1 = [fprarni1 rocarni.FPR]
	global tprarni1 = [tprarni1 rocarni.TPR]
	push!(aucarni1,AUC(rocarni))
	
	rocarni = roc(adjarni2 .+ ξ*rand(n,n),A2)
	x2 = rocarni.FPR
	y2 = rocarni.TPR
	global fprarni2 = [fprarni2 rocarni.FPR]
	global tprarni2 = [tprarni2 rocarni.TPR]
	push!(aucarni2,AUC(rocarni))
	
	rocarni = roc(adjarni3 .+ ξ*rand(n,n),A2)
	x3 = rocarni.FPR
	y3 = rocarni.TPR
	global fprarni3 = [fprarni3 rocarni.FPR]
	global tprarni3 = [tprarni3 rocarni.TPR]
	push!(aucarni3,AUC(rocarni))
	
	rocarni = roc(adjarni4 .+ ξ*rand(n,n),A2)
	x4 = rocarni.FPR
	y4 = rocarni.TPR
	global fprarni4 = [fprarni4 rocarni.FPR]
	global tprarni4 = [tprarni4 rocarni.TPR]
	push!(aucarni4,AUC(rocarni))
	
	Z = sortslices([x1 y1 -ones(length(y1));x0 -ones(length(y0)) y0],dims=1)
	x = Z[:,1]
	z1 = Z[:,2]
	temp = Int64[]
	for i in 1:length(z1)
		if z1[i] < 0
			push!(temp,i)
		else
			z1[temp] .= z1[i]
			temp = Int64[]
		end
	end
	z1[temp] .= 1.
	z0 = Z[:,3]
	temp = Int64[]
	for i in 1:length(z0)
		if z0[i] < 0
			push!(temp,i)
		else
			z0[temp] .= z0[i]
			temp = Int64[]
		end
	end
	z0[temp] .= 1.

	figure("ROCs")
	subplot(2,3,1)
	PyPlot.plot(rocme.FPR,rocme.TPR,color=cmapme(c0 + c1*iter/maximum(iters)))
	subplot(2,3,2)
	PyPlot.plot(fprarni1[:,end],tprarni1[:,end],color=cmaparni(c0 + c1*iter/maximum(iters)))
	subplot(2,3,3)
	PyPlot.plot(fprarni2[:,end],tprarni2[:,end],color=cmaparni(c0 + c1*iter/maximum(iters)))
	subplot(2,3,4)
	PyPlot.plot(fprarni3[:,end],tprarni3[:,end],color=cmaparni(c0 + c1*iter/maximum(iters)))
	subplot(2,3,5)
	PyPlot.plot(fprarni4[:,end],tprarni4[:,end],color=cmaparni(c0 + c1*iter/maximum(iters)))
	subplot(2,3,6)
	PyPlot.plot(x,z0-z1,color=cmapdiff(1-.9*(c0 + c1*iter/maximum(iters))))

	yyy = check_inference_bool(A2,A3,A4,xxx[1])
	push!(sen2,yyy[1][1])
	push!(spe2,yyy[1][2])
	push!(sen3,yyy[2][1])
	push!(spe3,yyy[2][2])
	push!(sen4,yyy[3][1])
	push!(spe4,yyy[3][2])
end


